/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2021 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/


#include "Sim3Solver.h"

#include <vector>
#include <cmath>
#include <opencv2/core/core.hpp>

#include "KeyFrame.h"
#include "ORBmatcher.h"

#include "Thirdparty/DBoW2/DUtils/Random.h"

namespace ORB_SLAM3
{


// ============================================================================
// Sim3Solver —— Horn 闭式法求两关键帧间相似变换 Sim(3)=(尺度 s, 旋转 R, 平移 t),RANSAC 鲁棒化
// 用途:回环/地图合并时,由两帧匹配的 3D 地图点对解出 KF2→KF1 的 (sR,t)。★单目需估尺度 s(mbFixScale=false),
//   双目/RGBD 尺度已知则固定 s=1(mbFixScale=true)。解出的 Sim3 作回环约束喂给本质图优化 OptimizeSim3/EssentialGraph。
// 流程:构造时把匹配点对转到各自相机系 3D(mvX3Dc1/2)→ iterate 跑 RANSAC:每轮取 3 对点 ComputeSim3
//   (Horn 闭式)→ CheckInliers(双向重投影卡方)→ 留最优。Horn 1987: Closed-form solution of absolute orientation。
// ============================================================================
/**
 * @brief 构造:从两 KF 的匹配地图点对建立 Sim3 求解所需数据(各自相机系 3D + 图像投影 + 误差阈值)
 * @details 遍历匹配 vpMatched12 过滤坏点;把 KF1/KF2 地图点世界坐标转到各自相机系存 mvX3Dc1/mvX3Dc2;
 *   记录卡方阈值(9.21×σ²,自由度 2 的 99%)。bFixScale 决定是否固定尺度。
 */
Sim3Solver::Sim3Solver(KeyFrame *pKF1, KeyFrame *pKF2, const vector<MapPoint *> &vpMatched12, const bool bFixScale,
                       vector<KeyFrame*> vpKeyFrameMatchedMP):
    mnIterations(0), mnBestInliers(0), mbFixScale(bFixScale),
    pCamera1(pKF1->mpCamera), pCamera2(pKF2->mpCamera)
{
    bool bDifferentKFs = false;
    if(vpKeyFrameMatchedMP.empty())
    {
        bDifferentKFs = true;
        vpKeyFrameMatchedMP = vector<KeyFrame*>(vpMatched12.size(), pKF2);
    }

    mpKF1 = pKF1;
    mpKF2 = pKF2;

    vector<MapPoint*> vpKeyFrameMP1 = pKF1->GetMapPointMatches();

    mN1 = vpMatched12.size();

    mvpMapPoints1.reserve(mN1);
    mvpMapPoints2.reserve(mN1);
    mvpMatches12 = vpMatched12;
    mvnIndices1.reserve(mN1);
    mvX3Dc1.reserve(mN1);
    mvX3Dc2.reserve(mN1);

    Eigen::Matrix3f Rcw1 = pKF1->GetRotation();
    Eigen::Vector3f tcw1 = pKF1->GetTranslation();
    Eigen::Matrix3f Rcw2 = pKF2->GetRotation();
    Eigen::Vector3f tcw2 = pKF2->GetTranslation();

    mvAllIndices.reserve(mN1);

    size_t idx=0;

    KeyFrame* pKFm = pKF2; //Default variable
    for(int i1=0; i1<mN1; i1++)
    {
        if(vpMatched12[i1])
        {
            MapPoint* pMP1 = vpKeyFrameMP1[i1];
            MapPoint* pMP2 = vpMatched12[i1];

            if(!pMP1)
                continue;

            if(pMP1->isBad() || pMP2->isBad())
                continue;

            if(bDifferentKFs)
                pKFm = vpKeyFrameMatchedMP[i1];

            int indexKF1 = get<0>(pMP1->GetIndexInKeyFrame(pKF1));
            int indexKF2 = get<0>(pMP2->GetIndexInKeyFrame(pKFm));

            if(indexKF1<0 || indexKF2<0)
                continue;

            const cv::KeyPoint &kp1 = pKF1->mvKeysUn[indexKF1];
            const cv::KeyPoint &kp2 = pKFm->mvKeysUn[indexKF2];

            const float sigmaSquare1 = pKF1->mvLevelSigma2[kp1.octave];
            const float sigmaSquare2 = pKFm->mvLevelSigma2[kp2.octave];

            mvnMaxError1.push_back(9.210*sigmaSquare1);
            mvnMaxError2.push_back(9.210*sigmaSquare2);

            mvpMapPoints1.push_back(pMP1);
            mvpMapPoints2.push_back(pMP2);
            mvnIndices1.push_back(i1);

            Eigen::Vector3f X3D1w = pMP1->GetWorldPos();
            mvX3Dc1.push_back(Rcw1*X3D1w+tcw1);

            Eigen::Vector3f X3D2w = pMP2->GetWorldPos();
            mvX3Dc2.push_back(Rcw2*X3D2w+tcw2);

            mvAllIndices.push_back(idx);
            idx++;
        }
    }

    FromCameraToImage(mvX3Dc1,mvP1im1,pCamera1);
    FromCameraToImage(mvX3Dc2,mvP2im2,pCamera2);

    SetRansacParameters();
}

/// @brief 据期望成功概率/最小内点数/对应点总数,推算 RANSAC 最大迭代次数(每次抽样 3 对点)
void Sim3Solver::SetRansacParameters(double probability, int minInliers, int maxIterations)
{
    mRansacProb = probability;
    mRansacMinInliers = minInliers;
    mRansacMaxIts = maxIterations;    

    N = mvpMapPoints1.size(); // number of correspondences

    mvbInliersi.resize(N);

    // Adjust Parameters according to number of correspondences
    float epsilon = (float)mRansacMinInliers/N;

    // Set RANSAC iterations according to probability, epsilon, and max iterations
    int nIterations;

    if(mRansacMinInliers==N)
        nIterations=1;
    else
        nIterations = ceil(log(1-mRansacProb)/log(1-pow(epsilon,3)));

    mRansacMaxIts = max(1,min(nIterations,mRansacMaxIts));

    mnIterations = 0;
}

/**
 * @brief RANSAC 主循环:反复抽 3 对点算 Sim3 + 统计内点,返回内点最多的 Sim3(T12)
 * @details 每轮随机取 3 对对应点 → ComputeSim3(Horn 闭式)→ CheckInliers;内点超过阈值即提前返回成功。
 *   (另有带 bConverge 的重载,语义相同,多回报是否收敛。)
 */
Eigen::Matrix4f Sim3Solver::iterate(int nIterations, bool &bNoMore, vector<bool> &vbInliers, int &nInliers)
{
    bNoMore = false;
    vbInliers = vector<bool>(mN1,false);
    nInliers=0;

    if(N<mRansacMinInliers)
    {
        bNoMore = true;
        return Eigen::Matrix4f::Identity();
    }

    vector<size_t> vAvailableIndices;

    Eigen::Matrix3f P3Dc1i;
    Eigen::Matrix3f P3Dc2i;

    int nCurrentIterations = 0;
    while(mnIterations<mRansacMaxIts && nCurrentIterations<nIterations)
    {
        nCurrentIterations++;
        mnIterations++;

        vAvailableIndices = mvAllIndices;

        // Get min set of points
        for(short i = 0; i < 3; ++i)
        {
            int randi = DUtils::Random::RandomInt(0, vAvailableIndices.size()-1);

            int idx = vAvailableIndices[randi];

            P3Dc1i.col(i) = mvX3Dc1[idx];
            P3Dc2i.col(i) = mvX3Dc2[idx];

            vAvailableIndices[randi] = vAvailableIndices.back();
            vAvailableIndices.pop_back();
        }

        ComputeSim3(P3Dc1i,P3Dc2i);

        CheckInliers();

        if(mnInliersi>=mnBestInliers)
        {
            mvbBestInliers = mvbInliersi;
            mnBestInliers = mnInliersi;
            mBestT12 = mT12i;
            mBestRotation = mR12i;
            mBestTranslation = mt12i;
            mBestScale = ms12i;

            if(mnInliersi>mRansacMinInliers)
            {
                nInliers = mnInliersi;
                for(int i=0; i<N; i++)
                    if(mvbInliersi[i])
                        vbInliers[mvnIndices1[i]] = true;
                return mBestT12;
            }
        }
    }

    if(mnIterations>=mRansacMaxIts)
        bNoMore=true;

    return Eigen::Matrix4f::Identity();
}

Eigen::Matrix4f Sim3Solver::iterate(int nIterations, bool &bNoMore, vector<bool> &vbInliers, int &nInliers, bool &bConverge)
{
    bNoMore = false;
    bConverge = false;
    vbInliers = vector<bool>(mN1,false);
    nInliers=0;

    if(N<mRansacMinInliers)
    {
        bNoMore = true;
        return Eigen::Matrix4f::Identity();
    }

    vector<size_t> vAvailableIndices;

    Eigen::Matrix3f P3Dc1i;
    Eigen::Matrix3f P3Dc2i;

    int nCurrentIterations = 0;

    Eigen::Matrix4f bestSim3;

    while(mnIterations<mRansacMaxIts && nCurrentIterations<nIterations)
    {
        nCurrentIterations++;
        mnIterations++;

        vAvailableIndices = mvAllIndices;

        // Get min set of points
        for(short i = 0; i < 3; ++i)
        {
            int randi = DUtils::Random::RandomInt(0, vAvailableIndices.size()-1);

            int idx = vAvailableIndices[randi];

            P3Dc1i.col(i) = mvX3Dc1[idx];
            P3Dc2i.col(i) = mvX3Dc2[idx];

            vAvailableIndices[randi] = vAvailableIndices.back();
            vAvailableIndices.pop_back();
        }

        ComputeSim3(P3Dc1i,P3Dc2i);

        CheckInliers();

        if(mnInliersi>=mnBestInliers)
        {
            mvbBestInliers = mvbInliersi;
            mnBestInliers = mnInliersi;
            mBestT12 = mT12i;
            mBestRotation = mR12i;
            mBestTranslation = mt12i;
            mBestScale = ms12i;

            if(mnInliersi>mRansacMinInliers)
            {
                nInliers = mnInliersi;
                for(int i=0; i<N; i++)
                    if(mvbInliersi[i])
                        vbInliers[mvnIndices1[i]] = true;
                bConverge = true;
                return mBestT12;
            }
            else
            {
                bestSim3 = mBestT12;
            }
        }
    }

    if(mnIterations>=mRansacMaxIts)
        bNoMore=true;

    return bestSim3;
}

Eigen::Matrix4f Sim3Solver::find(vector<bool> &vbInliers12, int &nInliers)
{
    bool bFlag;
    return iterate(mRansacMaxIts,bFlag,vbInliers12,nInliers);
}

void Sim3Solver::ComputeCentroid(Eigen::Matrix3f &P, Eigen::Matrix3f &Pr, Eigen::Vector3f &C)
{
    C = P.rowwise().sum();
    C = C / P.cols();
    for(int i=0; i<P.cols(); i++)
    Pr.col(i) = P.col(i) - C;
}


/**
 * @brief ★Horn 1987 闭式法:由 3+ 对 3D 点求相似变换 (s,R,t),使 P1 ≈ s·R·P2 + t
 * @details 8 步:① 各点集去质心(ComputeCentroid)② M = Pr2·Pr1ᵀ ③ 由 M 构造 4×4 对称阵 N
 *   ④ N 的最大特征值对应特征向量 = 最优旋转四元数 → mR12i ⑤ 旋转 set2 ⑥ 尺度 s = Σ(Pr1·P3)/Σ(P3·P3)
 *   (mbFixScale 时 s=1)⑦ 平移 t = O1 − s·R·O2 ⑧ 组装 T12 及其逆 T21。闭式解、无需迭代优化。
 */
void Sim3Solver::ComputeSim3(Eigen::Matrix3f &P1, Eigen::Matrix3f &P2)
{
    // Custom implementation of:
    // Horn 1987, Closed-form solution of absolute orientataion using unit quaternions

    // Step 1: Centroid and relative coordinates

    Eigen::Matrix3f Pr1; // Relative coordinates to centroid (set 1)
    Eigen::Matrix3f Pr2; // Relative coordinates to centroid (set 2)
    Eigen::Vector3f O1; // Centroid of P1
    Eigen::Vector3f O2; // Centroid of P2

    ComputeCentroid(P1,Pr1,O1);
    ComputeCentroid(P2,Pr2,O2);

    // Step 2: Compute M matrix

    Eigen::Matrix3f M = Pr2 * Pr1.transpose();

    // Step 3: Compute N matrix
    double N11, N12, N13, N14, N22, N23, N24, N33, N34, N44;

    Eigen::Matrix4f N;

    N11 = M(0,0)+M(1,1)+M(2,2);
    N12 = M(1,2)-M(2,1);
    N13 = M(2,0)-M(0,2);
    N14 = M(0,1)-M(1,0);
    N22 = M(0,0)-M(1,1)-M(2,2);
    N23 = M(0,1)+M(1,0);
    N24 = M(2,0)+M(0,2);
    N33 = -M(0,0)+M(1,1)-M(2,2);
    N34 = M(1,2)+M(2,1);
    N44 = -M(0,0)-M(1,1)+M(2,2);

    N << N11, N12, N13, N14,
         N12, N22, N23, N24,
         N13, N23, N33, N34,
         N14, N24, N34, N44;


    // Step 4: Eigenvector of the highest eigenvalue
    Eigen::EigenSolver<Eigen::Matrix4f> eigSolver;
    eigSolver.compute(N);

    Eigen::Vector4f eval = eigSolver.eigenvalues().real();
    Eigen::Matrix4f evec = eigSolver.eigenvectors().real(); //evec[0] is the quaternion of the desired rotation

    int maxIndex; // should be zero
    eval.maxCoeff(&maxIndex);

    Eigen::Vector3f vec = evec.block<3,1>(1,maxIndex); //extract imaginary part of the quaternion (sin*axis)

    // Rotation angle. sin is the norm of the imaginary part, cos is the real part
    double ang=atan2(vec.norm(),evec(0,maxIndex));

    vec = 2*ang*vec/vec.norm(); //Angle-axis representation. quaternion angle is the half
    mR12i = Sophus::SO3f::exp(vec).matrix();

    // Step 5: Rotate set 2
    Eigen::Matrix3f P3 = mR12i*Pr2;

    // Step 6: Scale

    if(!mbFixScale)
    {
        double cvnom = Converter::toCvMat(Pr1).dot(Converter::toCvMat(P3));
        double nom = (Pr1.array() * P3.array()).sum();
        if (abs(nom-cvnom)>1e-3)
            std::cout << "sim3 solver: " << abs(nom-cvnom) << std::endl << nom << std::endl;
        Eigen::Array<float,3,3> aux_P3;
        aux_P3 = P3.array() * P3.array();
        double den = aux_P3.sum();

        ms12i = nom/den;
    }
    else
        ms12i = 1.0f;

    // Step 7: Translation
    mt12i = O1 - ms12i * mR12i * O2;

    // Step 8: Transformation

    // Step 8.1 T12
    mT12i.setIdentity();

    Eigen::Matrix3f sR = ms12i*mR12i;
    mT12i.block<3,3>(0,0) = sR;
    mT12i.block<3,1>(0,3) = mt12i;


    // Step 8.2 T21
    mT21i.setIdentity();
    Eigen::Matrix3f sRinv = (1.0/ms12i)*mR12i.transpose();

    // sRinv.copyTo(mT21i.rowRange(0,3).colRange(0,3));
    mT21i.block<3,3>(0,0) = sRinv;

    Eigen::Vector3f tinv = -sRinv * mt12i;
    mT21i.block<3,1>(0,3) = tinv;
}


/// @brief 双向重投影检验内点:用 T12/T21 把点投到对方图像,两向重投影误差均小于卡方阈值即判内点
void Sim3Solver::CheckInliers()
{
    vector<Eigen::Vector2f> vP1im2, vP2im1;
    Project(mvX3Dc2,vP2im1,mT12i,pCamera1);
    Project(mvX3Dc1,vP1im2,mT21i,pCamera2);

    mnInliersi=0;

    for(size_t i=0; i<mvP1im1.size(); i++)
    {
        Eigen::Vector2f dist1 = mvP1im1[i] - vP2im1[i];
        Eigen::Vector2f dist2 = vP1im2[i] - mvP2im2[i];

        const float err1 = dist1.dot(dist1);
        const float err2 = dist2.dot(dist2);

        if(err1<mvnMaxError1[i] && err2<mvnMaxError2[i])
        {
            mvbInliersi[i]=true;
            mnInliersi++;
        }
        else
            mvbInliersi[i]=false;
    }
}

Eigen::Matrix4f Sim3Solver::GetEstimatedTransformation()
{
    return mBestT12;
}

Eigen::Matrix3f Sim3Solver::GetEstimatedRotation()
{
    return mBestRotation;
}

Eigen::Vector3f Sim3Solver::GetEstimatedTranslation()
{
    return mBestTranslation;
}

float Sim3Solver::GetEstimatedScale()
{
    return mBestScale;
}

void Sim3Solver::Project(const vector<Eigen::Vector3f> &vP3Dw, vector<Eigen::Vector2f> &vP2D, Eigen::Matrix4f Tcw, GeometricCamera* pCamera)
{
    Eigen::Matrix3f Rcw = Tcw.block<3,3>(0,0);
    Eigen::Vector3f tcw = Tcw.block<3,1>(0,3);

    vP2D.clear();
    vP2D.reserve(vP3Dw.size());

    for(size_t i=0, iend=vP3Dw.size(); i<iend; i++)
    {
        Eigen::Vector3f P3Dc = Rcw*vP3Dw[i]+tcw;
        Eigen::Vector2f pt2D = pCamera->project(P3Dc);
        vP2D.push_back(pt2D);
    }
}

void Sim3Solver::FromCameraToImage(const vector<Eigen::Vector3f> &vP3Dc, vector<Eigen::Vector2f> &vP2D, GeometricCamera* pCamera)
{
    vP2D.clear();
    vP2D.reserve(vP3Dc.size());

    for(size_t i=0, iend=vP3Dc.size(); i<iend; i++)
    {
        Eigen::Vector2f pt2D = pCamera->project(vP3Dc[i]);
        vP2D.push_back(pt2D);
    }
}

} //namespace ORB_SLAM
