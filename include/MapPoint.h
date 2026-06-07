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


#ifndef MAPPOINT_H
#define MAPPOINT_H

#include "KeyFrame.h"
#include "Frame.h"
#include "Map.h"
#include "Converter.h"

#include "SerializationUtils.h"

#include <opencv2/core/core.hpp>
#include <mutex>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>

namespace ORB_SLAM3
{

class KeyFrame;
class Map;
class Frame;

/**
 * @class MapPoint —— 地图点(3D 路标);完整职责见 MapPoint.cc 文件头。
 *        本头文件按"前生今世"标注各成员变量的生命周期(实例化 → 更新 → 使用 → 销毁)。
 */
class MapPoint
{

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & mnId;
        ar & mnFirstKFid;
        ar & mnFirstFrame;
        ar & nObs;
        // Variables used by the tracking
        //ar & mTrackProjX;
        //ar & mTrackProjY;
        //ar & mTrackDepth;
        //ar & mTrackDepthR;
        //ar & mTrackProjXR;
        //ar & mTrackProjYR;
        //ar & mbTrackInView;
        //ar & mbTrackInViewR;
        //ar & mnTrackScaleLevel;
        //ar & mnTrackScaleLevelR;
        //ar & mTrackViewCos;
        //ar & mTrackViewCosR;
        //ar & mnTrackReferenceForFrame;
        //ar & mnLastFrameSeen;

        // Variables used by local mapping
        //ar & mnBALocalForKF;
        //ar & mnFuseCandidateForKF;

        // Variables used by loop closing and merging
        //ar & mnLoopPointForKF;
        //ar & mnCorrectedByKF;
        //ar & mnCorrectedReference;
        //serializeMatrix(ar,mPosGBA,version);
        //ar & mnBAGlobalForKF;
        //ar & mnBALocalForMerge;
        //serializeMatrix(ar,mPosMerge,version);
        //serializeMatrix(ar,mNormalVectorMerge,version);

        // Protected variables
        ar & boost::serialization::make_array(mWorldPos.data(), mWorldPos.size());
        ar & boost::serialization::make_array(mNormalVector.data(), mNormalVector.size());
        //ar & BOOST_SERIALIZATION_NVP(mBackupObservationsId);
        //ar & mObservations;
        ar & mBackupObservationsId1;
        ar & mBackupObservationsId2;
        serializeMatrix(ar,mDescriptor,version);
        ar & mBackupRefKFId;
        //ar & mnVisible;
        //ar & mnFound;

        ar & mbBad;
        ar & mBackupReplacedId;

        ar & mfMinDistance;
        ar & mfMaxDistance;

    }


public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MapPoint();

    MapPoint(const Eigen::Vector3f &Pos, KeyFrame* pRefKF, Map* pMap);
    MapPoint(const double invDepth, cv::Point2f uv_init, KeyFrame* pRefKF, KeyFrame* pHostKF, Map* pMap);
    MapPoint(const Eigen::Vector3f &Pos,  Map* pMap, Frame* pFrame, const int &idxF);

    void SetWorldPos(const Eigen::Vector3f &Pos);
    Eigen::Vector3f GetWorldPos();

    Eigen::Vector3f GetNormal();
    void SetNormalVector(const Eigen::Vector3f& normal);

    KeyFrame* GetReferenceKeyFrame();

    std::map<KeyFrame*,std::tuple<int,int>> GetObservations();
    int Observations();

    void AddObservation(KeyFrame* pKF,int idx);
    void EraseObservation(KeyFrame* pKF);

    std::tuple<int,int> GetIndexInKeyFrame(KeyFrame* pKF);
    bool IsInKeyFrame(KeyFrame* pKF);

    void SetBadFlag();
    bool isBad();

    void Replace(MapPoint* pMP);    
    MapPoint* GetReplaced();

    void IncreaseVisible(int n=1);
    void IncreaseFound(int n=1);
    float GetFoundRatio();
    inline int GetFound(){
        return mnFound;
    }

    void ComputeDistinctiveDescriptors();

    cv::Mat GetDescriptor();

    void UpdateNormalAndDepth();

    float GetMinDistanceInvariance();
    float GetMaxDistanceInvariance();
    int PredictScale(const float &currentDist, KeyFrame*pKF);
    int PredictScale(const float &currentDist, Frame* pF);

    Map* GetMap();
    void UpdateMap(Map* pMap);

    void PrintObservations();

    void PreSave(set<KeyFrame*>& spKF,set<MapPoint*>& spMP);
    void PostLoad(map<long unsigned int, KeyFrame*>& mpKFid, map<long unsigned int, MapPoint*>& mpMPid);

public:
    // —— 身份与计数(构造时定,基本不变)——
    long unsigned int mnId;            ///< 全局唯一 id;前生:构造函数里 mnId=nNextId++(mMutexPointCreation 保护)
    static long unsigned int nNextId;  ///< 全局 id 发号器(所有 MapPoint 共享,单调递增)
    long int mnFirstKFid;              ///< 创建本点的关键帧 id(构造时定,溯源用)
    long int mnFirstFrame;             ///< 创建本点的帧 id
    int nObs;                          ///< 观测计数;AddObservation +1/+2、EraseObservation 减;今世:≤2 触发 SetBadFlag

    // —— 跟踪期临时量 ——
    // ★前生今世:每帧由 Frame::isInFrustum 投影本点时写入 → SearchByProjection 据此在"投影位置+预测层"找匹配
    //   → 下一帧再次 isInFrustum 时覆盖。纯每帧 scratch,不持久化(serialize 里全部注释掉)。
    float mTrackProjX;
    float mTrackProjY;
    float mTrackDepth;
    float mTrackDepthR;
    float mTrackProjXR;
    float mTrackProjYR;
    bool mbTrackInView, mbTrackInViewR;
    int mnTrackScaleLevel, mnTrackScaleLevelR;
    float mTrackViewCos, mTrackViewCosR;
    long unsigned int mnTrackReferenceForFrame;
    long unsigned int mnLastFrameSeen;

    // —— 防重复处理标记(★前生今世:被某线程以某 KF 的 id 标记"已处理过本点",避免同一轮重复入队/重复融合)——
    //   写点 = 各线程处理本点时填入当前 KF 的 mnId;读点 = 处理前比对"是否等于当前 KF id"决定跳过。无显式销毁(被新值覆盖)。
    // Variables used by local mapping
    long unsigned int mnBALocalForKF;       ///< 已被某 KF 的局部 BA 纳入
    long unsigned int mnFuseCandidateForKF; ///< 已被某 KF 选为融合候选(SearchInNeighbors)

    // Variables used by loop closing
    long unsigned int mnLoopPointForKF;
    long unsigned int mnCorrectedByKF;
    long unsigned int mnCorrectedReference;    
    Eigen::Vector3f mPosGBA;
    long unsigned int mnBAGlobalForKF;
    long unsigned int mnBALocalForMerge;

    // Variable used by merging
    Eigen::Vector3f mPosMerge;
    Eigen::Vector3f mNormalVectorMerge;


    // Fopr inverse depth optimization
    double mInvDepth;
    double mInitU;
    double mInitV;
    KeyFrame* mpHostKF;

    static std::mutex mGlobalMutex;

    unsigned int mnOriginMapId;

protected:    

     // ===== 核心状态(持久化)=====
     // 世界系 3D 坐标。前生:构造经 SetWorldPos 设入;今世:三角化新建 / BA 优化 / Map::ApplyScaledRotation 缩放时更新
     //   (GetWorldPos 读,重投影/SearchByProjection 到处用);受 mMutexPos 保护。
     Eigen::Vector3f mWorldPos;

     // 哪些 KF 在第几个特征点观测到本点(双目存左右索引)。前生:AddObservation 加;今世:UpdateNormalAndDepth /
     //   ComputeDistinctiveDescriptors / KeyFrame::UpdateConnections 遍历它;销毁:EraseObservation 删、SetBadFlag/Replace 清空。mMutexFeatures 保护。
     std::map<KeyFrame*,std::tuple<int,int> > mObservations;
     // For save relation without pointer, this is necessary for save/load function
     std::map<long unsigned int, int> mBackupObservationsId1;
     std::map<long unsigned int, int> mBackupObservationsId2;

     // 平均观测方向。前生/今世:UpdateNormalAndDepth 据各观测 KF 视线均值算出;用途:isInFrustum 视角检查。
     Eigen::Vector3f mNormalVector;

     // 代表性描述子。今世:ComputeDistinctiveDescriptors 取所有观测中位汉明距离最小者;用途:一切匹配/回环检索。
     cv::Mat mDescriptor;

     // 参考关键帧。前生:构造设为创建本点的 KF;今世:EraseObservation 若删的是它则迁到首个剩余观测 KF;UpdateNormalAndDepth 用它定尺度。
     KeyFrame* mpRefKF;
     long unsigned int mBackupRefKFId;

     // 跟踪统计。前生:构造为 1;今世:IncreaseVisible(进入视野次数)/IncreaseFound(实际匹配上次数)累加;
     //   用途:MapPointCulling 据 found/visible 比率 < 25% 剔除劣质点。
     int mnVisible;
     int mnFound;

     // 坏点标记 + 替换指针。前生:构造 mbBad=false / mpReplaced=null;今世:SetBadFlag 置 mbBad、Replace 置 mpReplaced(融合);
     //   读:isBad()/GetReplaced() 到处检查并跳过坏点。★软删除:不真正 delete 内存,只置 mbBad。
     bool mbBad;
     MapPoint* mpReplaced;
     // For save relation without pointer, this is necessary for save/load function
     long long int mBackupReplacedId;

     // 尺度不变距离范围。前生/今世:UpdateNormalAndDepth 据参考帧观测距离 × 金字塔尺度算出;
     //   用途:PredictScale / GetMin/MaxDistanceInvariance 在匹配时预测该点应在哪一金字塔层搜索。
     float mfMinDistance;
     float mfMaxDistance;

     Map* mpMap;

     // Mutex
     std::mutex mMutexPos;
     std::mutex mMutexFeatures;
     std::mutex mMutexMap;

};

} //namespace ORB_SLAM

#endif // MAPPOINT_H
