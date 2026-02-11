# ORB-SLAM3 完整算法分析文档

> 本文档深入分析 ORB-SLAM3 的核心算法原理、数学推导及代码实现，涵盖 ORB 特征提取、IMU 预积分、Sim3 求解、因子图优化等核心模块。

## 目录

1. [算法概述与系统架构](#1-算法概述与系统架构)
2. [ORB特征提取算法](#2-orb特征提取算法)
3. [IMU预积分理论与实现](#3-imu预积分理论与实现)
4. [李群李代数与SO(3)雅可比](#4-李群李代数与so3雅可比)
5. [Sim3相似变换求解](#5-sim3相似变换求解)
6. [因子图优化框架](#6-因子图优化框架)
7. [系统初始化流程](#7-系统初始化流程)
8. [回环检测与地图融合](#8-回环检测与地图融合)
9. [核心数据结构](#9-核心数据结构)
10. [算法对比分析](#10-算法对比分析)

---

## 1. 算法概述与系统架构

### 1.1 ORB-SLAM系列演进

| 版本 | 年份 | 传感器 | 核心贡献 |
|------|------|--------|----------|
| ORB-SLAM | 2015 | 单目 | 完整SLAM流程 + DBoW2回环 |
| ORB-SLAM2 | 2017 | 单目/双目/RGB-D | 多传感器支持 |
| **ORB-SLAM3** | 2021 | +IMU | Atlas多地图 + 视觉惯性融合 |

### 1.2 系统架构

```
┌─────────────────────────────────────────────────────────────────────┐
│                         ORB-SLAM3 系统架构                           │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  ┌─────────────┐    ┌─────────────────┐    ┌───────────────────┐  │
│  │  传感器输入  │───►│   Tracking线程   │───►│  LocalMapping线程  │  │
│  │ Camera+IMU  │    │  (实时位姿估计)  │    │   (局部建图)       │  │
│  └─────────────┘    └────────┬────────┘    └─────────┬─────────┘  │
│                              │                       │             │
│                              │     ┌─────────────────┘             │
│                              │     │                               │
│                              ▼     ▼                               │
│                        ┌───────────────────┐                       │
│                        │  LoopClosing线程  │                       │
│                        │  (回环检测/融合)   │                       │
│                        └─────────┬─────────┘                       │
│                                  │                                 │
│                                  ▼                                 │
│                           ┌───────────┐                            │
│                           │   Atlas   │  ← 多地图管理器             │
│                           │ (Map集合)  │                            │
│                           └───────────┘                            │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

### 1.3 状态向量定义

对于视觉惯性SLAM，关键帧 $i$ 的状态向量为：

$$\mathbf{x}_i = [\mathbf{R}_{wb}^i, \mathbf{p}_{wb}^i, \mathbf{v}_{w}^i, \mathbf{b}_g^i, \mathbf{b}_a^i]$$

- $\mathbf{R}_{wb}^i \in SO(3)$：IMU坐标系到世界系的旋转
- $\mathbf{p}_{wb}^i \in \mathbb{R}^3$：IMU在世界系中的位置
- $\mathbf{v}_{w}^i \in \mathbb{R}^3$：世界系中的速度
- $\mathbf{b}_g^i, \mathbf{b}_a^i \in \mathbb{R}^3$：陀螺仪和加速度计偏置

---

## 2. ORB特征提取算法

### 2.1 ORB特征概述

ORB (Oriented FAST and Rotated BRIEF) 是一种快速、鲁棒的特征点检测与描述算法，由 Rublee et al. (2011) 提出。

**核心组成**：
- **检测器**：FAST角点检测 + 方向计算
- **描述子**：旋转不变的BRIEF (rBRIEF)

### 2.2 图像金字塔构建

位置：`src/ORBextractor.cc:ComputePyramid()`

构建 $L$ 层图像金字塔，尺度因子为 $s$（默认1.2）：

$$I_l = \text{resize}(I_0, \frac{1}{s^l}), \quad l = 0, 1, ..., L-1$$

每层分配的特征点数量按面积比例递减：

$$N_l = N_{total} \cdot \frac{(1-s^{-2})}{(1-s^{-2L})} \cdot s^{-2l}$$

代码实现：
```cpp
// include/ORBextractor.h:101
std::vector<int> mnFeaturesPerLevel;

// 构造函数中计算
float factor = 1.0f / scaleFactor;
float nDesiredFeaturesPerScale = nfeatures * (1 - factor) / (1 - pow(factor, nlevels));
for (int level = 0; level < nlevels; level++) {
    mnFeaturesPerLevel[level] = cvRound(nDesiredFeaturesPerScale);
    nDesiredFeaturesPerScale *= factor;
}
```

### 2.3 FAST角点检测

FAST算法检测像素 $p$ 是否为角点的准则：

在以 $p$ 为圆心、半径为3像素的 Bresenham 圆上，检查是否存在连续 $n$（默认12）个像素满足：

$$I(x) > I(p) + t \quad \text{或} \quad I(x) < I(p) - t$$

其中 $t$ 为阈值。ORB-SLAM3使用两级阈值：
- `iniThFAST = 20`：初始阈值
- `minThFAST = 7`：最小阈值（特征不足时降级使用）

### 2.4 特征点方向计算（灰度质心法）

为实现旋转不变性，使用灰度质心法计算特征点方向：

$$m_{pq} = \sum_{x,y} x^p y^q I(x,y)$$

质心坐标：
$$C = \left(\frac{m_{10}}{m_{00}}, \frac{m_{01}}{m_{00}}\right)$$

方向角：
$$\theta = \arctan2(m_{01}, m_{10})$$

代码实现（使用预计算的圆形mask `umax`）：
```cpp
// src/ORBextractor.cc 中的 IC_Angle 函数
static float IC_Angle(const Mat& image, Point2f pt, const vector<int> & u_max) {
    int m_01 = 0, m_10 = 0;
    const uchar* center = &image.at<uchar>(cvRound(pt.y), cvRound(pt.x));

    // 逐行累加
    for (int u = -HALF_PATCH_SIZE; u <= HALF_PATCH_SIZE; ++u)
        m_10 += u * center[u];

    for (int v = 1; v <= HALF_PATCH_SIZE; ++v) {
        int v_sum = 0, d = u_max[v];
        for (int u = -d; u <= d; ++u) {
            int val_plus = center[u + v*step];
            int val_minus = center[u - v*step];
            v_sum += (val_plus - val_minus);
            m_10 += u * (val_plus + val_minus);
        }
        m_01 += v * v_sum;
    }

    return fastAtan2((float)m_01, (float)m_10);
}
```

### 2.5 rBRIEF描述子

BRIEF描述子通过比较预定义点对的灰度值生成256位二进制描述子：

$$d_i = \begin{cases} 1 & \text{if } I(\mathbf{p}_i) < I(\mathbf{q}_i) \\ 0 & \text{otherwise} \end{cases}$$

为实现旋转不变性，点对根据特征点方向 $\theta$ 旋转：

$$\begin{pmatrix} p'_x \\ p'_y \end{pmatrix} = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} \begin{pmatrix} p_x \\ p_y \end{pmatrix}$$

点对模式 (`pattern`) 在代码中预定义，共256对点（512个坐标）。

### 2.6 四叉树特征点分布

位置：`src/ORBextractor.cc:DistributeOctTree()`

为保证特征点均匀分布，使用四叉树（实际是四分法）进行分割：

```
┌─────────────────────────────────┐
│  初始：整个图像作为一个节点     │
└─────────────────────────────────┘
                │
                ▼
┌───────┬───────┬───────┬───────┐
│  UL   │  UR   │  BL   │  BR   │  ← 四分
└───────┴───────┴───────┴───────┘
                │
                ▼
       递归分割，直到：
       1. 节点数 >= 目标特征数
       2. 所有节点都不可再分
                │
                ▼
       每个节点保留响应最强的特征点
```

算法伪代码：
```
while (节点数 < 目标特征数):
    for 每个可分割节点:
        将节点分成4个子节点
        移除原节点
    if 无法再分割:
        break

for 每个节点:
    保留响应值最大的特征点
```

---

## 3. IMU预积分理论与实现

### 3.1 IMU测量模型

加速度计测量：
$$\tilde{\mathbf{a}}_t = \mathbf{R}_{wb}^T(\mathbf{a}_t - \mathbf{g}) + \mathbf{b}_a + \mathbf{n}_a$$

陀螺仪测量：
$$\tilde{\boldsymbol{\omega}}_t = \boldsymbol{\omega}_t + \mathbf{b}_g + \mathbf{n}_g$$

其中：
- $\mathbf{n}_a \sim \mathcal{N}(0, \sigma_a^2 \mathbf{I})$：加速度计白噪声
- $\mathbf{n}_g \sim \mathcal{N}(0, \sigma_g^2 \mathbf{I})$：陀螺仪白噪声
- $\mathbf{b}_a, \mathbf{b}_g$：缓慢变化的偏置（随机游走模型）

偏置随机游走：
$$\dot{\mathbf{b}}_a = \mathbf{n}_{ba}, \quad \dot{\mathbf{b}}_g = \mathbf{n}_{bg}$$

### 3.2 预积分量定义

在帧 $i$ 和 $j$ 之间的预积分量定义为（相对于帧 $i$ 的IMU坐标系）：

**旋转预积分**：
$$\Delta\mathbf{R}_{ij} = \prod_{k=i}^{j-1} \text{Exp}\big((\tilde{\boldsymbol{\omega}}_k - \mathbf{b}_g)\Delta t\big)$$

**速度预积分**：
$$\Delta\mathbf{v}_{ij} = \sum_{k=i}^{j-1} \Delta\mathbf{R}_{ik} (\tilde{\mathbf{a}}_k - \mathbf{b}_a) \Delta t$$

**位置预积分**：
$$\Delta\mathbf{p}_{ij} = \sum_{k=i}^{j-1} \bigg[\Delta\mathbf{v}_{ik}\Delta t + \frac{1}{2}\Delta\mathbf{R}_{ik}(\tilde{\mathbf{a}}_k - \mathbf{b}_a)\Delta t^2\bigg]$$

### 3.3 预积分递推公式

位置：`src/ImuTypes.cc:177-235`

离散时间递推（**欧拉积分**）：

```cpp
void Preintegrated::IntegrateNewMeasurement(
    const Eigen::Vector3f &acceleration,  // 加速度测量
    const Eigen::Vector3f &angVel,         // 角速度测量
    const float &dt)                       // 时间间隔
{
    // 去除偏置
    Eigen::Vector3f acc = acceleration - b.ba;   // ã - b_a
    Eigen::Vector3f accW = angVel - b.bg;        // ω̃ - b_g

    // 先更新位置和速度（依赖未更新的旋转）
    dP = dP + dV*dt + 0.5f*dR*acc*dt*dt;  // Δp += Δv·dt + 0.5·ΔR·a·dt²
    dV = dV + dR*acc*dt;                   // Δv += ΔR·a·dt

    // 更新旋转
    IntegratedRotation dRi(angVel, b, dt);  // 计算 Exp((ω-bg)·dt)
    dR = NormalizeRotation(dR * dRi.deltaR); // ΔR *= Exp(...)

    dT += dt;
}
```

**关键细节**：
1. ORB-SLAM3使用**欧拉法**（先p,v后R），而VINS使用**中点法**
2. 旋转在位置和速度之后更新，保证当前迭代使用旧的 $\Delta\mathbf{R}$

### 3.4 协方差传播

预积分量的协方差通过线性化传播：

$$\boldsymbol{\Sigma}_{k+1} = \mathbf{A}_k \boldsymbol{\Sigma}_k \mathbf{A}_k^T + \mathbf{B}_k \boldsymbol{\Sigma}_\eta \mathbf{B}_k^T$$

状态转移矩阵 $\mathbf{A}$（9×9，对应 $[\delta\boldsymbol{\phi}, \delta\mathbf{v}, \delta\mathbf{p}]$）：

$$\mathbf{A} = \begin{bmatrix} \Delta\mathbf{R}_k^T & \mathbf{0} & \mathbf{0} \\ -\Delta\mathbf{R}\lfloor\mathbf{a}\rfloor_\times dt & \mathbf{I} & \mathbf{0} \\ -\frac{1}{2}\Delta\mathbf{R}\lfloor\mathbf{a}\rfloor_\times dt^2 & \mathbf{I}dt & \mathbf{I} \end{bmatrix}$$

噪声雅可比 $\mathbf{B}$（9×6，对应 $[\mathbf{n}_g, \mathbf{n}_a]$）：

$$\mathbf{B} = \begin{bmatrix} \mathbf{J}_r dt & \mathbf{0} \\ \mathbf{0} & \Delta\mathbf{R} dt \\ \mathbf{0} & \frac{1}{2}\Delta\mathbf{R} dt^2 \end{bmatrix}$$

代码实现：
```cpp
// src/ImuTypes.cc:203-228
Eigen::Matrix<float,3,3> Wacc = Sophus::SO3f::hat(acc);  // [a]×

A.block<3,3>(3,0) = -dR*dt*Wacc;           // ∂Δv/∂δφ
A.block<3,3>(6,0) = -0.5f*dR*dt*dt*Wacc;   // ∂Δp/∂δφ
A.block<3,3>(6,3) = Eigen::DiagonalMatrix<float,3>(dt, dt, dt);  // ∂Δp/∂Δv

B.block<3,3>(3,3) = dR*dt;                  // ∂Δv/∂n_a
B.block<3,3>(6,3) = 0.5f*dR*dt*dt;         // ∂Δp/∂n_a

C.block<9,9>(0,0) = A * C.block<9,9>(0,0) * A.transpose() + B*Nga*B.transpose();
```

### 3.5 Bias修正的雅可比矩阵

为了在bias变化时避免重新积分，计算预积分量对bias的雅可比：

$$\frac{\partial \Delta\mathbf{R}_{ij}}{\partial \mathbf{b}_g} \approx -\sum_{k=i}^{j-1} \Delta\mathbf{R}_{k+1,j}^T \mathbf{J}_r^k \Delta t$$

$$\frac{\partial \Delta\mathbf{v}_{ij}}{\partial \mathbf{b}_g} = -\sum_{k=i}^{j-1} \Delta\mathbf{R}_{ik} \lfloor\mathbf{a}_k\rfloor_\times \frac{\partial \Delta\mathbf{R}_{ik}}{\partial \mathbf{b}_g} \Delta t$$

$$\frac{\partial \Delta\mathbf{v}_{ij}}{\partial \mathbf{b}_a} = -\sum_{k=i}^{j-1} \Delta\mathbf{R}_{ik} \Delta t$$

代码实现：
```cpp
// src/ImuTypes.cc:212-216
JPa = JPa + JVa*dt - 0.5f*dR*dt*dt;
JPg = JPg + JVg*dt - 0.5f*dR*dt*dt*Wacc*JRg;
JVa = JVa - dR*dt;
JVg = JVg - dR*dt*Wacc*JRg;

// src/ImuTypes.cc:231
JRg = dRi.deltaR.transpose()*JRg - dRi.rightJ*dt;
```

### 3.6 Bias一阶修正

当bias从 $\bar{\mathbf{b}}$ 变为 $\hat{\mathbf{b}}$ 时，预积分量可通过一阶泰勒展开修正：

$$\Delta\mathbf{R}_{ij}(\hat{\mathbf{b}}_g) \approx \Delta\mathbf{R}_{ij}(\bar{\mathbf{b}}_g) \cdot \text{Exp}\left(\frac{\partial \Delta\mathbf{R}_{ij}}{\partial \mathbf{b}_g} \delta\mathbf{b}_g\right)$$

$$\Delta\mathbf{v}_{ij}(\hat{\mathbf{b}}) \approx \Delta\mathbf{v}_{ij}(\bar{\mathbf{b}}) + \frac{\partial \Delta\mathbf{v}_{ij}}{\partial \mathbf{b}_g}\delta\mathbf{b}_g + \frac{\partial \Delta\mathbf{v}_{ij}}{\partial \mathbf{b}_a}\delta\mathbf{b}_a$$

代码实现：
```cpp
// src/ImuTypes.cc:283-307
Eigen::Matrix3f Preintegrated::GetDeltaRotation(const Bias &b_) {
    Eigen::Vector3f dbg;
    dbg << b_.bwx-b.bwx, b_.bwy-b.bwy, b_.bwz-b.bwz;  // δb_g
    return NormalizeRotation(dR * Sophus::SO3f::exp(JRg * dbg).matrix());
}

Eigen::Vector3f Preintegrated::GetDeltaVelocity(const Bias &b_) {
    Eigen::Vector3f dbg, dba;
    dbg << b_.bwx-b.bwx, b_.bwy-b.bwy, b_.bwz-b.bwz;
    dba << b_.bax-b.bax, b_.bay-b.bay, b_.baz-b.baz;
    return dV + JVg * dbg + JVa * dba;
}
```

---

## 4. 李群李代数与SO(3)雅可比

### 4.1 SO(3)指数映射

从李代数 $\boldsymbol{\phi} = \theta\mathbf{a}$（$\mathbf{a}$为单位向量）到旋转矩阵的映射（Rodrigues公式）：

$$\text{Exp}(\boldsymbol{\phi}) = \mathbf{I} + \frac{\sin\theta}{\theta}\lfloor\boldsymbol{\phi}\rfloor_\times + \frac{1-\cos\theta}{\theta^2}\lfloor\boldsymbol{\phi}\rfloor_\times^2$$

当 $\theta \to 0$ 时：
$$\text{Exp}(\boldsymbol{\phi}) \approx \mathbf{I} + \lfloor\boldsymbol{\phi}\rfloor_\times$$

代码实现：
```cpp
// src/ImuTypes.cc:84-105
IntegratedRotation::IntegratedRotation(const Eigen::Vector3f &angVel,
                                        const Bias &imuBias, const float &time) {
    const float x = (angVel(0)-imuBias.bwx)*time;  // φ = (ω - b_g) * dt
    const float y = (angVel(1)-imuBias.bwy)*time;
    const float z = (angVel(2)-imuBias.bwz)*time;

    const float d2 = x*x + y*y + z*z;  // θ²
    const float d = sqrt(d2);           // θ

    Eigen::Matrix3f W = Sophus::SO3f::hat(v);  // [φ]×

    if(d < eps) {
        deltaR = Eigen::Matrix3f::Identity() + W;  // 一阶近似
    } else {
        // Rodrigues公式
        deltaR = I + W*sin(d)/d + W*W*(1.0f-cos(d))/d2;
    }
}
```

### 4.2 SO(3)对数映射

从旋转矩阵到李代数：

$$\theta = \arccos\left(\frac{\text{tr}(\mathbf{R})-1}{2}\right)$$

$$\boldsymbol{\phi} = \frac{\theta}{2\sin\theta}(\mathbf{R} - \mathbf{R}^T)^\vee$$

### 4.3 SO(3)右雅可比矩阵

右雅可比矩阵定义了李代数增量与旋转矩阵增量的关系：

$$\text{Exp}(\boldsymbol{\phi} + \delta\boldsymbol{\phi}) \approx \text{Exp}(\boldsymbol{\phi}) \cdot \text{Exp}(\mathbf{J}_r(\boldsymbol{\phi}) \delta\boldsymbol{\phi})$$

**封闭形式**：

$$\mathbf{J}_r(\boldsymbol{\phi}) = \mathbf{I} - \frac{1-\cos\theta}{\theta^2}\lfloor\boldsymbol{\phi}\rfloor_\times + \frac{\theta - \sin\theta}{\theta^3}\lfloor\boldsymbol{\phi}\rfloor_\times^2$$

当 $\theta \to 0$ 时：$\mathbf{J}_r \approx \mathbf{I}$

代码实现：
```cpp
// src/ImuTypes.cc:39-54
Eigen::Matrix3f RightJacobianSO3(const float &x, const float &y, const float &z) {
    Eigen::Matrix3f I;
    I.setIdentity();
    const float d2 = x*x + y*y + z*z;  // θ²
    const float d = sqrt(d2);           // θ

    Eigen::Vector3f v;
    v << x, y, z;
    Eigen::Matrix3f W = Sophus::SO3f::hat(v);  // [φ]×

    if(d < eps) {
        return I;
    } else {
        return I - W*(1.0f-cos(d))/d2 + W*W*(d-sin(d))/(d2*d);
    }
}
```

**右雅可比的逆**：

$$\mathbf{J}_r^{-1}(\boldsymbol{\phi}) = \mathbf{I} + \frac{1}{2}\lfloor\boldsymbol{\phi}\rfloor_\times + \left(\frac{1}{\theta^2} - \frac{1+\cos\theta}{2\theta\sin\theta}\right)\lfloor\boldsymbol{\phi}\rfloor_\times^2$$

```cpp
// src/ImuTypes.cc:61-77
Eigen::Matrix3f InverseRightJacobianSO3(const float &x, const float &y, const float &z) {
    if(d < eps) {
        return I;
    } else {
        return I + W/2 + W*W*(1.0f/d2 - (1.0f+cos(d))/(2.0f*d*sin(d)));
    }
}
```

### 4.4 BCH公式（Baker-Campbell-Hausdorff）

两个李代数元素对应旋转的乘积：

$$\ln(\text{Exp}(\boldsymbol{\phi}_1) \cdot \text{Exp}(\boldsymbol{\phi}_2)) \approx \begin{cases} \mathbf{J}_l^{-1}(\boldsymbol{\phi}_1)\boldsymbol{\phi}_2 + \boldsymbol{\phi}_1 & \text{if } \boldsymbol{\phi}_2 \text{ small} \\ \mathbf{J}_r^{-1}(\boldsymbol{\phi}_2)\boldsymbol{\phi}_1 + \boldsymbol{\phi}_2 & \text{if } \boldsymbol{\phi}_1 \text{ small} \end{cases}$$

其中左雅可比 $\mathbf{J}_l(\boldsymbol{\phi}) = \mathbf{J}_r(-\boldsymbol{\phi})$

---

## 5. Sim3相似变换求解

### 5.1 Sim3群定义

Sim3（3D相似变换群）包含：旋转 $\mathbf{R}$、平移 $\mathbf{t}$、尺度 $s$

$$\text{Sim}(3) = \left\{ \begin{bmatrix} s\mathbf{R} & \mathbf{t} \\ \mathbf{0}^T & 1 \end{bmatrix} \middle| \mathbf{R} \in SO(3), \mathbf{t} \in \mathbb{R}^3, s > 0 \right\}$$

**作用于3D点**：
$$\mathbf{p}' = s\mathbf{R}\mathbf{p} + \mathbf{t}$$

### 5.2 Horn闭式解算法

位置：`src/Sim3Solver.cc:311-412`

给定两组匹配的3D点 $\{\mathbf{p}_i\}$ 和 $\{\mathbf{q}_i\}$，求解相似变换：

**Step 1: 计算质心并去中心化**

$$\bar{\mathbf{p}} = \frac{1}{n}\sum_{i=1}^n \mathbf{p}_i, \quad \bar{\mathbf{q}} = \frac{1}{n}\sum_{i=1}^n \mathbf{q}_i$$

$$\mathbf{p}'_i = \mathbf{p}_i - \bar{\mathbf{p}}, \quad \mathbf{q}'_i = \mathbf{q}_i - \bar{\mathbf{q}}$$

```cpp
void Sim3Solver::ComputeCentroid(Eigen::Matrix3f &P, Eigen::Matrix3f &Pr, Eigen::Vector3f &C) {
    C = P.rowwise().sum();
    C = C / P.cols();  // 质心
    for(int i=0; i<P.cols(); i++)
        Pr.col(i) = P.col(i) - C;  // 去中心化
}
```

**Step 2: 构造M矩阵**

$$\mathbf{M} = \sum_{i=1}^n \mathbf{q}'_i {\mathbf{p}'_i}^T = \mathbf{Q}' {\mathbf{P}'}^T$$

```cpp
Eigen::Matrix3f M = Pr2 * Pr1.transpose();
```

**Step 3: 构造N矩阵（4×4对称矩阵）**

将 $\mathbf{M}$ 的9个元素重新排列为用于四元数求解的矩阵：

$$\mathbf{N} = \begin{bmatrix} S_{xx}+S_{yy}+S_{zz} & S_{yz}-S_{zy} & S_{zx}-S_{xz} & S_{xy}-S_{yx} \\ S_{yz}-S_{zy} & S_{xx}-S_{yy}-S_{zz} & S_{xy}+S_{yx} & S_{zx}+S_{xz} \\ S_{zx}-S_{xz} & S_{xy}+S_{yx} & -S_{xx}+S_{yy}-S_{zz} & S_{yz}+S_{zy} \\ S_{xy}-S_{yx} & S_{zx}+S_{xz} & S_{yz}+S_{zy} & -S_{xx}-S_{yy}+S_{zz} \end{bmatrix}$$

其中 $S_{ab} = M_{ab}$

```cpp
N11 = M(0,0)+M(1,1)+M(2,2);
N12 = M(1,2)-M(2,1);
N13 = M(2,0)-M(0,2);
N14 = M(0,1)-M(1,0);
// ... 其他元素
```

**Step 4: 求最大特征值对应的特征向量**

$\mathbf{N}$ 的最大特征值对应的特征向量即为旋转四元数 $\mathbf{q} = [q_w, q_x, q_y, q_z]^T$

```cpp
Eigen::EigenSolver<Eigen::Matrix4f> eigSolver;
eigSolver.compute(N);

Eigen::Vector4f eval = eigSolver.eigenvalues().real();
Eigen::Matrix4f evec = eigSolver.eigenvectors().real();

int maxIndex;
eval.maxCoeff(&maxIndex);  // 找最大特征值

// 从四元数恢复旋转矩阵
Eigen::Vector3f vec = evec.block<3,1>(1,maxIndex);  // 虚部
double ang = atan2(vec.norm(), evec(0,maxIndex));   // 半角
vec = 2*ang*vec/vec.norm();  // 转为轴角
mR12i = Sophus::SO3f::exp(vec).matrix();  // 转为旋转矩阵
```

**Step 5: 计算尺度**

$$s = \frac{\sum_i \mathbf{p}'_i \cdot (\mathbf{R}\mathbf{q}'_i)}{\sum_i \|\mathbf{R}\mathbf{q}'_i\|^2}$$

```cpp
if(!mbFixScale) {
    Eigen::Matrix3f P3 = mR12i * Pr2;  // R * q'
    double nom = (Pr1.array() * P3.array()).sum();  // p' · (R*q')
    double den = (P3.array() * P3.array()).sum();   // ||R*q'||²
    ms12i = nom / den;
} else {
    ms12i = 1.0f;  // 双目/RGB-D固定尺度为1
}
```

**Step 6: 计算平移**

$$\mathbf{t} = \bar{\mathbf{p}} - s\mathbf{R}\bar{\mathbf{q}}$$

```cpp
mt12i = O1 - ms12i * mR12i * O2;
```

### 5.3 RANSAC框架

位置：`src/Sim3Solver.cc:149-216`

```
┌─────────────────────────────────────────────────────────────┐
│                   Sim3 RANSAC 流程                           │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│  while (迭代次数 < 最大迭代 && 当前迭代 < 本次迭代限制)     │
│  {                                                          │
│      1. 随机选择3对匹配点（最小集）                         │
│      2. Horn闭式解计算Sim3                                  │
│      3. 检查内点数（投影误差 < 阈值）                       │
│      4. 更新最佳结果                                        │
│      5. 如果内点数 > 最小内点数，提前返回                   │
│  }                                                          │
└─────────────────────────────────────────────────────────────┘
```

内点判断（双向重投影误差）：
```cpp
void Sim3Solver::CheckInliers() {
    Project(mvX3Dc2, vP2im1, mT12i, pCamera1);  // q通过T12投影到帧1
    Project(mvX3Dc1, vP1im2, mT21i, pCamera2);  // p通过T21投影到帧2

    for(size_t i=0; i<mvP1im1.size(); i++) {
        float err1 = (mvP1im1[i] - vP2im1[i]).squaredNorm();
        float err2 = (vP1im2[i] - mvP2im2[i]).squaredNorm();

        if(err1 < mvnMaxError1[i] && err2 < mvnMaxError2[i])
            mvbInliersi[i] = true;
    }
}
```

---

## 6. 因子图优化框架

### 6.1 g2o优化器概述

ORB-SLAM3使用g2o（General Graph Optimization）作为后端优化框架。核心概念：

- **顶点（Vertex）**：优化变量（位姿、速度、偏置、地图点）
- **边（Edge）**：约束/残差（重投影误差、IMU预积分误差）
- **信息矩阵**：残差的权重（协方差的逆）

### 6.2 顶点类型

位置：`include/G2oTypes.h`

| 顶点类 | 维度 | 状态量 | 用途 |
|--------|------|--------|------|
| `VertexPose` | 6 | $\mathbf{R}_{wb}, \mathbf{t}_{wb}$ | IMU坐标系位姿 |
| `VertexVelocity` | 3 | $\mathbf{v}_w$ | 世界系速度 |
| `VertexGyroBias` | 3 | $\mathbf{b}_g$ | 陀螺仪偏置 |
| `VertexAccBias` | 3 | $\mathbf{b}_a$ | 加速度计偏置 |
| `VertexGDir` | 2 | $\mathbf{R}_{wg}$ | 重力方向（2DoF） |
| `VertexScale` | 1 | $s$ | 尺度因子 |
| `VertexPose4DoF` | 4 | $\theta_{yaw}, \mathbf{t}$ | 4DoF位姿（回环） |

**位姿增量更新**（李代数参数化）：
```cpp
// include/G2oTypes.h:149-152
virtual void oplusImpl(const double* update_) {
    _estimate.Update(update_);  // 应用6DoF增量
    updateCache();
}

// ImuCamPose::Update - 李代数增量
void ImuCamPose::Update(const double *pu) {
    Eigen::Vector3d phi;
    phi << pu[0], pu[1], pu[2];  // 旋转增量
    Rwb = Rwb * ExpSO3(phi);      // R_new = R_old * Exp(δφ)

    twb += Eigen::Vector3d(pu[3], pu[4], pu[5]);  // t_new = t_old + δt
}
```

### 6.3 视觉重投影边

**单目重投影误差**：

$$\mathbf{e}_{reproj} = \mathbf{u} - \pi(\mathbf{T}_{cw} \mathbf{P}_w)$$

其中 $\pi$ 为相机投影函数，$\mathbf{u}$ 为观测像素坐标。

```cpp
// include/G2oTypes.h:353-358
class EdgeMono : public g2o::BaseBinaryEdge<2, Eigen::Vector2d,
                                             g2o::VertexSBAPointXYZ,
                                             VertexPose> {
    void computeError() {
        const g2o::VertexSBAPointXYZ* VPoint = ...;
        const VertexPose* VPose = ...;
        _error = obs - VPose->estimate().Project(VPoint->estimate(), cam_idx);
    }
};
```

**双目重投影误差**（3维：$u_l, v, u_r$）：

```cpp
// include/G2oTypes.h:435-440
class EdgeStereo : public g2o::BaseBinaryEdge<3, Eigen::Vector3d, ...> {
    void computeError() {
        _error = obs - VPose->estimate().ProjectStereo(VPoint->estimate(), cam_idx);
    }
};
```

### 6.4 IMU预积分边

位置：`include/G2oTypes.h:495-544`

**残差定义（9维）**：

$$\mathbf{r}_{IMU} = \begin{bmatrix} \mathbf{r}_R \\ \mathbf{r}_v \\ \mathbf{r}_p \end{bmatrix} = \begin{bmatrix} \text{Log}\big(\Delta\tilde{\mathbf{R}}_{ij}^T \mathbf{R}_i^T \mathbf{R}_j\big) \\ \mathbf{R}_i^T(\mathbf{v}_j - \mathbf{v}_i - \mathbf{g}\Delta t) - \Delta\tilde{\mathbf{v}}_{ij} \\ \mathbf{R}_i^T(\mathbf{p}_j - \mathbf{p}_i - \mathbf{v}_i\Delta t - \frac{1}{2}\mathbf{g}\Delta t^2) - \Delta\tilde{\mathbf{p}}_{ij} \end{bmatrix}$$

```cpp
class EdgeInertial : public g2o::BaseMultiEdge<9, Vector9d> {
    // 连接6个顶点：Pose_i, Vel_i, GyroBias_i, AccBias_i, Pose_j, Vel_j

    void computeError() {
        // 获取顶点状态
        const VertexPose* VP1 = ...;
        const VertexVelocity* VV1 = ...;
        // ...

        // 获取修正后的预积分量
        const Eigen::Matrix3d dR = mpInt->GetDeltaRotation(b).cast<double>();
        const Eigen::Vector3d dV = mpInt->GetDeltaVelocity(b).cast<double>();
        const Eigen::Vector3d dP = mpInt->GetDeltaPosition(b).cast<double>();

        // 计算残差
        const Eigen::Vector3d er = LogSO3(dR.transpose() * Rwb1.transpose() * Rwb2);
        const Eigen::Vector3d ev = Rwb1.transpose() * (Vwb2 - Vwb1 - g*dt) - dV;
        const Eigen::Vector3d ep = Rwb1.transpose() * (twb2 - twb1 - Vwb1*dt - 0.5*g*dt*dt) - dP;

        _error << er, ev, ep;
    }
};
```

### 6.5 偏置随机游走边

偏置约束为相邻帧偏置应该接近：

$$\mathbf{r}_{bias} = \mathbf{b}_{g,j} - \mathbf{b}_{g,i}$$

信息矩阵为 $\frac{1}{\sigma_{walk}^2 \Delta t}\mathbf{I}$

```cpp
// include/G2oTypes.h:635-668
class EdgeGyroRW : public g2o::BaseBinaryEdge<3, Eigen::Vector3d,
                                               VertexGyroBias,
                                               VertexGyroBias> {
    void computeError() {
        const VertexGyroBias* VG1 = ...;
        const VertexGyroBias* VG2 = ...;
        _error = VG2->estimate() - VG1->estimate();
    }

    void linearizeOplus() {
        _jacobianOplusXi = -Eigen::Matrix3d::Identity();
        _jacobianOplusXj.setIdentity();
    }
};
```

### 6.6 因子图结构

视觉惯性SLAM的完整因子图：

```
     b_g^0    b_g^1    b_g^2         (陀螺仪偏置)
       ●───────●───────●
       │       │       │
       │  RW   │  RW   │             (随机游走约束)
       │       │       │
     b_a^0    b_a^1    b_a^2         (加速度计偏置)
       ●───────●───────●
       │       │       │
       │       │       │
       ▼       ▼       ▼
    ┌─────┐ ┌─────┐ ┌─────┐
    │ x_0 │ │ x_1 │ │ x_2 │         (位姿)
    └──┬──┘ └──┬──┘ └──┬──┘
       │       │       │
       │  IMU  │  IMU  │             (IMU预积分约束)
       └───────┴───────┘
       │       │       │
    ┌──┴──┐ ┌──┴──┐ ┌──┴──┐
    │ v_0 │ │ v_1 │ │ v_2 │         (速度)
    └─────┘ └─────┘ └─────┘

       │       │       │
       ▼       ▼       ▼
       ●───────●───────●───────●     (重投影约束)
      P_1     P_2     P_3     P_4    (地图点)
```

---

## 7. 系统初始化流程

### 7.1 视觉初始化

**单目初始化**：两视图重建

位置：`src/Tracking.cc:MonocularInitialization()`

1. 提取两帧的ORB特征并匹配
2. 并行计算本质矩阵E和单应矩阵H
3. 选择更好的模型（根据对称重投影误差）
4. 从E/H恢复R,t
5. 三角化初始地图点

**双目/RGB-D初始化**：直接三角化

位置：`src/Tracking.cc:StereoInitialization()`

1. 提取ORB特征
2. 双目匹配或深度图恢复3D点
3. 创建初始地图

### 7.2 IMU初始化

位置：`src/LocalMapping.cc:1173-1427`

IMU初始化分三个阶段：

**阶段一：粗略初始化（5秒）**
- 使用视觉估计的速度来估计重力方向
- 使用强先验约束偏置

```cpp
// 估计重力方向
for (auto itKF = vpKF.begin(); itKF != vpKF.end(); itKF++) {
    // 从预积分速度累积重力方向
    dirG -= (*itKF)->mPrevKF->GetImuRotation()
          * (*itKF)->mpImuPreintegrated->GetUpdatedDeltaVelocity();
}
dirG = dirG / dirG.norm();

// 计算 IMU系 到 世界系 的旋转 Rwg
Eigen::Vector3f gI(0.0f, 0.0f, -1.0f);  // IMU系中的重力（假设z轴向上）
Eigen::Vector3f v = gI.cross(dirG);
float ang = acos(gI.dot(dirG));
Rwg = Sophus::SO3f::exp(v * ang / v.norm()).matrix();
```

**阶段二：优化（15秒）**
- 联合优化重力方向、尺度、速度、偏置
- 移除先验约束

**阶段三：尺度精化（持续）**
- 周期性重新优化以提高尺度精度

### 7.3 尺度和重力对齐

当IMU初始化完成后，需要将地图对齐到重力方向：

$$\mathbf{p}'_w = s \cdot \mathbf{R}_{wg}^T \mathbf{p}_w$$

```cpp
if (fabs(mScale - 1.f) > 0.00001 || !mbMonocular) {
    Sophus::SE3f Twg(mRwg.transpose(), Eigen::Vector3f::Zero());
    mpAtlas->GetCurrentMap()->ApplyScaledRotation(Twg, mScale, true);
}
```

---

## 8. 回环检测与地图融合

### 8.1 回环检测流程

位置：`src/LoopClosing.cc:NewDetectCommonRegions()`

```
┌─────────────────────────────────────────────────────────────┐
│              NewDetectCommonRegions() 主流程                 │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│  1. 检查新关键帧                                             │
│     - IMU未完成初始化时跳过                                  │
│     - 地图关键帧数 < 12 时跳过                               │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│  2. 检测候选关键帧                                           │
│     - 基于上一帧候选继续验证（时序一致性）                    │
│     - 或从词袋数据库查询新候选                               │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│  3. Sim3几何验证                                             │
│     - 计算当前帧与候选帧的Sim3变换                           │
│     - 使用RANSAC剔除外点                                     │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│  4. 时序一致性检验                                           │
│     - 连续3帧检测到同一回环才确认                            │
│     mnLoopNumCoincidences >= 3                              │
└─────────────────────────────────────────────────────────────┘
                         │
            ┌────────────┴────────────┐
            ▼                         ▼
    ┌───────────────┐         ┌───────────────┐
    │ 同一地图回环   │         │ 跨地图融合     │
    │ CorrectLoop() │         │ MergeLocal()  │
    └───────────────┘         └───────────────┘
```

### 8.2 回环矫正 CorrectLoop()

位置：`src/LoopClosing.cc:1270-1600`

1. **停止LocalMapping**
2. **计算Sim3位姿修正**
3. **传播位姿修正到当前地图所有关键帧**
4. **融合地图点**
5. **本质图优化（Essential Graph Optimization）**
6. **启动后台全局BA**

```cpp
void LoopClosing::CorrectLoop() {
    // 停止局部建图
    mpLocalMapper->RequestStop();

    // 计算位姿修正
    g2o::Sim3 g2oCorrectedScw = mg2oLoopScw;

    // 传播修正到所有共视关键帧
    for (auto pKFi : mpCurrentKF->GetVectorCovisibleKeyFrames()) {
        g2o::Sim3 g2oCorrectedSiw = g2oSiw * g2oSwc * g2oCorrectedScw;
        CorrectedSim3[pKFi] = g2oCorrectedSiw;

        // 更新地图点位置
        for (auto pMPi : pKFi->GetMapPointMatches()) {
            Eigen::Vector3d correctedP3Dw = g2oCorrectedSiw.map(
                g2oSiw.inverse().map(P3Dw));
            pMPi->SetWorldPos(correctedP3Dw);
        }
    }

    // 地图点融合
    SearchAndFuse(CorrectedSim3, mvpLoopMapPoints);

    // 本质图优化
    Optimizer::OptimizeEssentialGraph(...);

    // 后台全局BA
    mpThreadGBA = new thread(&LoopClosing::RunGlobalBundleAdjustment, ...);
}
```

### 8.3 地图融合 MergeLocal()

位置：`src/LoopClosing.cc:850-1200`

当检测到跨地图的公共区域时，进行地图融合：

1. **计算融合变换Sim3**
2. **传播变换到当前地图所有元素**
3. **融合关键帧和地图点**
4. **Welding BA**（焊接区域的联合优化）
5. **删除旧地图**

### 8.4 本质图 (Essential Graph)

本质图是共视图的稀疏化版本，用于高效的位姿优化：

- **生成树边**：关键帧之间的父子关系
- **回环边**：回环检测产生的约束
- **强共视边**：共视点数 > 100 的边

本质图优化仅优化关键帧位姿，固定地图点，使用Sim3/SE3变换。

---

## 9. 核心数据结构

### 9.1 Frame（普通帧）

```cpp
class Frame {
    // 时间戳
    double mTimeStamp;

    // ORB特征
    std::vector<cv::KeyPoint> mvKeys;        // 原始特征点
    std::vector<cv::KeyPoint> mvKeysUn;      // 畸变校正后
    cv::Mat mDescriptors;                     // 描述子

    // 位姿 (世界→相机)
    Sophus::SE3f mTcw;

    // 地图点对应
    std::vector<MapPoint*> mvpMapPoints;
    std::vector<bool> mvbOutlier;

    // IMU相关
    IMU::Preintegrated* mpImuPreintegrated;
    IMU::Bias mImuBias;
    Eigen::Vector3f mVw;  // 世界系速度
};
```

### 9.2 KeyFrame（关键帧）

在Frame基础上增加：

```cpp
class KeyFrame {
    // 共视图
    std::map<KeyFrame*, int> mConnectedKeyFrameWeights;
    std::vector<KeyFrame*> mvpOrderedConnectedKeyFrames;

    // 生成树
    KeyFrame* mpParent;
    std::set<KeyFrame*> mspChildrens;

    // 回环边
    std::set<KeyFrame*> mspLoopEdges;

    // IMU相邻关键帧
    KeyFrame* mPrevKF;
    KeyFrame* mNextKF;
};
```

### 9.3 MapPoint（地图点）

```cpp
class MapPoint {
    // 3D坐标（世界系）
    Eigen::Vector3f mWorldPos;

    // 观测信息
    std::map<KeyFrame*, std::tuple<int,int>> mObservations;  // KF→特征索引

    // 描述子（所有观测的中位数）
    cv::Mat mDescriptor;

    // 平均观测方向
    Eigen::Vector3f mNormalVector;

    // 尺度不变性范围
    float mfMinDistance;
    float mfMaxDistance;
};
```

### 9.4 Atlas（多地图管理）

```cpp
class Atlas {
    std::set<Map*> mspMaps;        // 所有地图
    Map* mpCurrentMap;             // 当前活跃地图
    std::set<Map*> mspBadMaps;     // 已删除的地图

    // 地图操作
    void CreateNewMap();           // 跟踪丢失时
    void ChangeMap(Map* pMap);     // 重定位成功时
    void SetBadMap(Map* pMap);     // 地图融合后
};
```

### 9.5 IMU::Preintegrated（预积分）

```cpp
class Preintegrated {
    // 预积分量
    float dT;                              // 总积分时间
    Eigen::Matrix3f dR;                    // ΔR
    Eigen::Vector3f dV, dP;                // Δv, Δp

    // 协方差
    Eigen::Matrix<float,15,15> C;          // 15×15协方差

    // Bias雅可比
    Eigen::Matrix3f JRg;                   // ∂ΔR/∂b_g
    Eigen::Matrix3f JVg, JVa;              // ∂Δv/∂b_g, ∂Δv/∂b_a
    Eigen::Matrix3f JPg, JPa;              // ∂Δp/∂b_g, ∂Δp/∂b_a

    // Bias
    Bias b;                                // 线性化点
    Bias bu;                               // 当前估计
    Eigen::Matrix<float,6,1> db;           // 偏置变化量
};
```

---

## 10. 算法对比分析

### 10.1 ORB-SLAM3 vs VINS

| 方面 | ORB-SLAM3 | VINS-Mono/Fusion |
|------|-----------|------------------|
| **前端特征** | ORB描述子（每帧重提取） | FAST + LK光流跟踪 |
| **后端优化** | g2o | Ceres |
| **地图管理** | Atlas多地图 | 单地图 |
| **回环检测** | DBoW2 | DBoW2 |
| **地图点参数化** | 3D世界坐标 | 逆深度 |
| **IMU积分** | 欧拉法 | 中点法 |
| **边缘化** | 关键帧剔除 | 滑动窗口边缘化 |

### 10.2 IMU预积分对比

**ORB-SLAM3（欧拉法）**：
```cpp
dP = dP + dV*dt + 0.5f*dR*acc*dt*dt;
dV = dV + dR*acc*dt;
// 然后更新 dR
```

**VINS（中点法）**：
```cpp
Vector3d un_acc_0 = delta_q * (_acc_0 - linearized_ba);
result_delta_q = delta_q * Quaterniond(1, un_gyr*dt/2);
Vector3d un_acc_1 = result_delta_q * (_acc_1 - linearized_ba);
Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);  // 中点
result_delta_p = delta_p + delta_v*dt + 0.5*un_acc*dt*dt;
result_delta_v = delta_v + un_acc*dt;
```

**理论差异**：
- 中点法精度更高（二阶收敛 vs 一阶收敛）
- 欧拉法计算量更小
- 在200Hz IMU采样率下，实际差异较小

### 10.3 特征点管理对比

| 对比项 | ORB-SLAM3 | VINS |
|--------|-----------|------|
| 检测 | 每帧ORB | 增量式FAST |
| 跟踪 | 描述子匹配 | LK光流 |
| 重定位 | BoW查询 | 依赖回环 |
| 优势 | 重定位鲁棒 | 跟踪效率高 |

### 10.4 优化框架对比

**g2o（ORB-SLAM3）**：
- 显式定义顶点和边类
- 稀疏矩阵优化
- 支持Sim3优化

**Ceres（VINS）**：
- 残差块 + 参数块
- 自动求导支持
- 更灵活的约束表达

---

## 参考文献

1. Campos, C., et al. "ORB-SLAM3: An Accurate Open-Source Library for Visual, Visual-Inertial and Multi-Map SLAM." IEEE T-RO, 2021.

2. Mur-Artal, R., et al. "ORB-SLAM: A Versatile and Accurate Monocular SLAM System." IEEE T-RO, 2015.

3. Forster, C., et al. "On-Manifold Preintegration for Real-Time Visual-Inertial Odometry." IEEE T-RO, 2017.

4. Horn, B. K. P. "Closed-form solution of absolute orientation using unit quaternions." JOSA A, 1987.

5. Rublee, E., et al. "ORB: An efficient alternative to SIFT or SURF." ICCV, 2011.

---

*文档版本：v2.0*
*更新日期：2026-02-01*
*基于 ORB-SLAM3 源码深度分析*
