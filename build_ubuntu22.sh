#!/bin/bash

# ORB-SLAM3 Build Script for Ubuntu 22.04
# 适配 OpenCV 4.x + C++14

set -e  # 遇到错误立即退出

echo "=========================================="
echo "ORB-SLAM3 Build Script for Ubuntu 22.04"
echo "=========================================="

# 获取 CPU 核心数
NUM_CORES=$(nproc)
echo "Using $NUM_CORES cores for compilation"

cd "$(dirname "$0")"
SCRIPT_DIR=$(pwd)

# 1. 编译 DBoW2
echo ""
echo "[1/4] Building Thirdparty/DBoW2 ..."
cd Thirdparty/DBoW2
rm -rf build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j${NUM_CORES}

# 2. 编译 g2o
echo ""
echo "[2/4] Building Thirdparty/g2o ..."
cd ../../g2o
rm -rf build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j${NUM_CORES}

# 3. 编译 Sophus (仅头文件库，但仍需 cmake)
echo ""
echo "[3/4] Building Thirdparty/Sophus ..."
cd ../../Sophus
rm -rf build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j${NUM_CORES}

# 4. 解压词汇表
echo ""
echo "[3.5/4] Extracting vocabulary ..."
cd ${SCRIPT_DIR}/Vocabulary
if [ ! -f ORBvoc.txt ]; then
    tar -xf ORBvoc.txt.tar.gz
    echo "Vocabulary extracted."
else
    echo "Vocabulary already exists."
fi

# 5. 编译 ORB-SLAM3 主库
echo ""
echo "[4/4] Building ORB-SLAM3 ..."
cd ${SCRIPT_DIR}
rm -rf build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j${NUM_CORES}

echo ""
echo "=========================================="
echo "ORB-SLAM3 Build Complete!"
echo "=========================================="
echo ""
echo "Executable examples in: ${SCRIPT_DIR}/Examples/"
echo ""
echo "Test with EuRoC dataset:"
echo "  ./Examples/Monocular-Inertial/mono_inertial_euroc ./Vocabulary/ORBvoc.txt ./Examples/Monocular-Inertial/EuRoC.yaml /path/to/MH01 ./Examples/Monocular-Inertial/EuRoC_TimeStamps/MH01.txt"
echo ""
