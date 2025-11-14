#!/bin/bash
# WSL2 环境设置脚本
# 用于在 WSL2 中快速设置 Uni-Mol 对接预测程序所需的环境

set -e  # 遇到错误立即退出

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Uni-Mol Docking WSL2 环境设置脚本 ===${NC}"
echo ""

# 检查是否在 WSL2 环境中
if [ ! -f /proc/version ] || ! grep -qi "microsoft\|wsl" /proc/version; then
    echo -e "${YELLOW}警告: 可能不在 WSL2 环境中运行${NC}"
    read -p "是否继续? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# 检查 conda 是否已安装
if ! command -v conda &> /dev/null; then
    echo -e "${YELLOW}未检测到 conda，开始安装 Miniconda...${NC}"
    
    # 下载 Miniconda
    echo -e "${BLUE}下载 Miniconda...${NC}"
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    
    # 安装 Miniconda
    echo -e "${BLUE}安装 Miniconda...${NC}"
    bash /tmp/miniconda.sh -b -p $HOME/miniconda3
    
    # 初始化 conda
    $HOME/miniconda3/bin/conda init bash
    source $HOME/.bashrc
    
    # 添加到 PATH
    export PATH="$HOME/miniconda3/bin:$PATH"
    
    echo -e "${GREEN}Miniconda 安装完成${NC}"
else
    echo -e "${GREEN}检测到 conda: $(conda --version)${NC}"
    # 确保 conda 在 PATH 中
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    fi
fi

# 检查 flash_dock 环境是否存在
if conda env list | grep -q "flash_dock"; then
    echo -e "${YELLOW}检测到 flash_dock 环境已存在${NC}"
    read -p "是否重新创建? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${BLUE}删除现有环境...${NC}"
        conda env remove -n flash_dock -y
    else
        echo -e "${GREEN}使用现有环境${NC}"
        conda activate flash_dock
        skip_install=true
    fi
fi

if [ "$skip_install" != "true" ]; then
    # 创建 conda 环境
    echo -e "${BLUE}创建 conda 环境 flash_dock (Python 3.10)...${NC}"
    conda create -n flash_dock python=3.10 -y
    
    # 激活环境
    echo -e "${BLUE}激活环境...${NC}"
    conda activate flash_dock
fi

# 检查 CUDA 是否可用（可选）
echo -e "${BLUE}检查 CUDA 可用性...${NC}"
if command -v nvidia-smi &> /dev/null; then
    echo -e "${GREEN}检测到 NVIDIA GPU${NC}"
    nvidia-smi
    USE_CUDA=true
else
    echo -e "${YELLOW}未检测到 NVIDIA GPU，将安装 CPU 版本的 PyTorch${NC}"
    USE_CUDA=false
fi

# 安装 PyTorch
echo -e "${BLUE}安装 PyTorch...${NC}"
if [ "$USE_CUDA" = true ]; then
    echo -e "${GREEN}安装 CUDA 版本的 PyTorch...${NC}"
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
else
    echo -e "${GREEN}安装 CPU 版本的 PyTorch...${NC}"
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
fi

# 安装项目依赖
echo -e "${BLUE}安装项目依赖...${NC}"
pip install lmdb numpy pandas rdkit-pypi tqdm scikit-learn

# 安装 unicore（可能需要从源码安装）
echo -e "${BLUE}安装 unicore...${NC}"
pip install unicore || {
    echo -e "${YELLOW}pip 安装 unicore 失败，尝试从源码安装...${NC}"
    # 如果 pip 安装失败，可能需要从源码安装
    # 这里假设用户已经下载了 unicore 源码
    echo -e "${YELLOW}请手动安装 unicore: pip install unicore 或从源码安装${NC}"
}

# 验证安装
echo -e "${BLUE}验证安装...${NC}"
python3 -c "import torch; print(f'PyTorch: {torch.__version__}')"
python3 -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
python3 -c "import lmdb; print('LMDB: OK')"
python3 -c "import numpy; print('NumPy: OK')"
python3 -c "import pandas; print('Pandas: OK')"
python3 -c "import rdkit; print('RDKit: OK')" || echo -e "${YELLOW}RDKit 可能未正确安装${NC}"

echo ""
echo -e "${GREEN}=== 环境设置完成 ===${NC}"
echo ""
echo -e "${BLUE}使用方法:${NC}"
echo "1. 激活环境: conda activate flash_dock"
echo "2. 进入项目目录: cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface"
echo "3. 运行程序: python3 demo.py [参数]"
echo ""
echo -e "${BLUE}或使用运行脚本:${NC}"
echo "bash run_in_wsl2.sh [参数]"
echo ""



