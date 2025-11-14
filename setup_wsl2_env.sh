#!/bin/bash
# WSL2 Ubuntu-24.04 环境设置脚本
# 用于自动安装 Miniconda 和设置项目环境

set -e  # 遇到错误立即退出

echo "=========================================="
echo "WSL2 环境设置脚本"
echo "=========================================="
echo ""

# 颜色定义
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# 检查是否已安装 Miniconda
if [ -d "$HOME/miniconda3" ] || command -v conda &> /dev/null; then
    echo -e "${GREEN}✓ Miniconda 已安装${NC}"
    if [ -d "$HOME/miniconda3" ]; then
        export PATH="$HOME/miniconda3/bin:$PATH"
    fi
    if command -v conda &> /dev/null; then
        conda --version
    fi
else
    echo -e "${YELLOW}正在安装 Miniconda...${NC}"
    
    # 设置代理（如果需要）
    if [ -n "$http_proxy" ]; then
        export http_proxy="$http_proxy"
        export https_proxy="$https_proxy"
    fi
    
    # 下载 Miniconda
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    MINICONDA_INSTALLER="/tmp/miniconda.sh"
    
    echo "下载 Miniconda..."
    wget -q "$MINICONDA_URL" -O "$MINICONDA_INSTALLER" || {
        echo -e "${RED}下载失败，请检查网络连接${NC}"
        exit 1
    }
    
    # 安装 Miniconda（非交互式）
    echo "安装 Miniconda..."
    bash "$MINICONDA_INSTALLER" -b -p "$HOME/miniconda3" || {
        echo -e "${RED}安装失败${NC}"
        exit 1
    }
    
    # 初始化 conda
    "$HOME/miniconda3/bin/conda" init bash
    
    # 添加到 PATH
    export PATH="$HOME/miniconda3/bin:$PATH"
    
    # 清理安装文件
    rm "$MINICONDA_INSTALLER"
    
    echo -e "${GREEN}✓ Miniconda 安装完成${NC}"
fi

# 重新加载 shell 配置
source "$HOME/.bashrc" 2>/dev/null || true

# 确保 conda 在 PATH 中
export PATH="$HOME/miniconda3/bin:$PATH"

# 接受 conda 服务条款
echo ""
echo -e "${YELLOW}接受 conda 服务条款...${NC}"
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

# 创建 conda 环境
ENV_NAME="flash_dock_wsl"
PYTHON_VERSION="3.10"

echo ""
echo -e "${YELLOW}检查 conda 环境: $ENV_NAME${NC}"

if conda env list | grep -q "^$ENV_NAME "; then
    echo -e "${GREEN}✓ 环境 $ENV_NAME 已存在${NC}"
else
    echo -e "${YELLOW}创建 conda 环境: $ENV_NAME (Python $PYTHON_VERSION)${NC}"
    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
    
    echo -e "${GREEN}✓ 环境创建完成${NC}"
fi

# 激活环境
echo ""
echo -e "${YELLOW}激活环境并安装依赖...${NC}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# 配置 pip 超时和重试
export PIP_DEFAULT_TIMEOUT=100
export PIP_RETRIES=10

# 升级 pip（增加超时和重试）
echo "升级 pip..."
pip install --upgrade pip --timeout=100 --retries=10 -q || {
    echo -e "${YELLOW}⚠ pip 升级失败，尝试使用 conda 安装 pip...${NC}"
    conda install pip -y -q
}

# 安装 PyTorch (CPU 版本)
echo "安装 PyTorch (CPU 版本)..."
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu --timeout=100 --retries=10 -q || {
    echo -e "${YELLOW}⚠ PyTorch 安装失败，尝试使用 conda 安装...${NC}"
    conda install pytorch torchvision torchaudio cpuonly -c pytorch -y -q
}

# 安装项目依赖
echo "安装项目依赖..."
PROJECT_ROOT="/mnt/e/Qinchaojun/AIDD-TRAIN"
REQUIREMENTS_FILE="$PROJECT_ROOT/requirements.txt"

if [ -f "$REQUIREMENTS_FILE" ]; then
    pip install -r "$REQUIREMENTS_FILE" --timeout=100 --retries=10 -q || {
        echo -e "${YELLOW}⚠ 部分依赖安装失败，继续安装其他依赖...${NC}"
    }
    echo -e "${GREEN}✓ 项目依赖安装完成${NC}"
else
    echo -e "${YELLOW}⚠ 未找到 requirements.txt，跳过${NC}"
fi

# 安装 Uni-Mol 特定依赖
echo "安装 Uni-Mol 特定依赖..."
pip install lmdb rdkit-pypi tqdm scikit-learn pandas --timeout=100 --retries=10 -q || {
    echo -e "${YELLOW}⚠ 部分依赖安装失败，尝试逐个安装...${NC}"
    for pkg in lmdb rdkit-pypi tqdm scikit-learn pandas; do
        pip install "$pkg" --timeout=100 --retries=10 -q || echo -e "${YELLOW}⚠ $pkg 安装失败${NC}"
    done
}

# 安装 Uni-Core
echo "安装 Uni-Core..."
pip install "git+https://github.com/dptech-corp/Uni-Core.git@stable" --timeout=200 --retries=10 -q || {
    echo -e "${YELLOW}⚠ Uni-Core 安装失败，可能需要手动安装${NC}"
    echo -e "${YELLOW}  手动安装命令: pip install git+https://github.com/dptech-corp/Uni-Core.git@stable${NC}"
}

echo ""
echo -e "${GREEN}=========================================="
echo "环境设置完成！"
echo "==========================================${NC}"
echo ""
echo "使用方法："
echo "1. 激活环境: conda activate $ENV_NAME"
echo "2. 进入项目目录: cd $PROJECT_ROOT/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface"
echo "3. 运行程序: python3 demo.py [参数]"
echo ""

