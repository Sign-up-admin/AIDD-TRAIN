#!/bin/bash
# FlashDock WSL 启动脚本
# 在 WSL 的 conda 虚拟环境中启动 FlashDock

set -e  # 遇到错误立即退出

# 颜色定义
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 配置
ENV_NAME="flash_dock_wsl"
PROJECT_ROOT="/mnt/e/Qinchaojun/AIDD-TRAIN"
FLASHDOCK_DIR="$PROJECT_ROOT/FLASH_DOCK-main"
PORT=8501

echo "=========================================="
echo "FlashDock WSL 启动脚本"
echo "=========================================="
echo ""

# 检查项目目录
if [ ! -d "$PROJECT_ROOT" ]; then
    echo -e "${RED}✗ 项目目录不存在: $PROJECT_ROOT${NC}"
    exit 1
fi

if [ ! -d "$FLASHDOCK_DIR" ]; then
    echo -e "${RED}✗ FlashDock 目录不存在: $FLASHDOCK_DIR${NC}"
    exit 1
fi

# 检查 FlashDock.py 文件
if [ ! -f "$FLASHDOCK_DIR/FlashDock.py" ]; then
    echo -e "${RED}✗ FlashDock.py 文件不存在: $FLASHDOCK_DIR/FlashDock.py${NC}"
    exit 1
fi

# 初始化 conda
if [ -d "$HOME/miniconda3" ]; then
    export PATH="$HOME/miniconda3/bin:$PATH"
    source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true
elif [ -d "$HOME/anaconda3" ]; then
    export PATH="$HOME/anaconda3/bin:$PATH"
    source "$HOME/anaconda3/etc/profile.d/conda.sh" 2>/dev/null || true
fi

# 检查 conda 是否可用
if ! command -v conda &> /dev/null; then
    echo -e "${RED}✗ Conda 未找到，请手动安装并配置 conda${NC}"
    exit 1
fi

# 检查环境是否存在
if ! conda env list | grep -q "^$ENV_NAME "; then
    echo -e "${RED}✗ Conda 环境 '$ENV_NAME' 不存在${NC}"
    echo -e "${YELLOW}请手动创建 conda 环境并安装依赖${NC}"
    exit 1
fi

# 激活环境
echo -e "${YELLOW}激活 conda 环境: $ENV_NAME${NC}"
conda activate "$ENV_NAME"

# 验证 Python 和关键依赖
echo -e "${YELLOW}验证环境...${NC}"
if ! python -c "import streamlit" 2>/dev/null; then
    echo -e "${RED}✗ Streamlit 未安装${NC}"
    echo -e "${YELLOW}请手动安装依赖包${NC}"
    exit 1
fi

if ! python -c "import rdkit" 2>/dev/null; then
    echo -e "${RED}✗ RDKit 未安装${NC}"
    echo -e "${YELLOW}请手动安装依赖包${NC}"
    exit 1
fi

echo -e "${GREEN}✓ 环境验证通过${NC}"
echo ""

# 设置 PYTHONPATH
export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"

# 进入 FlashDock 目录
cd "$FLASHDOCK_DIR"

# 显示启动信息
echo -e "${BLUE}启动信息：${NC}"
echo "  项目根目录: $PROJECT_ROOT"
echo "  FlashDock 目录: $FLASHDOCK_DIR"
echo "  端口: $PORT"
echo "  Python: $(which python)"
echo "  Python 版本: $(python --version)"
echo ""

# 检查端口是否被占用
if command -v netstat &> /dev/null; then
    if netstat -tuln 2>/dev/null | grep -q ":$PORT "; then
        echo -e "${YELLOW}⚠ 警告: 端口 $PORT 已被占用${NC}"
        echo -e "${YELLOW}  如果 FlashDock 已经在运行，请先停止它${NC}"
        echo ""
    fi
fi

# 启动 Streamlit
echo -e "${GREEN}正在启动 FlashDock...${NC}"
echo -e "${BLUE}访问地址: http://localhost:$PORT${NC}"
echo ""
echo -e "${YELLOW}按 Ctrl+C 停止服务${NC}"
echo "=========================================="
echo ""

# 运行 Streamlit
streamlit run FlashDock.py --server.port "$PORT" --server.address 0.0.0.0

