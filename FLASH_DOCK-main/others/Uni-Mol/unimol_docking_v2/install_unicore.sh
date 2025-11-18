#!/bin/bash
# Uni-Core 安装脚本
# 用于解决 ModuleNotFoundError: No module named 'unicore' 错误

set -e  # 遇到错误立即退出

echo "=========================================="
echo "Uni-Core 安装脚本"
echo "=========================================="
echo ""

# 检查是否在 conda 环境中
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "警告: 未检测到激活的 conda 环境"
    echo "请先激活 flash_dock 环境: conda activate flash_dock"
    read -p "是否继续? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
else
    echo "当前 conda 环境: $CONDA_DEFAULT_ENV"
fi

# 检查 Git 是否安装
echo ""
echo "检查 Git..."
if ! command -v git &> /dev/null; then
    echo "错误: Git 未安装，请先安装 Git"
    exit 1
fi
echo "Git 版本: $(git --version)"

# 检查 Python 和 pip
echo ""
echo "检查 Python 环境..."
python --version
pip --version

# 方法 1: 从 GitHub 直接安装（推荐）
echo ""
echo "=========================================="
echo "方法 1: 从 GitHub 直接安装 Uni-Core"
echo "=========================================="
if pip install git+https://github.com/dptech-corp/Uni-Core.git@stable; then
    echo ""
    echo "✓ Uni-Core 安装成功！"
    
    # 验证安装
    echo ""
    echo "验证安装..."
    if python -c "import unicore; print('✓ unicore 模块导入成功')" 2>/dev/null; then
        echo "✓ 核心模块检查..."
        if python -c "from unicore import checkpoint_utils, distributed_utils, options, utils; print('✓ 核心模块导入成功')" 2>/dev/null; then
            echo ""
            echo "=========================================="
            echo "✓ Uni-Core 安装并验证成功！"
            echo "=========================================="
            exit 0
        else
            echo "警告: 核心模块导入失败，但基础模块已安装"
        fi
    else
        echo "警告: unicore 模块导入失败"
    fi
    exit 0
fi

# 方法 2: 如果方法 1 失败，尝试使用 git:// 协议
echo ""
echo "方法 1 失败，尝试方法 2: 使用 git:// 协议"
echo "=========================================="
if pip install "git+git://github.com/dptech-corp/Uni-Core.git@stable"; then
    echo ""
    echo "✓ Uni-Core 安装成功！"
    exit 0
fi

# 方法 3: 克隆后安装
echo ""
echo "方法 2 失败，尝试方法 3: 克隆后安装"
echo "=========================================="
TEMP_DIR=$(mktemp -d)
echo "临时目录: $TEMP_DIR"

cd "$TEMP_DIR"
if git clone https://github.com/dptech-corp/Uni-Core.git; then
    cd Uni-Core
    git checkout stable 2>/dev/null || echo "注意: 无法切换到 stable 分支，使用默认分支"
    
    if pip install -e .; then
        echo ""
        echo "✓ Uni-Core 安装成功！"
        cd -
        rm -rf "$TEMP_DIR"
        exit 0
    else
        echo "错误: pip install -e . 失败"
        cd -
        rm -rf "$TEMP_DIR"
        exit 1
    fi
else
    echo "错误: Git clone 失败"
    echo "可能的原因:"
    echo "  1. 网络连接问题"
    echo "  2. GitHub 访问受限"
    echo "  3. 需要配置代理"
    rm -rf "$TEMP_DIR"
    exit 1
fi

