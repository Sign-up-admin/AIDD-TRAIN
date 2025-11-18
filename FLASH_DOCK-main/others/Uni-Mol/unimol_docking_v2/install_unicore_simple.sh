#!/bin/bash
# 简单的 Uni-Core 安装脚本 - 一键安装

echo "=========================================="
echo "正在安装 Uni-Core..."
echo "=========================================="

# 检查是否在 conda 环境中
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "警告: 未检测到激活的 conda 环境"
    echo "请先运行: conda activate flash_dock"
    exit 1
fi

echo "当前环境: $CONDA_DEFAULT_ENV"
echo ""

# 方法 1: 直接安装
echo "尝试方法 1: 从 GitHub 安装..."
if pip install git+https://github.com/dptech-corp/Uni-Core.git@stable; then
    echo ""
    echo "✓ 安装成功！"
    echo ""
    echo "验证安装..."
    if python -c "import unicore; from unicore import checkpoint_utils, distributed_utils, options, utils; print('✓ 验证成功！')" 2>/dev/null; then
        echo ""
        echo "=========================================="
        echo "✓ Uni-Core 已成功安装并验证！"
        echo "=========================================="
        exit 0
    fi
fi

# 方法 2: 克隆后安装
echo ""
echo "方法 1 失败，尝试方法 2: 克隆后安装..."
TEMP_DIR=$(mktemp -d)
cd "$TEMP_DIR"

if git clone https://github.com/dptech-corp/Uni-Core.git; then
    cd Uni-Core
    git checkout stable 2>/dev/null || echo "使用默认分支"
    if pip install -e .; then
        echo ""
        echo "✓ 安装成功！"
        cd -
        rm -rf "$TEMP_DIR"
        
        # 验证
        if python -c "import unicore; from unicore import checkpoint_utils, distributed_utils, options, utils; print('✓ 验证成功！')" 2>/dev/null; then
            echo ""
            echo "=========================================="
            echo "✓ Uni-Core 已成功安装并验证！"
            echo "=========================================="
            exit 0
        fi
    fi
    cd -
    rm -rf "$TEMP_DIR"
fi

echo ""
echo "✗ 安装失败，请检查网络连接和 Git 配置"
echo "详细说明请查看: INSTALL_UNICORE.md"
exit 1

