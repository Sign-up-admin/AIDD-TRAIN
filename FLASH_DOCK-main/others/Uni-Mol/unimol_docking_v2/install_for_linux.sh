#!/bin/bash
# 专为 Linux 系统设计的 Uni-Core 安装脚本

echo "=========================================="
echo "Uni-Core 安装脚本 (Linux)"
echo "=========================================="

# 检查是否在 conda 环境中
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "警告: 未检测到 conda 环境，请先运行:"
    echo "  conda activate flash_dock"
    exit 1
fi

echo "当前 conda 环境: $CONDA_DEFAULT_ENV"
echo ""

# 检查 Git
echo "检查 Git..."
if ! command -v git &> /dev/null; then
    echo "错误: Git 未安装，请先安装 Git"
    echo "Ubuntu/Debian: sudo apt-get install git"
    echo "CentOS/RHEL: sudo yum install git"
    exit 1
fi
echo "Git 版本: $(git --version)"
echo ""

# 方法 1: 直接安装
echo "方法 1: 从 GitHub 直接安装..."
echo "命令: pip install git+https://github.com/dptech-corp/Uni-Core.git@stable"

if pip install git+https://github.com/dptech-corp/Uni-Core.git@stable; then
    echo ""
    echo "✓ 安装成功！"
else
    echo ""
    echo "✗ 方法 1 失败，尝试方法 2..."

    # 方法 2: 克隆后安装
    echo ""
    echo "方法 2: 克隆后安装"
    TEMP_DIR=$(mktemp -d)
    echo "临时目录: $TEMP_DIR"

    cd "$TEMP_DIR"
    if git clone https://github.com/dptech-corp/Uni-Core.git; then
        cd Uni-Core
        git checkout stable 2>/dev/null || echo "使用默认分支"
        if pip install -e .; then
            echo ""
            echo "✓ 安装成功！"
        else
            echo ""
            echo "✗ 安装失败"
            echo "请手动运行:"
            echo "  cd $TEMP_DIR/Uni-Core"
            echo "  pip install -e ."
        fi
    else
        echo ""
        echo "✗ Git 克隆失败，请检查网络连接"
    fi
    cd - > /dev/null
    rm -rf "$TEMP_DIR"
fi

# 验证安装
echo ""
echo "验证安装..."
if python -c "import unicore; from unicore import checkpoint_utils, distributed_utils, options, utils; print('✓ Uni-Core 安装成功！')" 2>/dev/null; then
    echo ""
    echo "=========================================="
    echo "✓ 安装完成！现在可以运行 Uni-Mol Docking V2"
    echo "=========================================="
else
    echo ""
    echo "✗ 验证失败"
    echo "请检查安装过程或查看详细日志"
    echo ""
    echo "手动验证命令:"
    echo "  python -c \"import unicore; print('OK')\""
fi
