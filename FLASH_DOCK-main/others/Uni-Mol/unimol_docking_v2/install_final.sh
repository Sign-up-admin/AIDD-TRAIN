#!/bin/bash
# 最终 Uni-Core 安装脚本 - 一步到位

echo "=========================================="
echo "Uni-Core 最终安装脚本"
echo "=========================================="
echo ""

# 检查 conda 环境
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "错误: 未检测到 conda 环境"
    echo "请先运行: conda activate flash_dock"
    exit 1
fi

echo "当前 conda 环境: $CONDA_DEFAULT_ENV"
echo "Python 路径: $(which python)"
echo ""

# 检查 Git
echo "检查 Git..."
if ! command -v git &> /dev/null; then
    echo "错误: Git 未安装"
    echo "请安装 Git: sudo apt-get install git"
    exit 1
fi

echo "Git 版本: $(git --version)"
echo ""

# 升级 pip
echo "升级 pip..."
pip install --upgrade pip
echo ""

# 方法 1: 直接安装
echo "=========================================="
echo "方法 1: 直接安装 Uni-Core"
echo "=========================================="
echo "执行命令: pip install git+https://github.com/dptech-corp/Uni-Core.git@stable"

if pip install git+https://github.com/dptech-corp/Uni-Core.git@stable --timeout=300; then
    echo ""
    echo "✓ 方法 1 成功！"
else
    echo ""
    echo "✗ 方法 1 失败，尝试方法 2..."

    # 方法 2: 使用不同协议
    echo ""
    echo "=========================================="
    echo "方法 2: 使用 git 协议安装"
    echo "=========================================="
    echo "执行命令: pip install git+git://github.com/dptech-corp/Uni-Core.git@stable"

    if pip install git+git://github.com/dptech-corp/Uni-Core.git@stable --timeout=300; then
        echo ""
        echo "✓ 方法 2 成功！"
    else
        echo ""
        echo "✗ 方法 2 失败，尝试方法 3..."

        # 方法 3: 克隆后安装
        echo ""
        echo "=========================================="
        echo "方法 3: 克隆后安装"
        echo "=========================================="

        TEMP_DIR="/tmp/unicore_install_$(date +%s)"
        echo "临时目录: $TEMP_DIR"
        mkdir -p "$TEMP_DIR"
        cd "$TEMP_DIR"

        echo "克隆仓库..."
        if git clone https://github.com/dptech-corp/Uni-Core.git --depth=1; then
            cd Uni-Core
            echo "切换到 stable 分支..."
            git checkout stable 2>/dev/null || echo "使用默认分支"

            echo "安装包..."
            if pip install -e .; then
                echo ""
                echo "✓ 方法 3 成功！"
            else
                echo ""
                echo "✗ 安装失败"
                echo "请手动执行:"
                echo "  cd $TEMP_DIR/Uni-Core"
                echo "  pip install -e ."
                cd - >/dev/null
                exit 1
            fi
        else
            echo ""
            echo "✗ 克隆失败，请检查网络连接"
            cd - >/dev/null
            exit 1
        fi

        # 清理临时目录
        cd - >/dev/null
        rm -rf "$TEMP_DIR"
    fi
fi

# 验证安装
echo ""
echo "=========================================="
echo "验证安装"
echo "=========================================="

echo "测试 1: 导入 unicore..."
if python -c "import unicore; print('✓ unicore 导入成功')"; then
    echo ""
    echo "测试 2: 导入核心模块..."
    if python -c "from unicore import checkpoint_utils, distributed_utils, options, utils; print('✓ 核心模块导入成功')"; then
        echo ""
        echo "测试 3: 完整功能测试..."
        if python -c "
import unicore
from unicore import checkpoint_utils, distributed_utils, options, utils, tasks, models
print('✓ Uni-Core 完整安装成功！')
        "; then
            echo ""
            echo "🎉 ========================================="
            echo "🎉   Uni-Core 安装完成！"
            echo "🎉 ========================================="
            echo ""
            echo "现在可以运行 Uni-Mol Docking V2 了！"
            echo ""
            echo "使用方法:"
            echo "  cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface"
            echo "  python demo.py [参数]"
            exit 0
        fi
    fi
fi

echo ""
echo "❌ 安装验证失败"
echo "请检查错误信息并重试"
exit 1
