#!/bin/bash
# WSL2 运行脚本
# 用于在 WSL2 环境中运行 Uni-Mol 对接预测程序

set -e  # 遇到错误立即退出

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Uni-Mol Docking WSL2 运行脚本 ===${NC}"

# 检查是否在 WSL2 环境中
if [ ! -f /proc/version ] || ! grep -qi "microsoft\|wsl" /proc/version; then
    echo -e "${YELLOW}警告: 可能不在 WSL2 环境中运行${NC}"
fi

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# 检查 Python 环境
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}错误: 未找到 python3${NC}"
    exit 1
fi

# 检查 conda 环境
if [ -n "$CONDA_DEFAULT_ENV" ]; then
    echo -e "${GREEN}使用 Conda 环境: $CONDA_DEFAULT_ENV${NC}"
else
    echo -e "${YELLOW}警告: 未检测到 Conda 环境，建议激活 flash_dock 环境${NC}"
    echo -e "${YELLOW}运行: conda activate flash_dock${NC}"
fi

# 转换 Windows 路径参数为 WSL2 路径
convert_path() {
    local path="$1"
    if [[ "$path" =~ ^[A-Za-z]: ]]; then
        # Windows 路径格式: E:\path\to\file
        drive_letter=$(echo "$path" | cut -d: -f1 | tr '[:upper:]' '[:lower:]')
        rest_path=$(echo "$path" | cut -d: -f2 | tr '\\' '/')
        echo "/mnt/$drive_letter$rest_path"
    else
        echo "$path"
    fi
}

# 转换所有路径参数
CONVERTED_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-protein|--input-ligand|--input-docking-grid|--output-ligand-dir|--model-dir|--input-batch-file)
            CONVERTED_ARGS+=("$1")
            shift
            if [[ $# -gt 0 ]]; then
                CONVERTED_ARGS+=("$(convert_path "$1")")
            fi
            shift
            ;;
        *)
            CONVERTED_ARGS+=("$1")
            shift
            ;;
    esac
done

# 运行 demo.py
echo -e "${GREEN}运行命令: python3 demo.py ${CONVERTED_ARGS[*]}${NC}"
echo ""

python3 demo.py "${CONVERTED_ARGS[@]}"

exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo -e "${GREEN}=== 运行成功 ===${NC}"
else
    echo -e "${RED}=== 运行失败，退出码: $exit_code ===${NC}"
fi

exit $exit_code



