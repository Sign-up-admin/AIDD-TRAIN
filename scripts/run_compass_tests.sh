#!/bin/bash
# COMPASS自动化测试快速运行脚本

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 项目根目录
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}COMPASS自动化测试${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# 检查Python环境
if ! command -v python &> /dev/null; then
    echo -e "${RED}错误: 未找到Python${NC}"
    exit 1
fi

# 检查pytest
if ! python -m pytest --version &> /dev/null; then
    echo -e "${YELLOW}警告: pytest未安装，正在安装...${NC}"
    pip install pytest pytest-cov
fi

# 解析参数
TEST_TYPE="${1:-all}"
VERBOSE="${2:-true}"

case "$TEST_TYPE" in
    unit)
        echo -e "${GREEN}运行单元测试...${NC}"
        python -m pytest tests/ -m unit -v --cov=compass --cov-report=term-missing
        ;;
    integration)
        echo -e "${GREEN}运行集成测试...${NC}"
        python -m pytest tests/ -m integration -v --cov=compass --cov-report=term-missing
        ;;
    e2e)
        echo -e "${GREEN}运行端到端测试...${NC}"
        python -m pytest tests/ -m e2e -v
        ;;
    all|*)
        echo -e "${GREEN}运行所有测试...${NC}"
        python -m pytest tests/ -v --cov=compass --cov-report=term-missing --cov-report=html
        ;;
esac

echo ""
echo -e "${GREEN}测试完成！${NC}"

