#!/bin/bash
# WSL Java环境诊断脚本
# 用于检查WSL中Java环境配置

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "============================================================"
echo "Java环境诊断 (WSL版本)"
echo "============================================================"
echo

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "工作目录: $SCRIPT_DIR"
echo

all_ok=true

# 检查是否在WSL中
if [ -f /proc/version ] && grep -qi "microsoft\|wsl" /proc/version; then
    echo -e "${GREEN}✓ 检测到WSL环境${NC}"
else
    echo -e "${YELLOW}⚠ 警告: 可能不在WSL环境中运行${NC}"
fi
echo

# 检查java命令
if command -v java &> /dev/null; then
    java_path=$(which java)
    echo -e "${GREEN}✓ Java可执行文件找到: $java_path${NC}"
    
    # 获取Java版本
    java_version_output=$(java -version 2>&1)
    echo "Java版本信息:"
    echo "$java_version_output" | head -n 3
    
    # 检查版本号
    if echo "$java_version_output" | grep -q "version"; then
        version=$(echo "$java_version_output" | grep -oP 'version "(\d+)' | grep -oP '\d+' | head -n 1)
        if [ ! -z "$version" ]; then
            if [ "$version" -ge 17 ] && [ "$version" -le 23 ]; then
                echo -e "${GREEN}✓ Java版本 $version 符合要求（17-23）${NC}"
            else
                echo -e "${RED}✗ Java版本 $version 不符合要求（需要17-23）${NC}"
                all_ok=false
            fi
        fi
    fi
else
    echo -e "${RED}✗ Java可执行文件未找到${NC}"
    all_ok=false
fi
echo

# 检查JAVA_HOME
if [ ! -z "$JAVA_HOME" ]; then
    echo -e "${GREEN}✓ JAVA_HOME: $JAVA_HOME${NC}"
    if [ -f "$JAVA_HOME/bin/java" ]; then
        echo -e "${GREEN}✓ JAVA_HOME/bin/java 存在${NC}"
    else
        echo -e "${RED}✗ JAVA_HOME/bin/java 不存在${NC}"
        all_ok=false
    fi
else
    echo -e "${YELLOW}⚠ JAVA_HOME环境变量未设置（可选）${NC}"
fi
echo

# 检查P2Rank
p2rank_home="$SCRIPT_DIR/others/p2rank_2.5"
if [ -d "$p2rank_home" ]; then
    echo -e "${GREEN}✓ P2Rank目录存在: $p2rank_home${NC}"
    
    if [ -f "$p2rank_home/prank" ]; then
        echo -e "${GREEN}✓ prank 脚本存在${NC}"
        # 检查脚本是否有执行权限
        if [ -x "$p2rank_home/prank" ]; then
            echo -e "${GREEN}✓ prank 脚本有执行权限${NC}"
        else
            echo -e "${YELLOW}⚠ prank 脚本没有执行权限，尝试添加...${NC}"
            chmod +x "$p2rank_home/prank" 2>/dev/null && echo -e "${GREEN}✓ 已添加执行权限${NC}" || echo -e "${RED}✗ 无法添加执行权限${NC}"
        fi
    else
        echo -e "${RED}✗ prank 脚本不存在${NC}"
        all_ok=false
    fi
    
    if [ -f "$p2rank_home/bin/p2rank.jar" ]; then
        echo -e "${GREEN}✓ p2rank.jar 存在${NC}"
    else
        echo -e "${RED}✗ p2rank.jar 不存在${NC}"
        all_ok=false
    fi
else
    echo -e "${RED}✗ P2Rank目录不存在: $p2rank_home${NC}"
    all_ok=false
fi
echo

# 检查示例文件
example_file="$SCRIPT_DIR/examples/pocket/protein.pdb"
if [ -f "$example_file" ]; then
    echo -e "${GREEN}✓ 示例文件存在: $example_file${NC}"
    file_size=$(stat -f%z "$example_file" 2>/dev/null || stat -c%s "$example_file" 2>/dev/null)
    echo "  文件大小: $file_size 字节"
else
    echo -e "${YELLOW}⚠ 示例文件不存在: $example_file${NC}"
    # 示例文件不存在不影响Java环境检查
fi
echo

echo "============================================================"
if [ "$all_ok" = true ]; then
    echo -e "${GREEN}✓ 所有检查通过！Java环境配置正确。${NC}"
    exit 0
else
    echo -e "${RED}✗ 发现问题，请参考《Java安装指南-WSL版本.md》进行修复。${NC}"
    exit 1
fi

