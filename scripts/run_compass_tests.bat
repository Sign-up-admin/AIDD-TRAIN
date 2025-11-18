@echo off
REM COMPASS自动化测试快速运行脚本 (Windows)

setlocal enabledelayedexpansion

echo ========================================
echo COMPASS自动化测试
echo ========================================
echo.

REM 获取项目根目录
set "PROJECT_ROOT=%~dp0.."
cd /d "%PROJECT_ROOT%"

REM 检查Python
python --version >nul 2>&1
if errorlevel 1 (
    echo 错误: 未找到Python
    exit /b 1
)

REM 检查pytest
python -m pytest --version >nul 2>&1
if errorlevel 1 (
    echo 警告: pytest未安装，正在安装...
    pip install pytest pytest-cov
)

REM 解析参数
set "TEST_TYPE=%~1"
if "%TEST_TYPE%"=="" set "TEST_TYPE=all"

if "%TEST_TYPE%"=="unit" (
    echo 运行单元测试...
    python -m pytest tests/ -m unit -v --cov=compass --cov-report=term-missing
) else if "%TEST_TYPE%"=="integration" (
    echo 运行集成测试...
    python -m pytest tests/ -m integration -v --cov=compass --cov-report=term-missing
) else if "%TEST_TYPE%"=="e2e" (
    echo 运行端到端测试...
    python -m pytest tests/ -m e2e -v
) else (
    echo 运行所有测试...
    python -m pytest tests/ -v --cov=compass --cov-report=term-missing --cov-report=html
)

echo.
echo 测试完成！

endlocal

