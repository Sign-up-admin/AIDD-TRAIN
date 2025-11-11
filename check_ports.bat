@echo off
chcp 65001 >nul
echo ========================================
echo COMPASS服务系统端口检查工具
echo ========================================
echo.

cd /d %~dp0

set PYTHONPATH=%~dp0;%PYTHONPATH%

REM Check if Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo [错误] 未找到Python，请先安装Python
    pause
    exit /b 1
)

echo 正在检查端口...
echo.

python check_ports.py

echo.
pause






