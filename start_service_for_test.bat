@echo off
chcp 65001 >nul 2>&1
echo ========================================
echo 启动COMPASS服务用于测试
echo ========================================
echo.

cd /d %~dp0
set PYTHONPATH=%CD%

echo 正在启动服务...
echo 服务地址: http://localhost:8080
echo 按 Ctrl+C 停止服务
echo.

python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500

pause


