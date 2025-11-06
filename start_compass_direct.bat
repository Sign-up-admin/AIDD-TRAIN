@echo off
chcp 65001 >nul
cd /d %~dp0
echo 启动 COMPASS 服务...
echo 端口: 8080
echo 注册中心: http://localhost:8500
echo 工作目录: %CD%
echo.
set PYTHONPATH=%CD%;%PYTHONPATH%
C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
pause
