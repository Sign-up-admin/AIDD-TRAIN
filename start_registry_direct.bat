@echo off
chcp 65001 >nul
cd /d %~dp0
echo 启动服务注册中心...
echo 端口: 8500
echo 工作目录: %CD%
echo.
set PYTHONPATH=%CD%;%PYTHONPATH%
C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe services/registry/server.py --host 0.0.0.0 --port 8500
pause
