@echo off
chcp 65001 >nul
echo ========================================
echo 重新启动 COMPASS 微服务系统
echo ========================================
echo.

echo [步骤 1/4] 停止现有服务...
echo 正在查找并停止占用端口 8500 和 8080 的进程...

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8500" ^| findstr "LISTENING"') do (
    echo 停止进程 %%a (端口 8500)
    taskkill /F /PID %%a >nul 2>&1
)

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8080" ^| findstr "LISTENING"') do (
    echo 停止进程 %%a (端口 8080)
    taskkill /F /PID %%a >nul 2>&1
)

echo 等待服务完全停止...
timeout /t 2 /nobreak >nul

echo.
echo [步骤 2/4] 启动服务注册中心 (端口 8500)...
cd /d %~dp0
set PYTHONPATH=%CD%;%PYTHONPATH%
start "Service Registry - Port 8500" cmd /k "cd /d %~dp0 && set PYTHONPATH=%CD%;%PYTHONPATH% && C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe services/registry/server.py --host 0.0.0.0 --port 8500"
echo 等待服务注册中心启动...
timeout /t 5 /nobreak >nul

echo [步骤 3/4] 启动 COMPASS 服务 (端口 8080)...
start "COMPASS Service - Port 8080" cmd /k "cd /d %~dp0 && set PYTHONPATH=%CD%;%PYTHONPATH% && C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
echo 等待 COMPASS 服务启动...
timeout /t 5 /nobreak >nul

echo.
echo [步骤 4/4] 检查服务状态...
timeout /t 3 /nobreak >nul

echo.
echo ========================================
echo 服务重新启动完成！
echo ========================================
echo.
echo 服务地址：
echo   - 服务注册中心: http://localhost:8500
echo   - COMPASS 服务: http://localhost:8080
echo   - API 文档: http://localhost:8080/docs
echo   - 健康检查: http://localhost:8080/health
echo.
echo 请在新的服务窗口中查看启动日志
echo 如果服务未正常启动，请检查窗口中的错误信息
echo.
pause

