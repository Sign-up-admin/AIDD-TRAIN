@echo off
chcp 65001 >nul
echo ========================================
echo 重新启动 COMPASS 微服务系统
echo ========================================
echo.

REM 获取绝对路径
cd /d %~dp0
set PROJECT_ROOT=%CD%

REM 检查conda环境
set PYTHON_CMD=python
if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
    set PYTHON_CMD=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe
) else (
    if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
        set PYTHON_CMD=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe
    )
)

REM 检查FLASH-DOCK的conda环境
set FLASHDOCK_PYTHON=python
if exist "C:\ProgramData\Anaconda3\envs\flash_dock\python.exe" (
    set FLASHDOCK_PYTHON=C:\ProgramData\Anaconda3\envs\flash_dock\python.exe
) else (
    if exist "%USERPROFILE%\anaconda3\envs\flash_dock\python.exe" (
        set FLASHDOCK_PYTHON=%USERPROFILE%\anaconda3\envs\flash_dock\python.exe
    )
)

echo [步骤 1/5] 停止现有服务...
echo 正在查找并停止占用端口 8500、8080 和 8501 的进程...

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8500" ^| findstr "LISTENING"') do (
    echo 停止进程 %%a (端口 8500)
    taskkill /F /PID %%a >nul 2>&1
)

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8080" ^| findstr "LISTENING"') do (
    echo 停止进程 %%a (端口 8080)
    taskkill /F /PID %%a >nul 2>&1
)

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8501" ^| findstr "LISTENING"') do (
    echo 停止进程 %%a (端口 8501)
    taskkill /F /PID %%a >nul 2>&1
)

echo 等待服务完全停止...
timeout /t 3 /nobreak >nul

echo.
echo [步骤 2/4] 启动服务注册中心 (端口 8500)...
echo 项目目录: %PROJECT_ROOT%
echo Python: %PYTHON_CMD%
echo.

REM 设置PYTHONPATH - 只使用绝对路径
set PYTHONPATH=%PROJECT_ROOT%

start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_CMD% services\registry\server.py --host 0.0.0.0 --port 8500"

echo 等待服务注册中心启动...
timeout /t 5 /nobreak >nul

echo.
echo [步骤 3/4] 启动 COMPASS 服务 (端口 8080)...
start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_CMD% compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"

echo 等待 COMPASS 服务启动...
timeout /t 5 /nobreak >nul

echo.
echo [步骤 4/5] 启动 FLASH-DOCK (端口 8501)...
if exist "%FLASHDOCK_PYTHON%" (
    start "FLASH-DOCK - Port 8501" cmd /k "cd /d %PROJECT_ROOT%\FLASH_DOCK-main && set PYTHONPATH=%PROJECT_ROOT% && %FLASHDOCK_PYTHON% -m streamlit run FlashDock.py --server.port 8501"
) else (
    echo [WARNING] FLASH-DOCK Python not found, using default python
    start "FLASH-DOCK - Port 8501" cmd /k "cd /d %PROJECT_ROOT%\FLASH_DOCK-main && set PYTHONPATH=%PROJECT_ROOT% && python -m streamlit run FlashDock.py --server.port 8501"
)

echo 等待 FLASH-DOCK 启动...
timeout /t 5 /nobreak >nul

echo.
echo [步骤 5/5] 检查服务状态...
timeout /t 3 /nobreak >nul

echo.
echo ========================================
echo 服务重新启动完成！
echo ========================================
echo.
echo 服务地址：
echo   - 服务注册中心: http://localhost:8500
echo   - COMPASS 服务: http://localhost:8080
echo   - FLASH-DOCK: http://localhost:8501
echo   - API 文档: http://localhost:8080/docs
echo   - 健康检查: http://localhost:8080/health
echo.
echo 请在新的服务窗口中查看启动日志
echo 如果服务未正常启动，请检查窗口中的错误信息
echo.
pause

