@echo off
chcp 65001 >nul 2>&1
cd /d "%~dp0"

echo ========================================
echo Restarting COMPASS Microservices
echo ========================================
echo.

REM Get project root (absolute path)
cd /d %~dp0
set "PROJECT_ROOT=%CD%"

REM Find conda python executables (will be set later)
set PYTHON_AIDDTRAIN=python

REM FlashDock 使用 WSL 环境
set "WSL_DISTRO=Ubuntu-24.04"
set "WSL_ENV_NAME=flash_dock"
set "WSL_PROJECT_ROOT=/mnt/e/Qinchaojun/AIDD-TRAIN"
set "WSL_FLASHDOCK_DIR=%WSL_PROJECT_ROOT%/FLASH_DOCK-main"
set "WSL_PORT=8501"
set "WSL_FLASHDOCK_READY=0"

REM 检查 WSL 是否可用
wsl --status >nul 2>&1
if not errorlevel 1 (
    wsl -d %WSL_DISTRO% -- echo "WSL 发行版可用" >nul 2>&1
    if not errorlevel 1 (
        wsl -d %WSL_DISTRO% bash -c "conda env list | grep -q '^%WSL_ENV_NAME% '" >nul 2>&1
        if not errorlevel 1 (
            set "WSL_FLASHDOCK_READY=1"
        )
    )
)

echo [Step 1/5] Stopping existing services...
echo Finding and stopping processes on ports 8500, 8080, and 8501...

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8500" ^| findstr "LISTENING"') do (
    echo Stopping process %%a (port 8500)
    taskkill /F /PID %%a >nul 2>&1
)

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8080" ^| findstr "LISTENING"') do (
    echo Stopping process %%a (port 8080)
    taskkill /F /PID %%a >nul 2>&1
)

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8501" ^| findstr "LISTENING"') do (
    echo Stopping process %%a (port 8501)
    taskkill /F /PID %%a >nul 2>&1
)

echo Waiting for services to stop completely...
timeout /t 3 /nobreak >nul 2>&1

echo.
echo [Step 2/5] Starting Service Registry (port 8500)...
echo Project Root: %PROJECT_ROOT%
echo.

REM Set PYTHONPATH
set "PYTHONPATH=%PROJECT_ROOT%"

REM Find conda python executables
set PYTHON_AIDDTRAIN=python

if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
    set PYTHON_AIDDTRAIN=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe
) else (
    if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
        set PYTHON_AIDDTRAIN=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe
    )
)

if exist "%PYTHON_AIDDTRAIN%" (
    start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_AIDDTRAIN% services\registry\server.py --host 0.0.0.0 --port 8500"
) else (
    start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python services\registry\server.py --host 0.0.0.0 --port 8500"
)

echo Waiting for Service Registry to start...
timeout /t 5 /nobreak >nul 2>&1

echo.
echo [Step 3/5] Starting COMPASS Service (port 8080)...
if exist "%PYTHON_AIDDTRAIN%" (
    start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_AIDDTRAIN% compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
) else (
    start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
)

echo Waiting for COMPASS service to register...
timeout /t 5 /nobreak >nul 2>&1

echo.
echo [Step 4/5] Starting FLASH-DOCK (port 8501)...
if !WSL_FLASHDOCK_READY! EQU 1 (
    echo Using WSL environment: %WSL_ENV_NAME%
    start "FLASH-DOCK (WSL) - Port 8501" cmd /k "wsl -d %WSL_DISTRO% bash -c \"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %WSL_ENV_NAME% && export PYTHONPATH=%WSL_PROJECT_ROOT% && cd %WSL_FLASHDOCK_DIR% && streamlit run FlashDock.py --server.port %WSL_PORT% --server.address 0.0.0.0\""
) else (
    echo ========================================
    echo [ERROR] WSL FlashDock 环境不可用
    echo ========================================
    echo FlashDock 必须使用 WSL 环境运行
    echo.
    echo 请执行以下步骤设置 WSL 环境：
    echo 1. 确保 WSL 已安装并启用
    echo 2. 确保 WSL 发行版 '%WSL_DISTRO%' 存在
    echo 3. 在 WSL 中手动创建 conda 环境并安装依赖
    echo.
    echo [ERROR] FlashDock 服务启动失败
    echo.
)

echo.
echo [Step 5/5] All services restarted
echo ========================================
echo Service Registry: http://localhost:8500
echo COMPASS Service: http://localhost:8080
echo FLASH-DOCK: http://localhost:8501
echo.
pause
