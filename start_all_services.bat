@echo off
setlocal enabledelayedexpansion
chcp 65001 >nul 2>&1
cd /d %~dp0

echo ========================================
echo Starting COMPASS Microservices
echo ========================================
echo.

REM Find conda and python executables
set PYTHON_AIDDTRAIN=
set PYTHON_FLASHDOCK=

REM Try to find conda environments (check D:\conda_envs first as recommended)
if exist "D:\conda_envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=D:\conda_envs\AIDDTRAIN\python.exe"
    echo Found AIDDTRAIN at D:\conda_envs
) else if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe"
    echo Found AIDDTRAIN at C:\ProgramData\Anaconda3
) else if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe"
    echo Found AIDDTRAIN at %USERPROFILE%\anaconda3
) else (
    echo [WARNING] AIDDTRAIN environment not found, using default python
    set "PYTHON_AIDDTRAIN=python"
)

REM 注意：只有 FlashDock 使用 WSL 环境，其他服务（服务注册中心、COMPASS）使用 Windows 环境
REM FlashDock 使用 WSL 环境
REM 配置 WSL 相关参数
set "WSL_DISTRO=Ubuntu-24.04"
set "WSL_ENV_NAME=flash_dock"
set "WSL_PROJECT_ROOT=/mnt/e/Qinchaojun/AIDD-TRAIN"
set "WSL_FLASHDOCK_DIR=%WSL_PROJECT_ROOT%/FLASH_DOCK-main"
set "WSL_PORT=8501"
set "WSL_FLASHDOCK_READY=0"

REM 检查 WSL 是否可用
wsl --status >nul 2>&1
if errorlevel 1 (
    echo [WARNING] WSL 未安装或未启用，FlashDock 将无法在 WSL 中运行
    set "WSL_FLASHDOCK_READY=0"
) else (
    REM 检查 WSL 发行版是否存在
    wsl -d %WSL_DISTRO% -- echo "WSL 发行版可用" >nul 2>&1
    if errorlevel 1 (
        echo [WARNING] WSL 发行版 '%WSL_DISTRO%' 不存在
        echo 请检查 WSL 发行版名称，或修改脚本中的 WSL_DISTRO 变量
        set "WSL_FLASHDOCK_READY=0"
    ) else (
        REM 检查 conda 环境是否存在
        wsl -d %WSL_DISTRO% bash -c "conda env list | grep -q '^%WSL_ENV_NAME% '" >nul 2>&1
        if errorlevel 1 (
            echo [WARNING] WSL Conda 环境 '%WSL_ENV_NAME%' 不存在
            echo 请手动创建 conda 环境并安装依赖
            set "WSL_FLASHDOCK_READY=0"
        ) else (
            echo [OK] WSL FlashDock 环境已就绪
            set "WSL_FLASHDOCK_READY=1"
        )
    )
)

REM Get project root directory (absolute path)
cd /d %~dp0
set "PROJECT_ROOT=%CD%"

echo Project Root: %PROJECT_ROOT%
echo.

REM 注意：只有 FLASH-DOCK 使用 WSL 环境，其他服务使用 Windows 环境
echo Step 1: Starting Service Registry...
echo Using Python: !PYTHON_AIDDTRAIN!
if not "!PYTHON_AIDDTRAIN!"=="python" (
    if exist "!PYTHON_AIDDTRAIN!" (
        start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && \"!PYTHON_AIDDTRAIN!\" services\registry\server.py --host 0.0.0.0 --port 8500 && pause"
    ) else (
        echo [ERROR] AIDDTRAIN Python executable not found: !PYTHON_AIDDTRAIN!
        start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python services\registry\server.py --host 0.0.0.0 --port 8500 && pause"
    )
) else (
    start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python services\registry\server.py --host 0.0.0.0 --port 8500 && pause"
)
echo Waiting for Service Registry to start...
timeout /t 5 /nobreak >nul 2>&1

echo Step 2: Starting COMPASS Service...
echo Using Python: !PYTHON_AIDDTRAIN!
if not "!PYTHON_AIDDTRAIN!"=="python" (
    if exist "!PYTHON_AIDDTRAIN!" (
        start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && \"!PYTHON_AIDDTRAIN!\" compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500 && pause"
    ) else (
        echo [ERROR] AIDDTRAIN Python executable not found: !PYTHON_AIDDTRAIN!
        echo COMPASS service requires AIDDTRAIN environment with rdkit installed.
        start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500 && pause"
    )
) else (
    echo [WARNING] Using default Python for COMPASS (may fail without rdkit)
    start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500 && pause"
)
echo Waiting for COMPASS service to register...
timeout /t 5 /nobreak >nul 2>&1

echo Step 3: Starting FLASH-DOCK...
if !WSL_FLASHDOCK_READY! EQU 1 (
    echo Using WSL environment: %WSL_ENV_NAME%
    echo WSL Distribution: %WSL_DISTRO%
    echo Project Path: %WSL_PROJECT_ROOT%
    echo.
    start "FLASH-DOCK (WSL) - Port 8501" cmd /k "wsl -d %WSL_DISTRO% bash -c \"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %WSL_ENV_NAME% && export PYTHONPATH=%WSL_PROJECT_ROOT% && cd %WSL_FLASHDOCK_DIR% && streamlit run FlashDock.py --server.port %WSL_PORT% --server.address 0.0.0.0\" && pause"
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
    pause
    exit /b 1
)

echo.
echo ========================================
echo All services are starting in separate windows
echo ========================================
echo Service Registry: http://localhost:8500
echo COMPASS Service: http://localhost:8080
echo FLASH-DOCK: http://localhost:8501
echo.
echo Please check the opened windows for service status
echo If services fail to start, check the error messages in those windows
echo.
pause
