@echo off
chcp 65001 >nul 2>&1
REM FlashDock WSL 启动脚本（修复版）
REM 这个脚本从 Windows 调用 WSL 中的 conda 环境来运行 FlashDock

echo ============================================
echo FlashDock WSL 启动脚本（修复版）
echo ============================================
echo.

REM 配置
set WSL_DISTRO=Ubuntu-24.04
set ENV_NAME=flash_dock
set PROJECT_ROOT=/mnt/e/Qinchaojun/AIDD-TRAIN
set FLASHDOCK_DIR=%PROJECT_ROOT%/FLASH_DOCK-main
set PORT=8501

REM 检查 WSL 是否可用
wsl --status >nul 2>&1
if errorlevel 1 (
    echo [错误] WSL 未安装或未启用
    echo 请先安装并启用 WSL2
    pause
    exit /b 1
)

echo [1/5] 检查 WSL 发行版: %WSL_DISTRO%
wsl -d %WSL_DISTRO% -- echo "WSL 发行版可用" >nul 2>&1
if errorlevel 1 (
    echo [错误] WSL 发行版 '%WSL_DISTRO%' 不存在
    echo 请检查 WSL 发行版名称，或修改脚本中的 WSL_DISTRO 变量
    pause
    exit /b 1
)
echo [OK] WSL 发行版可用
echo.

echo [2/5] 检查 conda 环境: %ENV_NAME%
set "CHECK_ENV_CMD=source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda env list | grep -q '^%ENV_NAME% '"
wsl -d %WSL_DISTRO% bash -c "%CHECK_ENV_CMD%" >nul 2>&1
if errorlevel 1 (
    echo [错误] Conda 环境 '%ENV_NAME%' 不存在
    echo.
    echo 可用的环境:
    set "LIST_ENV_CMD=source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda env list"
    wsl -d %WSL_DISTRO% bash -c "%LIST_ENV_CMD%"
    echo.
    echo 请手动创建 conda 环境并安装依赖
    pause
    exit /b 1
)
echo [OK] Conda 环境存在
echo.

echo [3/5] 验证关键依赖
set "CHECK_STREAMLIT_CMD=source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && python -c 'import streamlit; print(\"streamlit OK\")' 2>&1"
wsl -d %WSL_DISTRO% bash -c "%CHECK_STREAMLIT_CMD%" | findstr /C:"streamlit OK" >nul 2>&1
if errorlevel 1 (
    echo [警告] Streamlit 可能未正确安装
    echo 尝试重新安装...
    set "INSTALL_STREAMLIT_CMD=source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && pip install streamlit -q"
    wsl -d %WSL_DISTRO% bash -c "%INSTALL_STREAMLIT_CMD%"
) else (
    echo [OK] Streamlit 已安装
)
echo.

echo [4/5] 检查项目文件
set "CHECK_FILE_CMD=test -f '%FLASHDOCK_DIR%/FlashDock.py' && echo 'FILE_EXISTS' || echo 'FILE_NOT_FOUND'"
wsl -d %WSL_DISTRO% bash -c "%CHECK_FILE_CMD%" | findstr /C:"FILE_EXISTS" >nul 2>&1
if errorlevel 1 (
    echo [错误] FlashDock.py 文件不存在: %FLASHDOCK_DIR%/FlashDock.py
    pause
    exit /b 1
)
echo [OK] 项目文件存在
echo.

echo [5/5] 检查端口占用
set "CHECK_PORT_CMD=netstat -tuln 2>/dev/null | grep -q ':%PORT% ' && echo 'IN_USE' || echo 'AVAILABLE'"
wsl -d %WSL_DISTRO% bash -c "%CHECK_PORT_CMD%" | findstr /C:"IN_USE" >nul 2>&1
if not errorlevel 1 (
    echo [警告] 端口 %PORT% 已被占用
    echo 正在尝试停止占用进程...
    set "KILL_PORT_CMD=lsof -ti :%PORT% | xargs kill -9 2>/dev/null || true"
    wsl -d %WSL_DISTRO% bash -c "%KILL_PORT_CMD%"
    timeout /t 2 /nobreak >nul 2>&1
)
echo [OK] 端口检查完成
echo.

echo ============================================
echo 启动 FlashDock
echo ============================================
echo.
echo 项目路径: %PROJECT_ROOT%
echo 端口: %PORT%
echo 访问地址: http://localhost:%PORT%
echo.
echo 提示: 服务将在新窗口中启动
echo 按 Ctrl+C 停止服务
echo ============================================
echo.

REM 启动 FlashDock（在新窗口中运行）
REM 构建WSL命令（避免引号嵌套问题）
set "WSL_CMD=source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && export PYTHONPATH=%PROJECT_ROOT% && cd %FLASHDOCK_DIR% && streamlit run FlashDock.py --server.port %PORT% --server.address 0.0.0.0"
start "FLASH-DOCK (WSL) - Port %PORT%" cmd /k wsl -d %WSL_DISTRO% bash -c "%WSL_CMD%"

echo [OK] FlashDock 已在新窗口中启动
echo.
echo 等待服务启动（约10-15秒）...
timeout /t 10 /nobreak >nul 2>&1

REM 检查服务是否启动
echo.
echo 检查服务状态...
powershell -Command "try { $response = Invoke-WebRequest -Uri 'http://localhost:%PORT%' -TimeoutSec 3 -UseBasicParsing; Write-Host '[OK] FLASH-DOCK 服务已启动' } catch { Write-Host '[WARNING] 服务可能还在启动中，请检查新窗口' }"

echo.
echo ============================================
echo 如果服务未启动，请检查新窗口中的错误信息
echo ============================================
pause

