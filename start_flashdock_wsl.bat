@echo off
chcp 65001 >nul 2>&1
REM FlashDock WSL 启动脚本（从 Windows 调用）
REM 这个脚本从 Windows 调用 WSL 中的 conda 环境来运行 FlashDock

echo ============================================
echo FlashDock WSL 启动脚本（从 Windows 调用）
echo ============================================
echo.

REM 配置
set WSL_DISTRO=Ubuntu-24.04
set ENV_NAME=flash_dock_wsl
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

echo [1/4] 检查 WSL 发行版: %WSL_DISTRO%
wsl -d %WSL_DISTRO% -- echo "WSL 发行版可用" >nul 2>&1
if errorlevel 1 (
    echo [错误] WSL 发行版 '%WSL_DISTRO%' 不存在
    echo 请检查 WSL 发行版名称，或修改脚本中的 WSL_DISTRO 变量
    pause
    exit /b 1
)
echo [OK] WSL 发行版可用
echo.

echo [2/4] 检查 conda 环境: %ENV_NAME%
wsl -d %WSL_DISTRO% bash -c "conda env list | grep -q '^%ENV_NAME% '" >nul 2>&1
if errorlevel 1 (
    echo [错误] Conda 环境 '%ENV_NAME%' 不存在
    echo 请手动创建 conda 环境并安装依赖
    pause
    exit /b 1
)
echo [OK] Conda 环境存在
echo.

echo [3/4] 验证关键依赖
wsl -d %WSL_DISTRO% bash -c "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && python -c 'import streamlit; import rdkit' 2>/dev/null" >nul 2>&1
if errorlevel 1 (
    echo [警告] 依赖验证失败，尝试重新安装...
    wsl -d %WSL_DISTRO% bash -c "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && pip install streamlit rdkit-pypi==2022.9.3 -q"
)
echo [OK] 依赖验证通过
echo.

echo [4/4] 启动 FlashDock
echo.
echo 项目路径: %PROJECT_ROOT%
echo 端口: %PORT%
echo 访问地址: http://localhost:%PORT%
echo.
echo 按 Ctrl+C 停止服务
echo ============================================
echo.

REM 启动 FlashDock（在 WSL 中运行）
wsl -d %WSL_DISTRO% bash -c "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && export PYTHONPATH=%PROJECT_ROOT% && cd %FLASHDOCK_DIR% && streamlit run FlashDock.py --server.port %PORT% --server.address 0.0.0.0"

echo.
echo ============================================
echo FlashDock 已停止
echo ============================================
pause

