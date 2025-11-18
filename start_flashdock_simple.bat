@echo off
chcp 65001 >nul 2>&1
REM FlashDock WSL 启动脚本（简化版，避免引号嵌套问题）

echo ============================================
echo FlashDock WSL 启动脚本
echo ============================================
echo.

REM 配置
set WSL_DISTRO=Ubuntu-24.04
set ENV_NAME=flash_dock
set PROJECT_ROOT=/mnt/e/Qinchaojun/AIDD-TRAIN
set FLASHDOCK_DIR=%PROJECT_ROOT%/FLASH_DOCK-main
set PORT=8501

REM 检查 WSL
wsl --status >nul 2>&1
if errorlevel 1 (
    echo [错误] WSL 未安装或未启用
    pause
    exit /b 1
)

echo [OK] WSL 可用
echo [OK] 环境: %ENV_NAME%
echo [OK] 端口: %PORT%
echo.

REM 构建启动命令（使用临时变量避免引号问题）
set "WSL_CMD=source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && export PYTHONPATH=%PROJECT_ROOT% && cd %FLASHDOCK_DIR% && streamlit run FlashDock.py --server.port %PORT% --server.address 0.0.0.0"

echo 正在启动 FlashDock...
echo 访问地址: http://localhost:%PORT%
echo.

REM 在新窗口中启动
start "FLASH-DOCK - Port %PORT%" cmd /k wsl -d %WSL_DISTRO% bash -c "%WSL_CMD%"

echo [OK] 启动命令已执行
echo 请查看新打开的窗口
echo.
pause

