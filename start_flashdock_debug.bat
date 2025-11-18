@echo off
chcp 65001 >nul 2>&1
REM FlashDock WSL 启动脚本（调试版）
REM 显示详细输出，便于诊断问题

echo ============================================
echo FlashDock WSL 启动脚本（调试版）
echo ============================================
echo.

REM 配置
set WSL_DISTRO=Ubuntu-24.04
set ENV_NAME=flash_dock
set PROJECT_ROOT=/mnt/e/Qinchaojun/AIDD-TRAIN
set FLASHDOCK_DIR=%PROJECT_ROOT%/FLASH_DOCK-main
set PORT=8501

echo 配置信息:
echo   WSL发行版: %WSL_DISTRO%
echo   环境名称: %ENV_NAME%
echo   项目路径: %PROJECT_ROOT%
echo   端口: %PORT%
echo.

REM 测试环境
echo [1/4] 测试conda环境...
wsl -d %WSL_DISTRO% bash -c "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && python -c 'import streamlit; print(\"OK\")' 2>&1"
if errorlevel 1 (
    echo [错误] Conda环境测试失败
    pause
    exit /b 1
)
echo [OK] 环境测试通过
echo.

REM 测试文件
echo [2/4] 测试文件...
wsl -d %WSL_DISTRO% bash -c "test -f '%FLASHDOCK_DIR%/FlashDock.py' && echo 'FILE_EXISTS' || echo 'FILE_NOT_FOUND'"
if errorlevel 1 (
    echo [错误] FlashDock.py文件不存在
    pause
    exit /b 1
)
echo [OK] 文件存在
echo.

REM 清理旧进程
echo [3/4] 清理旧进程...
wsl -d %WSL_DISTRO% bash -c "pkill -f 'streamlit.*FlashDock' 2>/dev/null || true"
timeout /t 2 /nobreak >nul 2>&1
echo [OK] 清理完成
echo.

REM 启动服务
echo [4/4] 启动FlashDock...
echo.
echo 注意: 服务将在当前窗口运行，可以看到所有输出
echo 按 Ctrl+C 停止服务
echo.
echo ============================================
echo.

REM 构建启动命令
set "WSL_CMD=source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && export PYTHONPATH=%PROJECT_ROOT% && export BROWSER=none && cd %FLASHDOCK_DIR% && echo 'Starting Streamlit...' && streamlit run FlashDock.py --server.port %PORT% --server.address 0.0.0.0 --server.headless true"

REM 在当前窗口运行（可以看到输出）
wsl -d %WSL_DISTRO% bash -c "%WSL_CMD%"

echo.
echo ============================================
echo FlashDock 已停止
echo ============================================
pause

