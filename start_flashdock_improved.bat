@echo off
chcp 65001 >nul 2>&1
REM FlashDock WSL 启动脚本（改进版）
REM 修复了启动问题，添加了更好的错误处理

echo ============================================
echo FlashDock WSL 启动脚本（改进版）
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

REM 检查端口是否被占用
netstat -ano | findstr ":%PORT%" | findstr "LISTENING" >nul 2>&1
if not errorlevel 1 (
    echo [警告] 端口 %PORT% 已被占用
    echo 正在尝试停止占用进程...
    for /f "tokens=5" %%p in ('netstat -ano ^| findstr ":%PORT%" ^| findstr "LISTENING"') do (
        taskkill /F /PID %%p >nul 2>&1
    )
    timeout /t 2 /nobreak >nul 2>&1
)

REM 停止WSL中可能存在的streamlit进程
echo 清理WSL中的旧进程...
wsl -d %WSL_DISTRO% bash -c "pkill -f 'streamlit.*FlashDock' 2>/dev/null || true" >nul 2>&1
timeout /t 1 /nobreak >nul 2>&1

REM 构建启动命令
REM 注意：禁用浏览器自动打开，避免WSL中的错误
set "WSL_CMD=source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate %ENV_NAME% && export PYTHONPATH=%PROJECT_ROOT% && export BROWSER=none && cd %FLASHDOCK_DIR% && streamlit run FlashDock.py --server.port %PORT% --server.address 0.0.0.0 --server.headless true 2>&1"

echo 正在启动 FlashDock...
echo 访问地址: http://localhost:%PORT%
echo.
echo 提示: 服务将在新窗口中启动
echo 如果看到 "gio: Operation not supported" 警告，可以忽略（这是WSL的正常现象）
echo.

REM 在新窗口中启动（不等待输出，避免阻塞）
start "FLASH-DOCK - Port %PORT%" cmd /k "wsl -d %WSL_DISTRO% bash -c \"%WSL_CMD%\""

echo [OK] 启动命令已执行
echo.
echo 等待服务启动（约10秒）...
timeout /t 10 /nobreak >nul 2>&1

REM 检查服务是否启动
echo.
echo 检查服务状态...
powershell -Command "try { $response = Invoke-WebRequest -Uri 'http://localhost:%PORT%' -TimeoutSec 3 -UseBasicParsing; Write-Host '[OK] FLASH-DOCK 服务已启动' -ForegroundColor Green } catch { Write-Host '[WARNING] 服务可能还在启动中，请检查新窗口' -ForegroundColor Yellow }"

echo.
echo ============================================
echo 如果服务未启动，请检查新窗口中的错误信息
echo 访问地址: http://localhost:%PORT%
echo ============================================
pause

