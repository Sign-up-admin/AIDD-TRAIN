@echo off
chcp 65001 >nul
echo ========================================
echo 停止 COMPASS 项目所有服务
echo ========================================
echo.

echo [步骤 1/3] 停止服务注册中心 (端口 8500)...
for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8500" ^| findstr "LISTENING"') do (
    echo 停止进程 %%a
    taskkill /F /PID %%a >nul 2>&1
    if !errorlevel! EQU 0 (
        echo   成功停止进程 %%a
    ) else (
        echo   进程 %%a 已停止或不存在
    )
)

echo [步骤 2/3] 停止 COMPASS 服务 (端口 8080)...
for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8080" ^| findstr "LISTENING"') do (
    echo 停止进程 %%a
    taskkill /F /PID %%a >nul 2>&1
    if !errorlevel! EQU 0 (
        echo   成功停止进程 %%a
    ) else (
        echo   进程 %%a 已停止或不存在
    )
)

echo [步骤 3/3] 停止 FLASH-DOCK 前端 (端口 8501)...
for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8501" ^| findstr "LISTENING"') do (
    echo 停止进程 %%a
    taskkill /F /PID %%a >nul 2>&1
    if !errorlevel! EQU 0 (
        echo   成功停止进程 %%a
    ) else (
        echo   进程 %%a 已停止或不存在
    )
)

echo.
echo 等待所有进程完全停止...
timeout /t 2 /nobreak >nul

echo.
echo 检查剩余进程...
netstat -ano | findstr ":8500 :8080 :8501" | findstr "LISTENING" >nul 2>&1
if errorlevel 1 (
    echo ========================================
    echo 所有服务已成功停止！
    echo ========================================
) else (
    echo 警告：仍有服务在运行，请检查：
    netstat -ano | findstr ":8500 :8080 :8501" | findstr "LISTENING"
)

echo.
pause


