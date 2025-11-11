@echo off
setlocal enabledelayedexpansion
chcp 65001 >nul 2>&1
cd /d %~dp0

echo ========================================
echo Stop all COMPASS services
echo ========================================
echo.

echo [Step 1/3] Stopping Service Registry (port 8500)...
call :stop_port 8500

echo [Step 2/3] Stopping COMPASS Service (port 8080)...
call :stop_port 8080

echo [Step 3/3] Stopping FLASH-DOCK Frontend (port 8501)...
call :stop_port 8501

echo.
echo Waiting for all processes to stop...
timeout /t 2 /nobreak >nul 2>&1

echo.
echo Checking remaining processes...
netstat -ano 2>nul | findstr ":8500 :8080 :8501" | findstr "LISTENING" >nul 2>&1
if errorlevel 1 (
    echo ========================================
    echo All services stopped successfully!
    echo ========================================
) else (
    echo Warning: Some services may still be running:
    netstat -ano 2>nul | findstr ":8500 :8080 :8501" | findstr "LISTENING"
)

echo.
pause
goto :eof

:stop_port
set PORT=%1
set FOUND=0

for /f "tokens=5" %%p in ('netstat -ano 2^>nul ^| findstr ":%PORT%" ^| findstr "LISTENING"') do (
    if not "%%p"=="" (
        echo Stopping process %%p
        taskkill /F /PID %%p >nul 2>&1
        if !errorlevel! EQU 0 (
            echo   Successfully stopped process %%p
        ) else (
            echo   Process %%p already stopped or not found
        )
        set FOUND=1
    )
)

if !FOUND! EQU 0 (
    echo   No service running on port %PORT%
)

exit /b 0
