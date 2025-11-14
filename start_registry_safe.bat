@echo off
chcp 65001 >nul 2>&1
cd /d "%~dp0"

REM Get project root (relative to script location)
set "PROJECT_ROOT=%~dp0"
set "PROJECT_ROOT=%PROJECT_ROOT:~0,-1%"

echo ========================================
echo Starting Service Registry
echo ========================================
echo.
echo Project Root: %PROJECT_ROOT%
echo.

REM Set PYTHONPATH
set "PYTHONPATH=%PROJECT_ROOT%"

REM Find conda environment AIDDTRAIN
set "PYTHON_AIDDTRAIN="
if exist "D:\conda_envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=D:\conda_envs\AIDDTRAIN\python.exe"
    echo Found AIDDTRAIN at D:\conda_envs
) else if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe"
    echo Found AIDDTRAIN at C:\ProgramData\Anaconda3
) else if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe"
    echo Found AIDDTRAIN at %USERPROFILE%\anaconda3
)

REM Check if conda is available and use it, otherwise use direct path
if defined PYTHON_AIDDTRAIN (
    echo Using conda environment AIDDTRAIN: %PYTHON_AIDDTRAIN%
    start "Service Registry - Port 8500" cmd /k "cd /d \"%PROJECT_ROOT%\" && set PYTHONPATH=%PROJECT_ROOT% && \"%PYTHON_AIDDTRAIN%\" services\registry\server.py --host 0.0.0.0 --port 8500"
) else (
    where conda >nul 2>&1
    if %ERRORLEVEL% EQU 0 (
        echo Using conda activate AIDDTRAIN...
        start "Service Registry - Port 8500" cmd /k "cd /d \"%PROJECT_ROOT%\" && conda activate AIDDTRAIN && set PYTHONPATH=%PROJECT_ROOT% && python services\registry\server.py --host 0.0.0.0 --port 8500"
    ) else (
        echo Using default python...
        start "Service Registry - Port 8500" cmd /k "cd /d \"%PROJECT_ROOT%\" && set PYTHONPATH=%PROJECT_ROOT% && python services\registry\server.py --host 0.0.0.0 --port 8500"
    )
)

echo.
echo Service started in new window
echo Port: 8500
echo.
pause
