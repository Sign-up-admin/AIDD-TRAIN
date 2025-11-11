@echo off
chcp 65001 >nul 2>&1
cd /d "%~dp0"

REM Get project root (parent directory)
cd /d %~dp0
set "FLASHDOCK_DIR=%CD%"
cd /d ..
set "PROJECT_ROOT=%CD%"
cd /d "%FLASHDOCK_DIR%"

echo ============================================
echo FlashDock Startup Script
echo ============================================
echo.
echo Project Root: %PROJECT_ROOT%
echo FlashDock Dir: %FLASHDOCK_DIR%
echo.

REM Set PYTHONPATH
set "PYTHONPATH=%PROJECT_ROOT%"

REM Find conda environments
set PYTHON_FLASHDOCK=python

if exist "%USERPROFILE%\anaconda3\envs\flash_dock\python.exe" (
    set PYTHON_FLASHDOCK=%USERPROFILE%\anaconda3\envs\flash_dock\python.exe
    echo Found conda environment at %USERPROFILE%\anaconda3\envs\flash_dock
) else (
    if exist "C:\ProgramData\Anaconda3\envs\flash_dock\python.exe" (
        set PYTHON_FLASHDOCK=C:\ProgramData\Anaconda3\envs\flash_dock\python.exe
        echo Found conda environment at C:\ProgramData\Anaconda3\envs\flash_dock
    ) else (
        echo Using default python...
    )
)

echo.
echo Starting FLASH-DOCK on port 8501...
echo.

cd /d "%FLASHDOCK_DIR%"
"%PYTHON_FLASHDOCK%" -m streamlit run FlashDock.py --server.port 8501

echo.
echo Browser will open automatically, or visit:
echo http://localhost:8501
echo.
echo Press Ctrl+C to stop server
echo ============================================
echo.
pause
