@echo off
chcp 65001 >nul 2>&1
cd /d "%~dp0"

echo ============================================
echo FlashDock Startup Script
echo ============================================
echo.

REM Check if conda is available
where conda >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo Using conda environment flash_dock...
    conda activate flash_dock
    python -m streamlit run FlashDock.py --server.port 8501
) else (
    echo Conda not found, using default python...
    python -m streamlit run FlashDock.py --server.port 8501
)

echo.
echo Browser will open automatically, or visit:
echo http://localhost:8501
echo.
echo Press Ctrl+C to stop server
echo ============================================
echo.
pause
