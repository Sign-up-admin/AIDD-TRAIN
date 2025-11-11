@echo off
REM Install code quality check tools in conda environment

set CONDA_ENV=AIDDTRAIN
set CONDA_BASE=C:\ProgramData\Anaconda3
set PYTHON_EXE=%CONDA_BASE%\envs\%CONDA_ENV%\python.exe

for %%I in ("%~dp0..") do set "PROJECT_ROOT=%%~fI"
set "PYTHONPATH=%PROJECT_ROOT%;%PYTHONPATH%"

echo [INFO] Installing code quality check tools in conda environment: %CONDA_ENV%
echo.

"%PYTHON_EXE%" -m pip install --quiet flake8 bandit black pylint mypy pytest pytest-cov pytest-html types-requests types-pyyaml

if errorlevel 1 (
    echo [ERROR] Failed to install tools
    exit /b 1
) else (
    echo [OK] Tools installed successfully
    exit /b 0
)



