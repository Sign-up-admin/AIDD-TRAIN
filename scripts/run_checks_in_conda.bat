@echo off
REM Run code quality checks in conda environment
REM This script activates the conda environment and runs the checks

setlocal enabledelayedexpansion

for %%I in ("%~dp0..") do set "PROJECT_ROOT=%%~fI"
set "PYTHONPATH=%PROJECT_ROOT%;%PYTHONPATH%"

REM Try to find conda environment
set CONDA_ENV=AIDDTRAIN
set CONDA_BASE=C:\ProgramData\Anaconda3

REM Check if conda environment exists
if exist "%CONDA_BASE%\envs\%CONDA_ENV%\python.exe" (
    echo [INFO] Found conda environment: %CONDA_ENV%
    echo [INFO] Running code checks in conda environment...
    echo.
    
    REM Check if --ci flag is already in arguments
    set HAS_CI=0
    for %%a in (%*) do (
        if "%%a"=="--ci" set HAS_CI=1
    )
    
    REM Add --ci flag if not present (to auto-install tools and non-interactive mode)
    if %HAS_CI%==0 (
        echo [INFO] Adding --ci flag for non-interactive mode (auto-install tools)
        "%CONDA_BASE%\envs\%CONDA_ENV%\python.exe" "%~dp0run_all_checks.py" --ci %*
    ) else (
        "%CONDA_BASE%\envs\%CONDA_ENV%\python.exe" "%~dp0run_all_checks.py" %*
    )
    
    if errorlevel 1 (
        echo.
        echo [ERROR] Code checks failed with exit code %errorlevel%
        exit /b %errorlevel%
    ) else (
        echo.
        echo [OK] Code checks completed successfully
        exit /b 0
    )
) else (
    echo [ERROR] Conda environment %CONDA_ENV% not found at %CONDA_BASE%\envs\%CONDA_ENV%
    echo [INFO] Please check conda environment name and path
    exit /b 1
)

