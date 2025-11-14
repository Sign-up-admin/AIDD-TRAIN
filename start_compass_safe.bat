@echo off
chcp 65001 >nul 2>&1
setlocal enabledelayedexpansion

REM Get project root (absolute path)
cd /d "%~dp0"
set "PROJECT_ROOT=%CD%"

echo ========================================
echo Starting COMPASS Service
echo ========================================
echo.
echo Project Root: %PROJECT_ROOT%
echo.

REM Check if project root exists
if not exist "%PROJECT_ROOT%" (
    echo [ERROR] Project root does not exist: %PROJECT_ROOT%
    pause
    exit /b 1
)

REM Check if service_main.py exists
if not exist "%PROJECT_ROOT%\compass\service_main.py" (
    echo [ERROR] Service entry point not found: %PROJECT_ROOT%\compass\service_main.py
    echo Please ensure you are running this script from the project root directory.
    pause
    exit /b 1
)

REM Check port 8080
echo Checking port 8080...
netstat -an | findstr ":8080" | findstr "LISTENING" >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo [WARNING] Port 8080 is already in use!
    echo Please stop the existing service or use a different port.
    echo.
    choice /C YN /M "Continue anyway"
    if errorlevel 2 exit /b 1
)

REM Set PYTHONPATH (use absolute path)
set "PYTHONPATH=%PROJECT_ROOT%"

REM Find conda environment AIDDTRAIN (check multiple locations)
set "PYTHON_AIDDTRAIN="
set "ENV_FOUND=0"

REM Check D:\conda_envs
if exist "D:\conda_envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=D:\conda_envs\AIDDTRAIN\python.exe"
    set "ENV_FOUND=1"
    echo [OK] Found AIDDTRAIN at D:\conda_envs
    goto :found_env
)

REM Check C:\ProgramData\Anaconda3
if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe"
    set "ENV_FOUND=1"
    echo [OK] Found AIDDTRAIN at C:\ProgramData\Anaconda3
    goto :found_env
)

REM Check %USERPROFILE%\anaconda3
if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe"
    set "ENV_FOUND=1"
    echo [OK] Found AIDDTRAIN at %USERPROFILE%\anaconda3
    goto :found_env
)

REM Check %USERPROFILE%\miniconda3
if exist "%USERPROFILE%\miniconda3\envs\AIDDTRAIN\python.exe" (
    set "PYTHON_AIDDTRAIN=%USERPROFILE%\miniconda3\envs\AIDDTRAIN\python.exe"
    set "ENV_FOUND=1"
    echo [OK] Found AIDDTRAIN at %USERPROFILE%\miniconda3
    goto :found_env
)

:found_env
REM Verify Python executable exists and is valid
if %ENV_FOUND% EQU 1 (
    if not exist "!PYTHON_AIDDTRAIN!" (
        echo [ERROR] Python executable not found: !PYTHON_AIDDTRAIN!
        set "ENV_FOUND=0"
    ) else (
        REM Test Python version
        "!PYTHON_AIDDTRAIN!" --version >nul 2>&1
        if !ERRORLEVEL! NEQ 0 (
            echo [WARNING] Python executable may be invalid: !PYTHON_AIDDTRAIN!
            set "ENV_FOUND=0"
        )
    )
)

REM Start service based on environment detection
if %ENV_FOUND% EQU 1 (
    echo.
    echo Using conda environment AIDDTRAIN: !PYTHON_AIDDTRAIN!
    echo.
    
    REM Verify critical dependencies before starting
    echo Checking critical dependencies...
    "!PYTHON_AIDDTRAIN!" -c "import fastapi, uvicorn, torch" >nul 2>&1
    if !ERRORLEVEL! NEQ 0 (
        echo [WARNING] Some critical dependencies may be missing.
        echo Please run: pip install -r requirements.txt
        echo.
        choice /C YN /M "Continue anyway"
        if errorlevel 2 exit /b 1
    ) else (
        echo [OK] Critical dependencies found
    )
    
    echo.
    echo Starting COMPASS service...
    echo.
    start "COMPASS Service - Port 8080" cmd /k "cd /d \"%PROJECT_ROOT%\" && set PYTHONPATH=%PROJECT_ROOT% && \"!PYTHON_AIDDTRAIN!\" compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
    
) else (
    REM Try to use conda activate
    where conda >nul 2>&1
    if %ERRORLEVEL% EQU 0 (
        echo.
        echo [INFO] AIDDTRAIN environment not found in standard locations.
        echo Attempting to use conda activate...
        echo.
        
        REM Test if conda activate works
        call conda activate AIDDTRAIN >nul 2>&1
        if %ERRORLEVEL% EQU 0 (
            python --version >nul 2>&1
            if %ERRORLEVEL% EQU 0 (
                echo [OK] Using conda activate AIDDTRAIN
                echo.
                start "COMPASS Service - Port 8080" cmd /k "cd /d \"%PROJECT_ROOT%\" && conda activate AIDDTRAIN && set PYTHONPATH=%PROJECT_ROOT% && python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
            ) else (
                goto :use_default
            )
        ) else (
            goto :use_default
        )
    ) else (
        :use_default
        echo.
        echo [WARNING] AIDDTRAIN conda environment not found!
        echo Using default Python interpreter.
        echo.
        echo [IMPORTANT] Make sure:
        echo   1. AIDDTRAIN environment is created: conda create -n AIDDTRAIN python=3.10
        echo   2. Dependencies are installed: pip install -r requirements.txt
        echo.
        choice /C YN /M "Continue with default Python"
        if errorlevel 2 exit /b 1
        
        REM Verify default Python
        python --version >nul 2>&1
        if %ERRORLEVEL% NEQ 0 (
            echo [ERROR] Python not found in PATH!
            echo Please install Python or activate the AIDDTRAIN environment.
            pause
            exit /b 1
        )
        
        echo.
        echo Starting COMPASS service with default Python...
        echo.
        start "COMPASS Service - Port 8080" cmd /k "cd /d \"%PROJECT_ROOT%\" && set PYTHONPATH=%PROJECT_ROOT% && python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
    )
)

echo.
echo ========================================
echo Service startup initiated
echo ========================================
echo Service URL: http://localhost:8080
echo Registry URL: http://localhost:8500
echo.
echo [TIP] Check the new window for service status and logs.
echo [TIP] If service fails to start, check the error messages in that window.
echo [TIP] Run diagnose_compass.py for detailed diagnostics.
echo.
pause
