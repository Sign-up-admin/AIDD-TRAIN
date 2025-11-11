@echo off
chcp 65001 >nul 2>&1
cd /d %~dp0

echo ========================================
echo Starting COMPASS Microservices
echo ========================================
echo.

REM Find conda and python executables
set PYTHON_AIDDTRAIN=python
set PYTHON_FLASHDOCK=python

REM Try to find conda environments
if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
    set PYTHON_AIDDTRAIN=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe
    set PYTHON_FLASHDOCK=%USERPROFILE%\anaconda3\envs\flash_dock\python.exe
    echo Found conda environments at %USERPROFILE%\anaconda3
) else (
    if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
        set PYTHON_AIDDTRAIN=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe
        set PYTHON_FLASHDOCK=C:\ProgramData\Anaconda3\envs\flash_dock\python.exe
        echo Found conda environments at C:\ProgramData\Anaconda3
    ) else (
        echo Conda environments not found, using default python
    )
)

REM Get project root directory (absolute path)
cd /d %~dp0
set "PROJECT_ROOT=%CD%"

REM Set PYTHONPATH
set "PYTHONPATH=%PROJECT_ROOT%"

echo Project Root: %PROJECT_ROOT%
echo.

echo Step 1: Starting Service Registry...
if exist "%PYTHON_AIDDTRAIN%" (
    start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_AIDDTRAIN% services\registry\server.py --host 0.0.0.0 --port 8500 && pause"
) else (
    start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python services\registry\server.py --host 0.0.0.0 --port 8500 && pause"
)
echo Waiting for Service Registry to start...
timeout /t 5 /nobreak >nul 2>&1

echo Step 2: Starting COMPASS Service...
if exist "%PYTHON_AIDDTRAIN%" (
    start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_AIDDTRAIN% compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500 && pause"
) else (
    start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500 && pause"
)
echo Waiting for COMPASS service to register...
timeout /t 5 /nobreak >nul 2>&1

echo Step 3: Starting FLASH-DOCK...
if exist "%PYTHON_FLASHDOCK%" (
    start "FLASH-DOCK - Port 8501" cmd /k "cd /d %PROJECT_ROOT%\FLASH_DOCK-main && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_FLASHDOCK% -m streamlit run FlashDock.py --server.port 8501 && pause"
) else (
    start "FLASH-DOCK - Port 8501" cmd /k "cd /d %PROJECT_ROOT%\FLASH_DOCK-main && set PYTHONPATH=%PROJECT_ROOT% && python -m streamlit run FlashDock.py --server.port 8501 && pause"
)

echo.
echo ========================================
echo All services are starting in separate windows
echo ========================================
echo Service Registry: http://localhost:8500
echo COMPASS Service: http://localhost:8080
echo FLASH-DOCK: http://localhost:8501
echo.
echo Please check the opened windows for service status
echo If services fail to start, check the error messages in those windows
echo.
pause
