@echo off
chcp 65001 >nul 2>&1
cd /d "%~dp0"

echo ========================================
echo Restarting COMPASS Microservices
echo ========================================
echo.

REM Get project root (absolute path)
cd /d %~dp0
set "PROJECT_ROOT=%CD%"

REM Find conda python executables (will be set later)
set PYTHON_AIDDTRAIN=python
set PYTHON_FLASHDOCK=python

echo [Step 1/5] Stopping existing services...
echo Finding and stopping processes on ports 8500, 8080, and 8501...

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8500" ^| findstr "LISTENING"') do (
    echo Stopping process %%a (port 8500)
    taskkill /F /PID %%a >nul 2>&1
)

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8080" ^| findstr "LISTENING"') do (
    echo Stopping process %%a (port 8080)
    taskkill /F /PID %%a >nul 2>&1
)

for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8501" ^| findstr "LISTENING"') do (
    echo Stopping process %%a (port 8501)
    taskkill /F /PID %%a >nul 2>&1
)

echo Waiting for services to stop completely...
timeout /t 3 /nobreak >nul 2>&1

echo.
echo [Step 2/5] Starting Service Registry (port 8500)...
echo Project Root: %PROJECT_ROOT%
echo.

REM Set PYTHONPATH
set "PYTHONPATH=%PROJECT_ROOT%"

REM Find conda python executables
set PYTHON_AIDDTRAIN=python
set PYTHON_FLASHDOCK=python

if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
    set PYTHON_AIDDTRAIN=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe
    set PYTHON_FLASHDOCK=%USERPROFILE%\anaconda3\envs\flash_dock\python.exe
) else (
    if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
        set PYTHON_AIDDTRAIN=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe
        set PYTHON_FLASHDOCK=C:\ProgramData\Anaconda3\envs\flash_dock\python.exe
    )
)

if exist "%PYTHON_AIDDTRAIN%" (
    start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_AIDDTRAIN% services\registry\server.py --host 0.0.0.0 --port 8500"
) else (
    start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python services\registry\server.py --host 0.0.0.0 --port 8500"
)

echo Waiting for Service Registry to start...
timeout /t 5 /nobreak >nul 2>&1

echo.
echo [Step 3/5] Starting COMPASS Service (port 8080)...
if exist "%PYTHON_AIDDTRAIN%" (
    start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_AIDDTRAIN% compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
) else (
    start "COMPASS Service - Port 8080" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
)

echo Waiting for COMPASS service to register...
timeout /t 5 /nobreak >nul 2>&1

echo.
echo [Step 4/5] Starting FLASH-DOCK (port 8501)...
if exist "%PYTHON_FLASHDOCK%" (
    start "FLASH-DOCK - Port 8501" cmd /k "cd /d %PROJECT_ROOT%\FLASH_DOCK-main && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_FLASHDOCK% -m streamlit run FlashDock.py --server.port 8501"
) else (
    start "FLASH-DOCK - Port 8501" cmd /k "cd /d %PROJECT_ROOT%\FLASH_DOCK-main && set PYTHONPATH=%PROJECT_ROOT% && python -m streamlit run FlashDock.py --server.port 8501"
)

echo.
echo [Step 5/5] All services restarted
echo ========================================
echo Service Registry: http://localhost:8500
echo COMPASS Service: http://localhost:8080
echo FLASH-DOCK: http://localhost:8501
echo.
pause
