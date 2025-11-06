@echo off
echo ========================================
echo Starting COMPASS Microservices
echo ========================================
echo.

echo Step 1: Starting Service Registry...
start "Service Registry" cmd /k "conda activate AIDDTRAIN && cd /d E:\Qinchaojun\AIDD-TRAIN && python services/registry/server.py --host 0.0.0.0 --port 8500"
timeout /t 3 /nobreak >nul

echo Step 2: Starting COMPASS Service...
start "COMPASS Service" cmd /k "conda activate AIDDTRAIN && cd /d E:\Qinchaojun\AIDD-TRAIN && python compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
timeout /t 3 /nobreak >nul

echo Step 3: Starting FLASH-DOCK...
start "FLASH-DOCK" cmd /k "conda activate flash_dock && cd /d E:\Qinchaojun\AIDD-TRAIN\FLASH_DOCK-main && streamlit run FlashDock.py --server.port 8501"

echo.
echo ========================================
echo All services are starting in separate windows
echo ========================================
echo Service Registry: http://localhost:8500
echo COMPASS Service: http://localhost:8080
echo FLASH-DOCK: http://localhost:8501
echo.
pause

