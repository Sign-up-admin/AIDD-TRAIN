@echo off
REM Start COMPASS Service
echo Starting COMPASS Service...
conda activate AIDDTRAIN
python compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
pause

