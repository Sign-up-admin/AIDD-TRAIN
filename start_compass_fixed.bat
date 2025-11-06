@echo off
echo Starting COMPASS Service...
cd /d E:\Qinchaojun\AIDD-TRAIN
C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
pause

