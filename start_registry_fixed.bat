@echo off
echo Starting Service Registry...
cd /d E:\Qinchaojun\AIDD-TRAIN
C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe services/registry/server.py --host 0.0.0.0 --port 8500
pause

