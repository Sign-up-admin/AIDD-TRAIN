@echo off
REM Start Service Registry
echo Starting Service Registry...
conda activate AIDDTRAIN
python services/registry/server.py --host 0.0.0.0 --port 8500
pause

