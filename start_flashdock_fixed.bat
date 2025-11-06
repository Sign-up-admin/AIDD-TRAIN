@echo off
echo Starting FLASH-DOCK...
cd /d E:\Qinchaojun\AIDD-TRAIN\FLASH_DOCK-main
C:\ProgramData\Anaconda3\envs\flash_dock\python.exe -m streamlit run FlashDock.py --server.port 8501
pause

