@echo off
chcp 65001 >nul
cd /d %~dp0FLASH_DOCK-main
echo 启动 FLASH-DOCK 前端...
echo 端口: 8501
echo.
C:\ProgramData\Anaconda3\envs\flash_dock\python.exe -m streamlit run FlashDock.py --server.port 8501
pause

