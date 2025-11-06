@echo off
chcp 65001 >nul
echo ============================================
echo FlashDock 项目启动脚本（修复版）
echo ============================================
echo.

REM 获取脚本所在目录
cd /d "%~dp0"

REM 设置环境变量
set PYTHON_EXE=C:\ProgramData\Anaconda3\envs\flash_dock\python.exe
set STREAMLIT_SCRIPT=C:\ProgramData\Anaconda3\envs\flash_dock\Scripts\streamlit.exe

REM 检查环境是否存在
if not exist "%PYTHON_EXE%" (
    echo 错误：找不到 flash_dock 环境
    echo Python 路径: %PYTHON_EXE%
    echo.
    echo 请确保 flash_dock 环境已创建
    pause
    exit /b 1
)

REM 检查 streamlit 是否安装
if not exist "%STREAMLIT_SCRIPT%" (
    echo 警告：找不到 streamlit，尝试使用 python -m streamlit...
    echo 正在启动 FlashDock...
    echo.
    "%PYTHON_EXE%" -m streamlit run FlashDock.py --server.port 8501
) else (
    echo 使用 streamlit 启动 FlashDock...
    echo.
    "%STREAMLIT_SCRIPT%" run FlashDock.py --server.port 8501
)

echo.
echo 浏览器将自动打开，如果没有自动打开，请访问：
echo http://localhost:8501
echo.
echo 按 Ctrl+C 停止服务器
echo ============================================
echo.

pause

