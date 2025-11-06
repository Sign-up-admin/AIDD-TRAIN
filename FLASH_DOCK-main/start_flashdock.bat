@echo off
echo ============================================
echo FlashDock 项目启动脚本
echo ============================================
echo.

REM 获取脚本所在目录
cd /d "%~dp0"

REM 初始化conda（如果需要）
call "C:\ProgramData\Anaconda3\Scripts\activate.bat" 2>nul
if errorlevel 1 (
    call "%USERPROFILE%\anaconda3\Scripts\activate.bat" 2>nul
    if errorlevel 1 (
        call "%USERPROFILE%\miniconda3\Scripts\activate.bat" 2>nul
        if errorlevel 1 (
            echo 正在查找conda安装路径...
            where conda >nul 2>&1
            if errorlevel 1 (
                echo 错误：未找到conda，请确保已安装Anaconda或Miniconda
                pause
                exit /b 1
            )
        )
    )
)

echo 正在激活flash_dock虚拟环境...
call conda activate flash_dock
if errorlevel 1 (
    echo.
    echo 错误：无法激活flash_dock环境
    echo 请确保环境已创建，可以使用以下命令创建：
    echo   conda create -n flash_dock python=3.10
    echo   conda activate flash_dock
    echo   pip install -r requirements.txt
    echo.
    pause
    exit /b 1
)

echo.
echo 环境激活成功！
echo 正在启动FlashDock项目...
echo.
echo 浏览器将自动打开，如果没有自动打开，请访问：
echo http://localhost:8501
echo.
echo 按 Ctrl+C 停止服务器
echo ============================================
echo.

streamlit run FlashDock.py

pause

