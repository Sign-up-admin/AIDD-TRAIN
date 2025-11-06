@echo off
chcp 65001 >nul
cd /d %~dp0

REM 获取绝对路径
set PROJECT_ROOT=%CD%

REM 检查conda环境
set PYTHON_CMD=python
if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
    set PYTHON_CMD=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe
) else (
    if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
        set PYTHON_CMD=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe
    )
)

echo ========================================
echo 启动服务注册中心
echo ========================================
echo.
echo 项目目录: %PROJECT_ROOT%
echo Python: %PYTHON_CMD%
echo.

REM 设置PYTHONPATH - 只使用绝对路径
set PYTHONPATH=%PROJECT_ROOT%

REM 切换到项目目录
cd /d %PROJECT_ROOT%

REM 启动服务
echo 正在启动服务注册中心...
echo 端口: 8500
echo.

start "Service Registry - Port 8500" cmd /k "cd /d %PROJECT_ROOT% && set PYTHONPATH=%PROJECT_ROOT% && %PYTHON_CMD% services\registry\server.py --host 0.0.0.0 --port 8500"

echo.
echo 服务已在新窗口中启动
echo 请查看新窗口中的服务状态
echo.
pause

