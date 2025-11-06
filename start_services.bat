@echo off
chcp 65001 >nul
cd /d %~dp0
echo ========================================
echo 启动 COMPASS 微服务系统
echo ========================================
echo.
echo 工作目录: %CD%
echo.

REM 设置PYTHONPATH
set PYTHONPATH=%CD%;%PYTHONPATH%

REM 设置conda环境路径
set CONDA_ENV_PATH=
if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
    set CONDA_ENV_PATH=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe
)

if "%CONDA_ENV_PATH%"=="" (
    if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
        set CONDA_ENV_PATH=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe
    )
)

REM 决定使用哪个Python
if not "%CONDA_ENV_PATH%"=="" (
    echo 使用 conda 环境: %CONDA_ENV_PATH%
    set PYTHON_CMD=%CONDA_ENV_PATH%
) else (
    echo 使用系统 Python
    set PYTHON_CMD=python
)

echo.
echo [1/3] 启动服务注册中心 (端口 8500)...
start "Service Registry - Port 8500" cmd /k "cd /d %~dp0 && set PYTHONPATH=%CD%;%PYTHONPATH% && %PYTHON_CMD% services/registry/server.py --host 0.0.0.0 --port 8500"
echo 等待服务注册中心启动...
timeout /t 5 /nobreak >nul

echo [2/3] 启动 COMPASS 服务 (端口 8080)...
start "COMPASS Service - Port 8080" cmd /k "cd /d %~dp0 && set PYTHONPATH=%CD%;%PYTHONPATH% && %PYTHON_CMD% compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
echo 等待 COMPASS 服务启动...
timeout /t 5 /nobreak >nul

echo.
echo ========================================
echo 服务启动完成！
echo ========================================
echo.
echo 服务地址：
echo   - 服务注册中心: http://localhost:8500
echo   - COMPASS 服务: http://localhost:8080
echo   - API 文档: http://localhost:8080/docs
echo   - 健康检查: http://localhost:8080/health
echo   - 性能指标: http://localhost:8080/metrics
echo.
echo 提示：
echo   - 每个服务都在独立的窗口中运行
echo   - 关闭窗口即可停止对应服务
echo   - 查看日志请查看各服务窗口输出
echo.
pause
