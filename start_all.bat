@echo off
chcp 65001 >nul
echo ========================================
echo 启动 COMPASS 完整项目
echo ========================================
echo.
echo 此脚本将启动以下服务：
echo   1. 服务注册中心 (端口 8500)
echo   2. COMPASS 服务 (端口 8080)
echo   3. FLASH-DOCK 前端 (端口 8501)
echo.

cd /d %~dp0
set PYTHONPATH=%CD%;%PYTHONPATH%

echo [1/3] 启动服务注册中心 (端口 8500)...
start "Service Registry - Port 8500" cmd /k "cd /d %~dp0 && set PYTHONPATH=%CD%;%PYTHONPATH% && C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe services/registry/server.py --host 0.0.0.0 --port 8500"
echo 等待服务注册中心启动...
timeout /t 5 /nobreak >nul

echo [2/3] 启动 COMPASS 服务 (端口 8080)...
start "COMPASS Service - Port 8080" cmd /k "cd /d %~dp0 && set PYTHONPATH=%CD%;%PYTHONPATH% && C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
echo 等待 COMPASS 服务启动...
timeout /t 5 /nobreak >nul

echo [3/3] 启动 FLASH-DOCK 前端 (端口 8501)...
start "FLASH-DOCK - Port 8501" cmd /k "cd /d %~dp0FLASH_DOCK-main && C:\ProgramData\Anaconda3\envs\flash_dock\python.exe -m streamlit run FlashDock.py --server.port 8501"
echo 等待 FLASH-DOCK 启动...
timeout /t 3 /nobreak >nul

echo.
echo ========================================
echo 所有服务启动完成！
echo ========================================
echo.
echo 服务地址：
echo   - 服务注册中心: http://localhost:8500
echo   - COMPASS 服务: http://localhost:8080
echo   - FLASH-DOCK 前端: http://localhost:8501
echo.
echo API 文档：
echo   - COMPASS API: http://localhost:8080/docs
echo   - 健康检查: http://localhost:8080/health
echo.
echo 提示：
echo   - 每个服务都在独立的窗口中运行
echo   - 关闭窗口即可停止对应服务
echo   - 浏览器应自动打开 FLASH-DOCK 界面
echo.
pause
