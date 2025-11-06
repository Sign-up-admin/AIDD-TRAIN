@echo off
chcp 65001 >nul
echo ========================================
echo 启动 COMPASS 微服务系统（简化版）
echo ========================================
echo.
echo 此脚本直接使用 Python 启动服务
echo 请确保已正确配置 Python 环境
echo.

echo [1/3] 启动服务注册中心 (端口 8500)...
start "Service Registry - Port 8500" cmd /k "python services/registry/server.py --host 0.0.0.0 --port 8500"
echo 等待服务注册中心启动...
timeout /t 5 /nobreak >nul

echo [2/3] 启动 COMPASS 服务 (端口 8080)...
start "COMPASS Service - Port 8080" cmd /k "python compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500"
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
echo.
pause

