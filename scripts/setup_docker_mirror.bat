@echo off
REM Docker 镜像加速器配置脚本（Windows）
REM 此脚本会帮助配置 Docker Desktop 的镜像加速器

echo ========================================
echo Docker 镜像加速器配置助手
echo ========================================
echo.
echo 此脚本将帮助您配置 Docker Desktop 使用国内镜像加速器
echo 以解决 Docker Hub 连接问题
echo.
pause

echo.
echo 请按照以下步骤操作：
echo.
echo 1. 打开 Docker Desktop
echo 2. 点击右上角的设置图标（齿轮）
echo 3. 选择 "Docker Engine"
echo 4. 在 JSON 配置中添加以下内容：
echo.
echo    "registry-mirrors": [
echo      "https://docker.mirrors.ustc.edu.cn",
echo      "https://hub-mirror.c.163.com",
echo      "https://mirror.baidubce.com"
echo    ]
echo.
echo 5. 点击 "Apply & Restart"
echo.
echo 配置完成后，Docker 将使用国内镜像源，速度会更快
echo.
pause



