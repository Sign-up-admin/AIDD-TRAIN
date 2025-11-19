@echo off
chcp 65001 >nul 2>&1
echo ========================================
echo Java安装快速指南
echo ========================================
echo.
echo 当前系统未检测到Java环境
echo.
echo 请按照以下步骤安装Java：
echo.
echo 【方法一：推荐】使用Adoptium（免费开源）
echo 1. 访问：https://adoptium.net/zh-CN/temurin/releases/
echo 2. 选择：Version=17或21, OS=Windows, Arch=x64, Package=JDK
echo 3. 下载并安装，勾选"Set JAVA_HOME"和"Add to PATH"
echo 4. 安装后重启命令行和Streamlit应用
echo.
echo 【方法二】使用Oracle JDK
echo 1. 访问：https://www.oracle.com/java/technologies/downloads/
echo 2. 选择Java 17或21的Windows x64 Installer
echo 3. 安装后需要手动配置环境变量
echo.
echo ========================================
echo.
echo 是否要打开Adoptium下载页面？
echo 按 Y 打开浏览器，按 N 退出
choice /C YN /M "请选择"
if errorlevel 2 exit /b 0
if errorlevel 1 (
    start https://adoptium.net/zh-CN/temurin/releases/
    echo.
    echo 浏览器已打开，请按照页面提示下载Java
    echo.
    echo 安装完成后，请运行以下命令验证：
    echo   java -version
    echo.
    echo 然后重新运行诊断脚本：
    echo   python check_java_environment.py
)
pause

