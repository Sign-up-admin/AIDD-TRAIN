@echo off
chcp 65001 >nul 2>&1
echo ========================================
echo 在D盘创建 COMPASS (AIDDTRAIN) 环境
echo ========================================
echo.

REM 检查conda是否可用
where conda >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [错误] 未找到conda命令，请先安装Anaconda或Miniconda
    pause
    exit /b 1
)

echo [1/4] 配置conda环境目录到D盘...
if not exist D:\conda_envs (
    mkdir D:\conda_envs
    echo 已创建目录: D:\conda_envs
)
conda config --add envs_dirs D:\conda_envs
echo 已添加D:\conda_envs到conda环境目录列表
echo.

echo [2/4] 创建AIDDTRAIN环境 (Python 3.12)...
conda create -n AIDDTRAIN python=3.12 -y
if %ERRORLEVEL% NEQ 0 (
    echo [错误] 环境创建失败
    pause
    exit /b 1
)
echo.

echo [3/4] 激活环境并升级pip...
call conda activate AIDDTRAIN
python -m pip install --upgrade pip
echo.

echo [4/4] 安装PyTorch (CPU版本，如需GPU请修改)...
pip install torch torchvision
echo.

echo ========================================
echo 环境创建完成！
echo ========================================
echo.
echo 环境名称: AIDDTRAIN
echo 环境路径: D:\conda_envs\AIDDTRAIN
echo.
echo 激活环境命令:
echo   conda activate AIDDTRAIN
echo.
echo 接下来请安装其他依赖:
echo   pip install -r requirements.txt
echo.
pause

