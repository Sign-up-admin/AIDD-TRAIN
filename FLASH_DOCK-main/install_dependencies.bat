@echo off
REM Flash_Dock 依赖安装脚本
REM 使用flash_dock环境安装依赖

echo ========================================
echo Flash_Dock 依赖安装脚本
echo ========================================
echo.

REM 设置conda路径
set CONDA_PATH=C:\ProgramData\Anaconda3
set ENV_NAME=flash_dock
set PYTHON_PATH=%CONDA_PATH%\envs\%ENV_NAME%\python.exe

echo 检查环境: %ENV_NAME%
if not exist "%PYTHON_PATH%" (
    echo 错误: 环境 %ENV_NAME% 不存在或Python路径不正确
    echo 请确保已创建flash_dock环境
    pause
    exit /b 1
)

echo.
echo 当前Python版本:
"%PYTHON_PATH%" --version

echo.
echo ========================================
echo 开始安装依赖包...
echo ========================================
echo.

REM 使用清华镜像源安装
echo [1/6] 安装基础依赖...
"%PYTHON_PATH%" -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn streamlit pandas "numpy<2.0.0" pyyaml tqdm scikit-learn biopandas==0.4.1

if errorlevel 1 (
    echo 警告: 基础依赖安装失败，可能是网络问题
    echo 请检查网络连接或使用代理
)

echo.
echo [2/6] 安装RDKit...
"%PYTHON_PATH%" -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn rdkit-pypi==2022.9.3

if errorlevel 1 (
    echo 警告: RDKit安装失败，尝试使用conda安装...
    conda install -c conda-forge rdkit -y
)

echo.
echo [3/6] 安装可视化工具...
"%PYTHON_PATH%" -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn py3dmol stmol streamlit-molstar streamlit-ketcher

echo.
echo [4/6] 安装PyTorch (CPU版本)...
"%PYTHON_PATH%" -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn torch

echo.
echo [5/6] 安装Uni-Mol相关...
"%PYTHON_PATH%" -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn huggingface_hub

REM 尝试安装unimol_tools（可能失败，需要手动处理）
echo.
echo [6/6] 尝试安装unimol_tools...
"%PYTHON_PATH%" -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn unimol_tools

if errorlevel 1 (
    echo.
    echo 警告: unimol_tools安装失败，可能需要从源码安装
    echo 请参考: https://github.com/deepmodeling/Uni-Mol
)

echo.
echo ========================================
echo 安装完成！
echo ========================================
echo.
echo 验证安装:
"%PYTHON_PATH%" -c "import streamlit; print('Streamlit:', streamlit.__version__)"
"%PYTHON_PATH%" -c "import rdkit; print('RDKit:', rdkit.__version__)"
"%PYTHON_PATH%" -c "import torch; print('PyTorch:', torch.__version__)"

echo.
echo 下一步:
echo 1. 下载模型文件: unimol_docking_v2_240517.pt (464MB)
echo 2. 放置到: ./others/Uni-Mol/unimol_docking_v2/
echo 3. 运行: streamlit run FlashDock.py
echo.
pause

