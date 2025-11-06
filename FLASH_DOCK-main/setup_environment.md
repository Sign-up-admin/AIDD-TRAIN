# Flash_Dock 环境配置指南

## 已创建的Conda环境

环境名称：`flash_dock`
- Python版本：3.8
- 位置：`C:\ProgramData\Anaconda3\envs\flash_dock`

## 激活环境

```bash
# 使用完整路径激活（如果conda不在PATH中）
C:\ProgramData\Anaconda3\Scripts\activate.bat flash_dock

# 或者如果conda已添加到PATH
conda activate flash_dock
```

## 安装依赖

### 方法1：使用requirements.txt（推荐）

```bash
# 激活环境后
pip install -r requirements.txt
```

### 方法2：分步安装

```bash
# 1. 安装基础依赖
pip install streamlit pandas "numpy<2.0.0" pyyaml tqdm scikit-learn biopandas==0.4.1

# 2. 安装RDKit（注意：需要numpy<2.0.0）
pip install rdkit-pypi==2022.9.3 -i https://pypi.tuna.tsinghua.edu.cn/simple/ --trusted-host pypi.tuna.tsinghua.edu.cn

# 3. 安装可视化工具
pip install py3dmol stmol streamlit-molstar streamlit-ketcher

# 4. 安装PyTorch（根据系统选择）
# CPU版本：
pip install torch

# GPU版本（CUDA 11.3）：
pip install torch -f https://download.pytorch.org/whl/torch_stable.html

# 5. 安装Uni-Mol相关
pip install huggingface_hub
pip install unimol_tools

# 或者从源码安装Uni-Core（如果pip安装失败）
pip install git+https://github.com/dptech-corp/Uni-Core.git@stable
```

## 重要依赖说明

1. **numpy版本限制**：必须使用 `numpy<2.0.0`，因为RDKit需要这个版本
2. **Java环境**：需要安装Java 17-23用于P2Rank（系统级别安装）
3. **模型文件**：需要下载 `unimol_docking_v2_240517.pt` (464MB) 并放到 `./others/Uni-Mol/unimol_docking_v2/`

## 验证安装

```bash
# 激活环境
conda activate flash_dock

# 检查Python版本
python --version  # 应该显示 Python 3.8.x

# 检查关键包
python -c "import streamlit; print('Streamlit:', streamlit.__version__)"
python -c "import rdkit; print('RDKit:', rdkit.__version__)"
python -c "import torch; print('PyTorch:', torch.__version__)"
```

## 启动应用

```bash
# 激活环境后
streamlit run FlashDock.py
```

## 常见问题

### 1. RDKit安装失败
```bash
# 尝试使用conda安装
conda install -c conda-forge rdkit
```

### 2. numpy版本冲突
```bash
# 确保numpy版本正确
pip install "numpy<2.0.0"
```

### 3. 找不到conda命令
```bash
# 使用完整路径
C:\ProgramData\Anaconda3\Scripts\conda.exe activate flash_dock
```


