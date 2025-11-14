# COMPASS 环境清单

本文档列出了 COMPASS 项目运行所需的所有环境依赖和版本要求。

## 环境名称

### Conda 环境名称
- **环境名称**: `AIDDTRAIN`
- **默认环境路径**: `C:\ProgramData\Anaconda3\envs\AIDDTRAIN` (Windows)
- **D盘环境路径**: `D:\conda_envs\AIDDTRAIN` (推荐，节省C盘空间)
- **激活命令**: 
  ```bash
  conda activate AIDDTRAIN
  ```

### 在D盘创建环境

#### 方法1: 使用自动配置脚本（推荐）
运行项目根目录下的 `create_compass_env_d.bat` 脚本，会自动：
1. 创建 `D:\conda_envs` 目录
2. 配置conda使用D盘作为环境目录
3. 创建 `AIDDTRAIN` 环境
4. 安装基础依赖

#### 方法2: 手动配置
```bash
# 1. 创建D盘环境目录
mkdir D:\conda_envs

# 2. 添加D盘到conda环境目录列表
conda config --add envs_dirs D:\conda_envs

# 3. 创建AIDDTRAIN环境（会自动在D盘创建）
conda create -n AIDDTRAIN python=3.12 -y

# 4. 激活环境
conda activate AIDDTRAIN

# 5. 升级pip
python -m pip install --upgrade pip
```

#### 验证环境位置
```bash
# 查看所有conda环境及其路径
conda info --envs

# 查看conda环境目录配置
conda config --show envs_dirs
```

### 虚拟环境
如果使用 Python venv，环境通常命名为 `venv` 或自定义名称。

## 系统要求

### 操作系统
- Windows 10/11
- Linux (Ubuntu 20.04+)
- macOS 10.15+

### Python 版本
- **推荐版本**: Python 3.12
- **支持版本**: Python 3.8, 3.9, 3.10, 3.11
- **当前环境**: Python 3.10.19

## 核心依赖

### 深度学习框架
- **PyTorch**: >= 2.8.0
  - 当前安装: 2.9.0+cpu
  - 安装命令（CUDA 12.8）:
    ```bash
    pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu128
    ```
  - 安装命令（CPU版本）:
    ```bash
    pip install torch torchvision
    ```

- **PyTorch Geometric**: >= 2.7.0
  - 当前安装: 2.7.0
  - 安装命令:
    ```bash
    pip install torch_geometric
    ```

- **PyTorch Geometric 扩展库**:
  - pyg_lib
  - torch_scatter
  - torch_sparse
  - torch_cluster
  - torch_spline_conv
  - 安装命令（CUDA 12.8）:
    ```bash
    pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.8.0+cu128.html
    ```

### 科学计算库
- **NumPy**: >= 1.26.0
  - 当前安装: 1.26.1

- **SciPy**: >= 1.10.0
  - 安装命令:
    ```bash
    pip install scipy
    ```

### 化学信息学库
- **RDKit**: 通过 rdkit-pypi 安装
  - 安装命令:
    ```bash
    pip install rdkit-pypi
    ```

- **BioPython**: >= 1.81
  - 安装命令:
    ```bash
    pip install biopython
    ```

### Web 框架（API 服务）
- **FastAPI**: >= 0.104.0
  - 当前安装: 0.121.0
  - 最低要求: 0.104.0

- **Uvicorn**: >= 0.24.0
  - 当前安装: 0.38.0
  - 安装命令:
    ```bash
    pip install uvicorn[standard]
    ```

- **Pydantic**: >= 2.0.0
  - 当前安装: 2.12.4
  - 最低要求: 2.0.0

- **python-multipart**: >= 0.0.6
  - 安装命令:
    ```bash
    pip install python-multipart
    ```

### 工具库
- **tqdm**: >= 4.66.0
  - 当前安装: 4.67.1
  - 安装命令:
    ```bash
    pip install tqdm
    ```

## 开发依赖

### 代码质量工具
- **black**: >= 23.0.0
  - 代码格式化工具

- **flake8**: >= 6.0.0
  - 代码风格检查

- **pylint**: >= 3.0.0
  - 代码质量检查

- **mypy**: >= 1.5.0
  - 类型检查工具

- **bandit**: >= 1.7.0
  - 安全漏洞扫描

### 测试工具
- **pytest**: >= 7.4.0
  - 测试框架

- **pytest-cov**: >= 4.1.0
  - 测试覆盖率

- **pytest-html**: >= 4.0.0
  - HTML 测试报告

### 类型存根
- **types-requests**: >= 2.31.0
- **types-pyyaml**: >= 6.0.15

## 完整安装命令

### 基础环境安装（推荐顺序）

```bash
# 1. 创建虚拟环境
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 2. 升级 pip
python -m pip install --upgrade pip

# 3. 安装 PyTorch（根据您的 CUDA 版本选择）
# CUDA 12.8
pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu128

# 或 CPU 版本
pip install torch torchvision

# 4. 安装 PyTorch Geometric 及其扩展
pip install torch_geometric
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.8.0+cu128.html

# 5. 安装其他核心依赖
pip install rdkit-pypi biopython tqdm scipy

# 6. 安装 Web 框架依赖
pip install fastapi>=0.104.0 uvicorn[standard]>=0.24.0 pydantic>=2.0.0 python-multipart>=0.0.6

# 7. 安装项目依赖
pip install -r requirements.txt

# 8. 安装开发依赖（可选）
pip install -r requirements-dev.txt
```

### 使用 requirements.txt 安装

```bash
# 安装生产依赖
pip install -r requirements.txt

# 安装开发依赖
pip install -r requirements-dev.txt
```

## 版本兼容性说明

### PyTorch 与 CUDA 版本对应关系
- PyTorch 2.8+ 支持 CUDA 11.8, 12.1, 12.4, 12.8
- 请根据您的 GPU 驱动版本选择合适的 CUDA 版本

### Python 版本兼容性
- Python 3.8-3.12 均支持
- 推荐使用 Python 3.10 或 3.12
- CI/CD 使用 Python 3.12

### ViSNet 模型约束
- **重要**: ViSNet 模型要求隐藏通道数必须能被注意力头数整除
- 配置模型时请确保满足此约束条件

## 验证安装

运行以下命令验证关键依赖是否正确安装：

```bash
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import torch_geometric; print(f'PyG: {torch_geometric.__version__}')"
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"
python -c "import fastapi; print(f'FastAPI: {fastapi.__version__}')"
python -c "from rdkit import Chem; print('RDKit: OK')"
```

## 当前环境状态

- **Python**: 3.10.19 (Anaconda)
- **PyTorch**: 2.9.0+cpu
- **torch-geometric**: 2.7.0
- **NumPy**: 1.26.1
- **FastAPI**: 0.121.0
- **Pydantic**: 2.12.4
- **Uvicorn**: 0.38.0
- **tqdm**: 4.67.1

## 注意事项

1. **CUDA 版本**: 如果使用 GPU，请确保 PyTorch 的 CUDA 版本与系统 CUDA 驱动兼容
2. **内存要求**: 训练大型模型需要足够的系统内存和显存
3. **数据准备**: 使用前需要下载 PDBbind 数据集并配置路径
4. **硬件优化**: 建议运行 `python -m compass.optimizer` 进行硬件配置优化

## 更新日期

最后更新: 2025年1月

