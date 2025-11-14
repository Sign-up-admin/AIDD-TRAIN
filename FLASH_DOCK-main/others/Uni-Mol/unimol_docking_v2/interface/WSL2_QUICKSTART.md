# WSL2 快速开始指南

## 一键设置（推荐）

### 1. 进入 WSL2

```bash
wsl -d Ubuntu-24.04
```

### 2. 运行环境设置脚本

```bash
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface
bash setup_wsl2_env.sh
```

脚本会自动：
- 检查/安装 Miniconda
- 创建 `flash_dock` conda 环境
- 安装 PyTorch 和所有依赖
- 验证安装

### 3. 激活环境并运行

```bash
conda activate flash_dock
python3 demo.py --mode single \
  --input-protein /mnt/e/data/protein.pdb \
  --input-ligand /mnt/e/data/ligand.sdf \
  --input-docking-grid /mnt/e/data/grid.json \
  --output-ligand-name result \
  --output-ligand-dir /mnt/e/output
```

## 手动设置

如果自动脚本失败，可以手动设置：

### 1. 安装 Miniconda（如果未安装）

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### 2. 创建环境

```bash
conda create -n flash_dock python=3.10 -y
conda activate flash_dock
```

### 3. 安装依赖

```bash
# PyTorch (CPU)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# 项目依赖
pip install lmdb numpy pandas rdkit-pypi tqdm scikit-learn unicore
```

## 从 Windows 调用 WSL2

如果你在 Windows 上，可以使用包装脚本：

```powershell
cd E:\Qinchaojun\AIDD-TRAIN\FLASH_DOCK-main\others\Uni-Mol\unimol_docking_v2\interface

python run_in_wsl2.py \
  --mode single \
  --input-protein E:\data\protein.pdb \
  --input-ligand E:\data\ligand.sdf \
  --input-docking-grid E:\data\grid.json \
  --output-ligand-name result \
  --output-ligand-dir E:\output
```

## 路径说明

程序支持两种路径格式：

1. **WSL2 路径**（推荐）:
   ```
   /mnt/e/data/protein.pdb
   ```

2. **Windows 路径**（自动转换）:
   ```
   E:\data\protein.pdb
   ```

## 验证安装

```bash
conda activate flash_dock
python3 -c "import torch; print('PyTorch:', torch.__version__)"
python3 -c "import lmdb; print('LMDB: OK')"
python3 -c "from predictor import UnimolPredictor; print('UnimolPredictor: OK')"
```

## 常见问题

### Q: 如何检查 WSL2 是否已安装？

```powershell
wsl --list --verbose
```

### Q: 如何进入 WSL2？

```bash
wsl -d Ubuntu-24.04
```

### Q: 路径转换失败怎么办？

确保路径格式正确：
- Windows: `E:\path\to\file` 或 `E:/path/to/file`
- WSL2: `/mnt/e/path/to/file`

### Q: CUDA 不可用？

WSL2 中使用 GPU 需要：
1. Windows 上安装 NVIDIA 驱动
2. WSL2 中安装 CUDA toolkit

参考: https://docs.nvidia.com/cuda/wsl-user-guide/index.html

## 下一步

设置完成后，参考 `WSL2_MIGRATION_GUIDE.md` 了解详细使用方法。



