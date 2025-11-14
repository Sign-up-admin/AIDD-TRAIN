# WSL2 迁移指南

## 概述

本指南帮助你将 Uni-Mol 对接预测程序从 Windows 迁移到 WSL2，以解决访问违规错误（0xC0000005）。

## 为什么使用 WSL2？

Windows 上的访问违规错误很可能是由于：
- Windows 的 C 扩展库兼容性问题
- LMDB 在 Windows 上的多线程支持有限
- PyTorch DataLoader 在 Windows 上的 multiprocessing 问题
- unicore 库在 Windows 上的兼容性问题

WSL2 运行在真正的 Linux 内核上，可以完全避免这些问题。

## 快速开始

### 方法 1: 直接在 WSL2 中运行（推荐）

1. **进入 WSL2**
   ```bash
   wsl -d Ubuntu-24.04
   ```

2. **激活 conda 环境**
   ```bash
   conda activate flash_dock
   ```

3. **进入项目目录**
   ```bash
   cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface
   ```

4. **运行程序**
   ```bash
   # 使用 Linux 路径
   python3 demo.py \
     --mode single \
     --input-protein /mnt/e/data/protein.pdb \
     --input-ligand /mnt/e/data/ligand.sdf \
     --input-docking-grid /mnt/e/data/grid.json \
     --output-ligand-name result \
     --output-ligand-dir /mnt/e/output
   
   # 或使用 Windows 路径（会自动转换）
   python3 demo.py \
     --mode single \
     --input-protein E:\data\protein.pdb \
     --input-ligand E:\data\ligand.sdf \
     --input-docking-grid E:\data\grid.json \
     --output-ligand-name result \
     --output-ligand-dir E:\output
   ```

### 方法 2: 从 Windows 调用 WSL2

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

## 环境设置

### 1. 检查 WSL2 是否已安装

```powershell
wsl --list --verbose
```

如果看到 `Ubuntu-24.04`，说明已安装。

### 2. 在 WSL2 中安装 Miniconda（如果未安装）

```bash
# 进入 WSL2
wsl -d Ubuntu-24.04

# 下载 Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 安装
bash Miniconda3-latest-Linux-x86_64.sh

# 重新加载 shell
source ~/.bashrc
```

### 3. 创建 conda 环境

```bash
# 创建环境
conda create -n flash_dock python=3.10 -y
conda activate flash_dock

# 安装 PyTorch（CPU 版本）
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# 如果有 GPU，使用 CUDA 版本
# pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# 安装项目依赖
pip install lmdb numpy pandas rdkit-pypi tqdm scikit-learn unicore
```

## 路径转换

程序会自动处理路径转换：

- **Windows 路径** → **WSL2 路径**
  - `E:\data\file.pdb` → `/mnt/e/data/file.pdb`
  - `C:\Users\file.pdb` → `/mnt/c/Users/file.pdb`

- **WSL2 路径** → **Windows 路径**
  - `/mnt/e/data/file.pdb` → `E:\data\file.pdb`

## 性能优化建议

1. **将项目复制到 WSL2 文件系统**
   ```bash
   # 在 WSL2 中
   cp -r /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main ~/AIDD-TRAIN/
   cd ~/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface
   ```

2. **使用更大的 batch_size**
   - Windows: 自动限制为 1
   - WSL2: 可以使用 4 或更大

3. **启用多线程数据加载**
   - Windows: `num_workers=0`
   - WSL2: `num_workers=8`（默认）

## 故障排除

### 问题 1: 无法访问 Windows 文件系统

```bash
# 检查文件系统挂载
ls /mnt/e

# 如果无法访问，检查 WSL2 状态
wsl --status
```

### 问题 2: 路径转换失败

确保路径格式正确：
- Windows: `E:\path\to\file` 或 `E:/path/to/file`
- WSL2: `/mnt/e/path/to/file`

### 问题 3: CUDA 不可用

WSL2 中使用 GPU 需要：
1. 在 Windows 上安装 NVIDIA 驱动
2. 在 WSL2 中安装 CUDA toolkit

参考: https://docs.nvidia.com/cuda/wsl-user-guide/index.html

### 问题 4: 权限问题

```bash
# 修复文件权限
sudo chown -R $USER:$USER /mnt/e/Qinchaojun/AIDD-TRAIN
```

## 优势总结

✅ **完全避免 Windows 兼容性问题**
✅ **可以使用更大的 batch_size**（Linux 上更稳定）
✅ **更好的多线程支持**
✅ **LMDB 等库在 Linux 上更稳定**
✅ **可以使用 GPU**（如果安装了 NVIDIA 驱动）

## 下一步

1. 按照本指南设置 WSL2 环境
2. 在 WSL2 中测试运行程序
3. 验证程序不再崩溃
4. 享受稳定的运行体验！

## 相关文件

- `run_in_wsl2.sh` - WSL2 运行脚本（在 WSL2 中使用）
- `run_in_wsl2.py` - Windows 调用 WSL2 的包装脚本
- `demo.py` - 主程序（已支持 WSL2 路径转换）



