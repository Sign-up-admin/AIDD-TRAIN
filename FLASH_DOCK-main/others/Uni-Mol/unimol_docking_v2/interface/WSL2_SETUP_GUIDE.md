# WSL2 使用指南

## 为什么使用 WSL2？

Windows 上的访问违规错误（0xC0000005）很可能是由于：
- Windows 的 C 扩展库兼容性问题
- Windows 的线程模型与 Linux 不同
- LMDB 等库在 Windows 上的实现问题

WSL2 运行在真正的 Linux 内核上，可以完全避免这些问题。

## 安装 WSL2

### 1. 检查 WSL 是否已安装

在 PowerShell（管理员权限）中运行：

```powershell
wsl --list --verbose
```

如果已安装，会显示已安装的发行版。

### 2. 安装 WSL2（如果未安装）

```powershell
# 启用 WSL 功能
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart

# 启用虚拟机平台
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart

# 重启计算机（必需）

# 下载并安装 WSL2 内核更新
# 访问：https://aka.ms/wsl2kernel

# 设置 WSL2 为默认版本
wsl --set-default-version 2
```

### 3. 安装 Ubuntu（推荐）

```powershell
# 从 Microsoft Store 安装 Ubuntu，或使用命令行：
wsl --install -d Ubuntu-22.04
```

## 在 WSL2 中设置环境

### 1. 进入 WSL2

```bash
wsl
# 或
ubuntu
```

### 2. 安装 Miniconda/Anaconda

```bash
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

# 安装依赖
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
# 如果有 GPU，使用：
# pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

pip install lmdb numpy pandas rdkit-pypi tqdm scikit-learn
pip install unicore
```

### 4. 访问 Windows 文件系统

WSL2 中可以通过 `/mnt/` 访问 Windows 驱动器：

```bash
# Windows 的 E: 盘在 WSL2 中对应 /mnt/e/
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2
```

## 在 WSL2 中运行程序

### 方法 1: 直接在 WSL2 中运行

```bash
# 进入 WSL2
wsl

# 激活环境
conda activate flash_dock

# 进入项目目录
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2

# 运行程序（使用 Linux 路径）
python interface/demo.py \
    --mode single \
    --model-dir /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2 \
    --input-protein <蛋白质文件路径> \
    --input-ligand <配体文件路径> \
    --input-docking-grid <网格文件路径> \
    --output-ligand-name <输出名称> \
    --output-ligand-dir <输出目录路径> \
    --conf-size 10 \
    --batch-size 4  # 在 Linux 上可以使用更大的 batch_size
```

### 方法 2: 从 Windows 调用 WSL2

在 Windows PowerShell 中：

```powershell
wsl -e bash -c "cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2 && conda activate flash_dock && python interface/demo.py --mode single ..."
```

## 路径转换

### Windows 路径 → WSL2 路径

- `E:\Qinchaojun\AIDD-TRAIN` → `/mnt/e/Qinchaojun/AIDD-TRAIN`
- `C:\Users\...` → `/mnt/c/Users/...`

### 注意事项

1. **路径分隔符**: Windows 使用 `\`，Linux 使用 `/`
2. **大小写敏感**: Linux 路径是大小写敏感的
3. **性能**: 访问 `/mnt/` 下的文件可能比原生 Linux 文件系统慢一些

## 性能优化

### 1. 将项目复制到 WSL2 文件系统

为了更好的性能，可以将项目复制到 WSL2 的文件系统：

```bash
# 在 WSL2 中
cp -r /mnt/e/Qinchaojun/AIDD-TRAIN ~/AIDD-TRAIN
cd ~/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2
```

### 2. 使用符号链接（如果数据文件很大）

```bash
# 只复制代码，数据文件使用符号链接
ln -s /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/Result ~/AIDD-TRAIN/FLASH_DOCK-main/Result
```

## GPU 支持（如果适用）

如果需要在 WSL2 中使用 GPU：

1. **安装 NVIDIA 驱动**（在 Windows 上）
   - 下载并安装最新的 NVIDIA 驱动
   - 确保支持 WSL2

2. **在 WSL2 中安装 CUDA**

```bash
# 安装 CUDA Toolkit（如果使用 GPU）
# 参考 NVIDIA 官方文档
```

## 故障排除

### 问题 1: WSL2 无法启动

```powershell
# 检查 WSL 状态
wsl --status

# 如果问题持续，尝试重置
wsl --shutdown
wsl --unregister Ubuntu
wsl --install -d Ubuntu-22.04
```

### 问题 2: 文件权限问题

```bash
# 修复文件权限
sudo chown -R $USER:$USER /mnt/e/Qinchaojun/AIDD-TRAIN
```

### 问题 3: 性能问题

- 将项目复制到 WSL2 文件系统（`~/` 目录）
- 避免频繁访问 `/mnt/` 下的文件

## 优势总结

✅ **完全避免 Windows 兼容性问题**
✅ **可以使用更大的 batch_size**（Linux 上更稳定）
✅ **更好的多线程支持**
✅ **LMDB 等库在 Linux 上更稳定**
✅ **可以使用 GPU**（如果安装了 NVIDIA 驱动）

## 总结

WSL2 是解决 Windows 上访问违规错误的最佳方案。虽然需要一些设置，但一旦配置完成，程序应该能够稳定运行，不再出现崩溃问题。

