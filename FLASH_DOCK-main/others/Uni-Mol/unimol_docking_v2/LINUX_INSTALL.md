# Linux 系统 Uni-Core 安装指南

## 问题

运行 Uni-Mol Docking V2 时出现错误：
```
ModuleNotFoundError: No module named 'unicore'
```

## 解决方案

在 Linux 系统（如 `root@Print150`）上执行以下命令：

### 步骤 1: 激活 conda 环境
```bash
conda activate flash_dock
```

### 步骤 2: 安装 Uni-Core

**推荐方法（一键安装）：**
```bash
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2
bash install_for_linux.sh
```

**手动安装（如果脚本失败）：**
```bash
pip install git+https://github.com/dptech-corp/Uni-Core.git@stable
```

**备用方法（如果网络问题）：**
```bash
git clone https://github.com/dptech-corp/Uni-Core.git
cd Uni-Core
pip install -e .
cd ..
```

### 步骤 3: 验证安装
```bash
python -c "import unicore; print('✓ 安装成功！')"
```

## 常见问题

### 1. 网络连接问题
如果无法访问 GitHub，可以：
- 配置代理
- 使用备用方法（克隆后安装）

### 2. 权限问题
如果遇到权限错误，可以使用 `--user`：
```bash
pip install --user git+https://github.com/dptech-corp/Uni-Core.git@stable
```

### 3. Git 未安装
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install git

# CentOS/RHEL
sudo yum install git
```

## 安装完成后

验证安装成功后，就可以正常使用 Uni-Mol Docking V2 了：

```bash
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface
python demo.py [参数]
```

## 相关文件

- `install_for_linux.sh` - Linux 专用安装脚本
- `INSTALL_UNICORE.md` - 详细安装指南
- `install_unicore.py` - 跨平台安装脚本
