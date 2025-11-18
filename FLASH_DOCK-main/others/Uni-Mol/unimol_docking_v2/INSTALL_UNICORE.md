# Uni-Core 安装指南

## 问题描述

如果遇到以下错误：
```
ModuleNotFoundError: No module named 'unicore'
```

这是因为 `unicore` 不是 PyPI 上的独立包，而是 Uni-Core 项目的一部分，需要从 GitHub 安装。

## 快速安装

### 方法 1: 使用安装脚本（推荐）

在 Linux/WSL2 系统上：

```bash
# 激活 conda 环境
conda activate flash_dock

# 运行安装脚本
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2
python install_unicore.py
```

或者使用 bash 脚本：

```bash
bash install_unicore.sh
```

### 方法 2: 手动安装

#### 步骤 1: 从 GitHub 直接安装（最简单）

```bash
conda activate flash_dock
pip install git+https://github.com/dptech-corp/Uni-Core.git@stable
```

#### 步骤 2: 如果方法 1 失败（网络问题）

```bash
# 使用 git:// 协议
pip install "git+git://github.com/dptech-corp/Uni-Core.git@stable"
```

#### 步骤 3: 如果仍有问题，克隆后安装

```bash
# 克隆仓库
git clone https://github.com/dptech-corp/Uni-Core.git
cd Uni-Core
git checkout stable  # 切换到 stable 分支

# 安装
pip install -e .
cd ..
```

## 验证安装

安装完成后，验证是否成功：

```bash
python -c "import unicore; print('✓ unicore 模块导入成功')"
python -c "from unicore import checkpoint_utils, distributed_utils, options, utils; print('✓ 核心模块导入成功')"
```

## 常见问题

### 1. 网络连接问题

如果无法访问 GitHub，可以：

- 配置 Git 代理：
  ```bash
  git config --global http.proxy http://proxy.example.com:8080
  git config --global https.proxy https://proxy.example.com:8080
  ```

- 使用镜像源（如果有）

### 2. Git 未安装

**Linux:**
```bash
sudo apt-get update
sudo apt-get install git
```

**Windows:**
从 https://git-scm.com/download/win 下载安装

### 3. 权限问题

如果遇到权限问题，可以使用 `--user` 标志：

```bash
pip install --user git+https://github.com/dptech-corp/Uni-Core.git@stable
```

### 4. 安装后仍无法导入

- 确保在正确的 conda 环境中
- 重新启动 Python 解释器
- 检查 Python 路径：`python -c "import sys; print(sys.path)"`

## 相关资源

- Uni-Core GitHub: https://github.com/dptech-corp/Uni-Core
- Uni-Mol Docking V2 README: [README.md](README.md)

