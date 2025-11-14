# FlashDock WSL 使用指南

本指南说明如何在 WSL 的 conda 虚拟环境中运行 FlashDock。

## 前提条件

1. **WSL2 已安装并运行**
   - 推荐使用 Ubuntu 24.04 或更高版本
   - 检查 WSL 状态：`wsl --status`

2. **项目文件可访问**
   - 项目路径：`/mnt/e/Qinchaojun/AIDD-TRAIN`
   - 确保 FlashDock 目录存在：`/mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main`

## 快速开始

### 步骤 1: 设置环境

在 WSL 中运行环境设置脚本：

```bash
# 进入项目目录
cd /mnt/e/Qinchaojun/AIDD-TRAIN

# 运行设置脚本
bash setup_flashdock_wsl.sh
```

这个脚本会：
- 检查并安装 Miniconda（如果未安装）
- 创建名为 `flash_dock_wsl` 的 conda 环境（Python 3.8）
- 安装所有必需的依赖包
- 验证安装是否成功

**注意**：首次运行可能需要较长时间（10-30分钟），取决于网络速度。

### 步骤 2: 启动 FlashDock

环境设置完成后，使用启动脚本运行 FlashDock：

```bash
# 在项目根目录
bash start_flashdock_wsl.sh
```

或者手动启动：

```bash
# 激活环境
conda activate flash_dock_wsl

# 进入 FlashDock 目录
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main

# 设置 PYTHONPATH
export PYTHONPATH=/mnt/e/Qinchaojun/AIDD-TRAIN:$PYTHONPATH

# 启动 Streamlit
streamlit run FlashDock.py --server.port 8501 --server.address 0.0.0.0
```

### 步骤 3: 访问 FlashDock

启动成功后，在浏览器中访问：
- **本地访问**: http://localhost:8501
- **从 Windows 访问**: http://localhost:8501（WSL2 自动端口转发）

## 环境信息

- **Conda 环境名称**: `flash_dock_wsl`
- **Python 版本**: 3.8
- **主要依赖**:
  - Streamlit >= 1.28.0
  - RDKit (rdkit-pypi==2022.9.3)
  - NumPy < 2.0.0
  - PyTorch (CPU 版本)
  - streamlit-molstar, py3dmol, stmol
  - 其他依赖见 `FLASH_DOCK-main/requirements.txt`

## 常见问题

### 问题 1: Conda 命令未找到

**解决方案**:
```bash
# 如果使用 Miniconda
export PATH="$HOME/miniconda3/bin:$PATH"
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# 如果使用 Anaconda
export PATH="$HOME/anaconda3/bin:$PATH"
source "$HOME/anaconda3/etc/profile.d/conda.sh"
```

### 问题 2: 端口已被占用

**解决方案**:
```bash
# 检查端口占用
netstat -tuln | grep 8501

# 或者使用其他端口启动
streamlit run FlashDock.py --server.port 8502 --server.address 0.0.0.0
```

### 问题 3: RDKit 安装失败

**解决方案**:
```bash
# 激活环境
conda activate flash_dock_wsl

# 使用 conda 安装 RDKit
conda install -c conda-forge rdkit -y
```

### 问题 4: 依赖安装缓慢

**解决方案**:
脚本已配置使用清华大学镜像源加速。如果仍然缓慢，可以手动设置：

```bash
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
pip config set global.trusted-host pypi.tuna.tsinghua.edu.cn
```

### 问题 5: 无法从 Windows 访问

**解决方案**:
- WSL2 会自动转发端口，确保使用 `--server.address 0.0.0.0`
- 检查 Windows 防火墙设置
- 尝试使用 WSL 的 IP 地址访问（在 WSL 中运行 `hostname -I` 获取 IP）

## 手动管理环境

### 激活环境
```bash
conda activate flash_dock_wsl
```

### 查看已安装的包
```bash
conda list
# 或
pip list
```

### 更新依赖
```bash
conda activate flash_dock_wsl
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main
pip install -r requirements.txt --upgrade
```

### 删除环境（重新开始）
```bash
conda env remove -n flash_dock_wsl
# 然后重新运行 setup_flashdock_wsl.sh
```

## 性能优化建议

1. **使用 WSL2 原生文件系统**（可选）
   - 将项目复制到 `~/AIDD-TRAIN` 可能获得更好的性能
   - 但需要确保所有路径引用正确

2. **使用 GPU 版本 PyTorch**（如果有 NVIDIA GPU）
   ```bash
   conda activate flash_dock_wsl
   pip uninstall torch torchvision torchaudio
   pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
   ```

3. **配置 Streamlit 缓存**
   - 在 `~/.streamlit/config.toml` 中配置缓存设置

## 与 Windows 环境的区别

| 特性 | Windows 环境 | WSL 环境 |
|------|-------------|----------|
| 环境名称 | `flash_dock` | `flash_dock_wsl` |
| Python 路径 | `C:\ProgramData\Anaconda3\envs\flash_dock\python.exe` | `~/miniconda3/envs/flash_dock_wsl/bin/python` |
| 项目路径 | `E:\Qinchaojun\AIDD-TRAIN` | `/mnt/e/Qinchaojun/AIDD-TRAIN` |
| 启动方式 | `.bat` 脚本 | `.sh` 脚本 |

## 验证安装

运行以下命令验证环境是否正确设置：

```bash
conda activate flash_dock_wsl
python -c "import streamlit; print('Streamlit:', streamlit.__version__)"
python -c "import rdkit; print('RDKit:', rdkit.__version__)"
python -c "import torch; print('PyTorch:', torch.__version__)"
python -c "import numpy; print('NumPy:', numpy.__version__)"
```

所有命令应该都能成功执行并显示版本号。

## 停止服务

在运行 FlashDock 的终端中按 `Ctrl+C` 停止服务。

## 日志和调试

如果遇到问题，可以查看 Streamlit 的详细日志：

```bash
streamlit run FlashDock.py --server.port 8501 --server.address 0.0.0.0 --logger.level=debug
```

## 相关文件

- `setup_flashdock_wsl.sh` - 环境设置脚本
- `start_flashdock_wsl.sh` - 启动脚本
- `FLASH_DOCK-main/requirements.txt` - 依赖列表
- `FLASH_DOCK-main/FlashDock.py` - 主程序文件

---

*最后更新: 2025-01-XX*

