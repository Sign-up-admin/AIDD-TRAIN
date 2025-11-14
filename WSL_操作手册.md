# WSL 操作手册

## 当前WSL状态

### WSL版本信息
- **WSL版本**: 2.5.10.0
- **内核版本**: 6.6.87.2-1
- **WSLg版本**: 1.0.66
- **Windows版本**: 10.0.26100.7171

### 已安装的发行版
1. **docker-desktop** (默认，运行中，WSL2)
2. **Ubuntu-24.04** (已停止，WSL2)

### Ubuntu-24.04 详细信息
- **发行版**: Ubuntu 24.04.3 LTS (noble)
- **内核**: Linux 6.6.87.2-microsoft-standard-WSL2
- **当前用户**: root
- **当前目录**: /mnt/e/Qinchaojun/AIDD-TRAIN
- **Python版本**: Python 3.12.3
- **Conda**: 未安装
- **项目目录**: 存在（/mnt/e/Qinchaojun/AIDD-TRAIN）

---

## 常用操作

### 检查WSL状态

使用项目提供的检查脚本：
```bash
python check_wsl_status.py
```

或使用WSL命令：
```bash
# 查看WSL版本
wsl --version

# 查看WSL状态
wsl --status

# 列出所有发行版
wsl --list --verbose
```

### 启动Ubuntu-24.04

```bash
# 启动Ubuntu-24.04
wsl -d Ubuntu-24.04

# 或设置为默认后直接启动
wsl --set-default Ubuntu-24.04
wsl
```

### 设置默认发行版

```bash
# 设置Ubuntu-24.04为默认发行版
wsl --set-default Ubuntu-24.04
```

### 停止WSL发行版

```bash
# 停止特定发行版
wsl --terminate Ubuntu-24.04

# 停止所有发行版
wsl --shutdown
```

### 在WSL中访问Windows文件

WSL中可以通过 `/mnt/` 访问Windows驱动器：
- `E:\Qinchaojun\AIDD-TRAIN` → `/mnt/e/Qinchaojun/AIDD-TRAIN`
- `C:\Users\...` → `/mnt/c/Users/...`

### 从Windows访问WSL文件

在Windows资源管理器中输入：
```
\\wsl$\Ubuntu-24.04\home\用户名
```

或在PowerShell中：
```powershell
cd \\wsl$\Ubuntu-24.04\home\用户名
```

---

## 在WSL中使用项目

### 进入项目目录

```bash
# 在WSL中
cd /mnt/e/Qinchaojun/AIDD-TRAIN
```

### 安装Conda（如需要）

```bash
# 下载Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 安装
bash Miniconda3-latest-Linux-x86_64.sh

# 重新加载shell
source ~/.bashrc

# 验证安装
conda --version
```

### 创建Conda环境

```bash
# 创建环境
conda create -n AIDDTRAIN python=3.10 -y

# 激活环境
conda activate AIDDTRAIN

# 安装依赖
pip install -r requirements.txt
```

---

## 故障排除

### 问题1: WSL无法启动

**解决方案**:
```bash
# 检查WSL服务状态
wsl --status

# 重启WSL
wsl --shutdown
wsl -d Ubuntu-24.04
```

### 问题2: 无法访问Windows文件

**解决方案**:
- 确保路径使用 `/mnt/` 前缀
- 检查文件权限：`ls -la /mnt/e/`

### 问题3: 网络连接问题

**解决方案**:
- WSL使用Windows的网络配置
- 如需代理，在WSL中设置环境变量：
```bash
export http_proxy=http://127.0.0.1:10808
export https_proxy=http://127.0.0.1:10808
```

### 问题4: 权限问题

**解决方案**:
- 在Windows文件系统上，WSL以当前用户权限运行
- 如需root权限，使用 `sudo` 命令

---

## 性能优化建议

1. **使用WSL2原生文件系统**（推荐）
   - 将项目复制到WSL2原生文件系统（如 `~/AIDD-TRAIN`）
   - 性能比访问 `/mnt/` 下的Windows文件系统更好

2. **避免在Windows和WSL之间频繁切换**
   - 在WSL中完成所有操作，减少文件系统切换

3. **使用WSL2而不是WSL1**
   - 当前已使用WSL2，性能更好

---

## 相关工具

- **检查WSL状态**: `python check_wsl_status.py`
- **WSL帮助**: `wsl --help`

---

## 注意事项

1. WSL2使用虚拟化技术，需要启用Windows的虚拟化功能
2. 在WSL中修改Windows文件系统上的文件可能影响性能
3. 建议将项目文件放在WSL2原生文件系统中以获得最佳性能
4. 当前默认发行版是 `docker-desktop`，如需使用Ubuntu-24.04，请设置默认发行版

---

*最后更新: 2025-11-14*

