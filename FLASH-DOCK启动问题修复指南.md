# FLASH-DOCK 启动问题修复指南

## 🔍 诊断结果

根据诊断工具的结果，发现以下问题：

### ✅ 正常的部分
- WSL已安装并启用
- WSL发行版 `Ubuntu-24.04` 可用
- Conda环境 `flash_dock` 存在
- 项目文件完整
- 端口8501可用
- 启动命令测试通过

### ⚠️ 发现的问题
1. **环境名称不一致**: 
   - `start_flashdock_wsl.bat` 使用 `flash_dock_wsl`
   - 实际环境名称是 `flash_dock`
   - **已修复**: 已将脚本中的环境名称改为 `flash_dock`

2. **依赖检查网络问题**: 
   - WSL网络连接可能有问题，但实际依赖已安装
   - Streamlit可以正常导入

---

## 🚀 启动方法

### 方法1: 使用修复后的启动脚本（推荐）

```bash
start_flashdock_fixed.bat
```

这个脚本会：
- ✅ 自动检查WSL环境
- ✅ 验证conda环境
- ✅ 检查项目文件
- ✅ 在新窗口中启动服务
- ✅ 自动检查服务状态

### 方法2: 使用原始启动脚本（已修复环境名称）

```bash
start_flashdock_wsl.bat
```

### 方法3: 使用Python启动脚本

```bash
python fix_and_start_services.py
```

这会启动所有服务，包括FLASH-DOCK。

### 方法4: 手动在WSL中启动

如果自动启动失败，可以手动在WSL中启动：

```bash
# 1. 进入WSL
wsl -d Ubuntu-24.04

# 2. 激活conda环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate flash_dock

# 3. 设置环境变量
export PYTHONPATH=/mnt/e/Qinchaojun/AIDD-TRAIN

# 4. 进入FlashDock目录
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main

# 5. 启动Streamlit
streamlit run FlashDock.py --server.port 8501 --server.address 0.0.0.0
```

---

## 🔧 常见问题解决

### 问题1: 环境名称不匹配

**症状**: 提示 "Conda 环境 'flash_dock_wsl' 不存在"

**解决方案**: 
- ✅ 已修复 `start_flashdock_wsl.bat` 中的环境名称
- 或使用 `start_flashdock_fixed.bat`

### 问题2: 端口被占用

**症状**: 提示 "端口 8501 已被占用"

**解决方案**:
```bash
# 在WSL中停止占用端口的进程
wsl -d Ubuntu-24.04 bash -c "lsof -ti :8501 | xargs kill -9 2>/dev/null || true"

# 或在Windows中停止
netstat -ano | findstr :8501
taskkill /F /PID <PID>
```

### 问题3: Streamlit启动失败

**症状**: Streamlit无法启动或报错

**解决方案**:
```bash
# 在WSL中重新安装streamlit
wsl -d Ubuntu-24.04 bash -c "source ~/miniconda3/etc/profile.d/conda.sh && conda activate flash_dock && pip install --upgrade streamlit"
```

### 问题4: 依赖缺失

**症状**: 导入模块失败

**解决方案**:
```bash
# 在WSL中安装依赖
wsl -d Ubuntu-24.04 bash -c "source ~/miniconda3/etc/profile.d/conda.sh && conda activate flash_dock && pip install -r /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/requirements.txt"
```

### 问题5: WSL网络问题

**症状**: WSL无法连接网络或执行命令超时

**解决方案**:
```bash
# 重启WSL
wsl --shutdown
wsl -d Ubuntu-24.04

# 检查网络
ping google.com
```

---

## 📋 启动检查清单

启动前请确认：

- [ ] WSL已安装并启用
- [ ] WSL发行版 `Ubuntu-24.04` 存在
- [ ] Conda环境 `flash_dock` 存在
- [ ] 项目目录 `/mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main` 存在
- [ ] `FlashDock.py` 文件存在
- [ ] 端口8501未被占用
- [ ] Streamlit已安装

运行诊断工具检查：
```bash
python diagnose_flashdock.py
```

---

## 🎯 验证服务运行

### 方法1: 浏览器访问
打开浏览器访问: http://localhost:8501

### 方法2: 命令行检查
```bash
# 使用Python检查
python -c "import requests; r = requests.get('http://localhost:8501', timeout=3); print('FLASH-DOCK运行中' if r.status_code < 500 else '响应异常')"
```

### 方法3: 检查进程
```bash
# 在WSL中检查
wsl -d Ubuntu-24.04 bash -c "ps aux | grep streamlit"
```

---

## 📝 启动日志

如果启动失败，请检查：

1. **WSL窗口中的错误信息**
   - 启动脚本会在新窗口中显示错误
   - 查看具体的错误消息

2. **Streamlit日志**
   - 在WSL中查看: `~/.streamlit/logs/`

3. **系统日志**
   - Windows事件查看器
   - WSL日志: `wsl --list --verbose`

---

## 🆘 获取帮助

如果以上方法都无法解决问题：

1. 运行诊断工具: `python diagnose_flashdock.py`
2. 查看错误日志
3. 检查WSL环境配置
4. 确认所有依赖已正确安装

---

## ✅ 已修复的问题

- ✅ 环境名称不一致 (`flash_dock_wsl` → `flash_dock`)
- ✅ 创建了改进的启动脚本 (`start_flashdock_fixed.bat`)
- ✅ 添加了更详细的错误检查
- ✅ 改进了服务状态验证

---

**最后更新**: 2025-01-XX
**状态**: 环境名称问题已修复 ✅

