# COMPASS 微服务系统启动指南

## 目录
- [快速开始](#快速开始)
- [启动方法](#启动方法)
- [服务管理](#服务管理)
- [服务地址](#服务地址)
- [故障排除](#故障排除)

---

## 快速开始

### 推荐方式：一键启动（最简单）

**启动所有服务：**
```bash
# 双击运行或命令行执行
start_all_services.bat
```

**重启服务：**
```bash
restart_services_safe.bat
```

**停止所有服务：**
```bash
stop_all.bat
```

---

## 启动方法

### 方法1：使用安全启动脚本（推荐，已修复路径问题）

#### 单独启动服务

**1. 启动服务注册中心**
```bash
start_registry_safe.bat
```

**2. 启动COMPASS服务**
```bash
start_compass_safe.bat
```

**3. 启动FLASH-DOCK前端（可选）**
```bash
cd FLASH_DOCK-main
start_flashdock_fixed.bat
```

#### 启动顺序
1. 先启动服务注册中心（端口8500）
2. 等待3-5秒后启动COMPASS服务（端口8080）
3. 最后启动FLASH-DOCK前端（端口8501）

### 方法2：使用Python脚本自动启动

**检查并启动服务：**
```bash
python check_and_start_services.py
```

这个脚本会：
- 自动检查服务状态
- 如果服务未运行，自动启动
- 显示服务状态信息

**检查服务状态：**
```bash
python check_service_status.py
```

**检查端口占用：**
```bash
python check_ports.py
# 或
check_ports.bat
```

### 方法3：手动启动（用于调试）

#### 启动服务注册中心
```bash
cd E:\Qinchaojun\AIDD-TRAIN
set PYTHONPATH=E:\Qinchaojun\AIDD-TRAIN
python services\registry\server.py --host 0.0.0.0 --port 8500
```

#### 启动COMPASS服务
```bash
cd E:\Qinchaojun\AIDD-TRAIN
set PYTHONPATH=E:\Qinchaojun\AIDD-TRAIN
python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
```

#### 使用Conda环境
```bash
conda activate AIDDTRAIN
cd E:\Qinchaojun\AIDD-TRAIN
set PYTHONPATH=%CD%
python services\registry\server.py --host 0.0.0.0 --port 8500
```

---

## 服务管理

### 启动服务

| 脚本 | 功能 | 说明 |
|------|------|------|
| `start_all_services.bat` | 启动所有服务 | 一键启动注册中心、COMPASS和FLASH-DOCK |
| `start_registry_safe.bat` | 启动注册中心 | 安全启动脚本（已修复路径问题） |
| `start_compass_safe.bat` | 启动COMPASS服务 | 安全启动脚本（已修复路径问题） |
| `check_and_start_services.py` | 检查并启动 | Python脚本，自动检查并启动缺失的服务 |

### 重启服务

| 脚本 | 功能 | 说明 |
|------|------|------|
| `restart_services_safe.bat` | 重启所有服务 | 先停止再启动（推荐） |
| `restart_services.bat` | 重启服务 | 旧版本（可能有路径问题） |

### 停止服务

| 脚本 | 功能 | 说明 |
|------|------|------|
| `stop_all.bat` | 停止所有服务 | 停止端口8500、8080、8501的服务 |

### 检查服务

| 脚本 | 功能 | 说明 |
|------|------|------|
| `check_service_status.py` | 检查服务状态 | 检查所有服务的运行状态 |
| `check_ports.py` | 检查端口占用 | 检查端口是否被占用 |
| `check_and_start_services.py` | 检查并启动 | 检查状态，自动启动缺失的服务 |

---

## 服务地址

启动成功后，可以访问以下地址：

### 核心服务
- **服务注册中心**: http://localhost:8500
  - 健康检查: http://localhost:8500/health
  
- **COMPASS服务**: http://localhost:8080
  - API文档: http://localhost:8080/docs
  - 健康检查: http://localhost:8080/health
  - 性能指标: http://localhost:8080/metrics

### 前端界面
- **FLASH-DOCK**: http://localhost:8501

---

## 故障排除

### 1. 路径错误（已修复）

**问题**: `OSError: failed to make path absolute`

**解决方案**:
- ✅ 使用 `*_safe.bat` 脚本（已修复路径问题）
- ✅ 使用 `check_and_start_services.py`（已处理路径问题）

**原因**: PYTHONPATH包含无效路径或相对路径

### 2. 端口被占用

**检查端口占用**:
```bash
python check_ports.py
```

**Windows查看占用进程**:
```bash
netstat -ano | findstr :8500
netstat -ano | findstr :8080
```

**停止占用进程**:
```bash
# 查看PID后停止
taskkill /F /PID <PID>
```

**解决方案**:
1. 停止占用端口的进程
2. 或修改服务配置使用其他端口
3. 使用 `stop_all.bat` 停止所有服务后重新启动

### 3. 服务启动失败

**检查清单**:
1. ✅ Python环境是否正确（conda环境AIDDTRAIN）
2. ✅ 依赖包是否安装完整
   ```bash
   pip install -r requirements_service.txt
   ```
3. ✅ 查看服务窗口中的错误信息
4. ✅ 确保PYTHONPATH设置正确
5. ✅ 检查日志文件（位于 `logs/` 目录）

### 4. 连接超时

**原因**:
- 服务注册中心未启动
- 防火墙阻止连接
- 端口配置错误

**解决方案**:
1. 确保服务注册中心先启动
2. 检查防火墙设置
3. 验证端口是否正确
4. 使用 `check_service_status.py` 检查服务状态

### 5. Conda环境问题

**激活conda环境**:
```bash
conda activate AIDDTRAIN
```

**检查Python路径**:
```bash
where python
```

**如果脚本找不到conda环境**:
- 修改脚本中的conda路径
- 或使用完整路径：`C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe`

---

## 服务启动顺序

正确的启动顺序非常重要：

1. **服务注册中心** (端口8500) - 必须最先启动
2. **等待3-5秒** - 让注册中心完全启动
3. **COMPASS服务** (端口8080) - 依赖注册中心
4. **FLASH-DOCK前端** (端口8501) - 可选，依赖COMPASS服务

---

## 验证服务运行

### 快速验证
```bash
python check_service_status.py
```

### 手动验证
- 浏览器访问: http://localhost:8500/health
- 浏览器访问: http://localhost:8080/health
- 浏览器访问: http://localhost:8080/docs

### Python验证
```python
import requests

# 检查注册中心
try:
    r = requests.get('http://localhost:8500/health', timeout=3)
    print(f"注册中心: {'OK' if r.status_code == 200 else 'Failed'}")
except Exception as e:
    print(f"注册中心: 不可用 - {e}")

# 检查COMPASS服务
try:
    r = requests.get('http://localhost:8080/health', timeout=3)
    print(f"COMPASS服务: {'OK' if r.status_code == 200 else 'Failed'}")
except Exception as e:
    print(f"COMPASS服务: 不可用 - {e}")
```

---

## 日志文件

服务日志位于以下位置：

- 服务注册中心: `logs/registry.log`
- COMPASS服务: `logs/compass-service.log`
- 训练任务日志: `logs/task_<task_id>/`

---

## 注意事项

1. **启动顺序**: 必须按顺序启动（注册中心 → COMPASS → FLASH-DOCK）
2. **路径设置**: 使用 `*_safe.bat` 脚本会自动设置正确的路径
3. **环境变量**: 安全脚本会自动清理无效的PYTHONPATH路径
4. **端口冲突**: 如果端口被占用，先停止占用进程或使用其他端口
5. **服务窗口**: 不要关闭服务运行窗口，关闭窗口会停止服务
6. **等待时间**: 启动服务后等待几秒让服务完全初始化

---

## 相关文档

- [STOP_GUIDE.md](STOP_GUIDE.md) - 停止服务指南
- [README_SERVICES.md](README_SERVICES.md) - 服务架构说明
- [README.md](README.md) - 项目主文档

---

## 更新日志

### 2025-01-XX
- ✅ 修复了PYTHONPATH路径问题
- ✅ 创建了安全启动脚本（`*_safe.bat`）
- ✅ 添加了服务状态检查工具
- ✅ 改进了错误处理和提示
- ✅ 统一了启动脚本和文档

---

## 获取帮助

如果遇到问题：
1. 查看服务窗口中的错误信息
2. 检查日志文件
3. 运行 `check_service_status.py` 检查服务状态
4. 运行 `check_ports.py` 检查端口占用
5. 参考本文档的故障排除部分

