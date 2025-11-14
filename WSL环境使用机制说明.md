# WSL 环境使用机制说明

本文档详细说明如何让项目使用 WSL 中的 conda 环境。

## 核心机制

### 1. 路径映射

WSL 可以通过 `/mnt/` 目录访问 Windows 文件系统：

```
Windows 路径: E:\Qinchaojun\AIDD-TRAIN
WSL 路径:     /mnt/e/Qinchaojun/AIDD-TRAIN
```

**转换规则**：
- Windows: `E:\path\to\file` 
- WSL: `/mnt/e/path/to/file`
- 驱动器字母转为小写
- 反斜杠 `\` 转为正斜杠 `/`

### 2. 从 Windows 调用 WSL

使用 `wsl` 命令可以从 Windows 执行 WSL 中的命令：

```batch
REM Windows 批处理文件 (.bat)
wsl -d Ubuntu-24.04 bash -c "命令"
```

**关键参数**：
- `-d Ubuntu-24.04`: 指定 WSL 发行版
- `bash -c "命令"`: 在 WSL 中执行 bash 命令

### 3. 环境隔离

WSL 中的 conda 环境与 Windows 中的 conda 环境完全独立：

```
Windows 环境:
  C:\ProgramData\Anaconda3\envs\flash_dock\python.exe

WSL 环境:
  ~/miniconda3/envs/flash_dock_wsl/bin/python
```

## 实现方式

### 方式 1: 直接在 WSL 中运行（推荐）

**优点**：
- 性能最好
- 环境完全隔离
- 不需要路径转换

**使用步骤**：
```bash
# 1. 进入 WSL
wsl -d Ubuntu-24.04

# 2. 运行设置脚本
cd /mnt/e/Qinchaojun/AIDD-TRAIN
bash setup_flashdock_wsl.sh

# 3. 运行启动脚本
bash start_flashdock_wsl.sh
```

**脚本位置**：
- `setup_flashdock_wsl.sh` - 在 WSL 中运行的环境设置脚本
- `start_flashdock_wsl.sh` - 在 WSL 中运行的启动脚本

### 方式 2: 从 Windows 调用 WSL（便捷）

**优点**：
- 可以从 Windows 直接启动
- 不需要手动进入 WSL
- 适合集成到 Windows 工作流

**使用步骤**：
```batch
REM 在 Windows 中双击运行
start_flashdock_wsl.bat
```

**工作原理**：
```batch
@echo off
REM 1. 检查 WSL 是否可用
wsl --status

REM 2. 在 WSL 中激活 conda 环境
wsl -d Ubuntu-24.04 bash -c "conda activate flash_dock_wsl && ..."

REM 3. 在 WSL 中运行 Python 程序
wsl -d Ubuntu-24.04 bash -c "streamlit run FlashDock.py ..."
```

**脚本位置**：
- `start_flashdock_wsl.bat` - 从 Windows 调用 WSL 的启动脚本

## 详细工作流程

### 场景：从 Windows 启动 FlashDock（使用 WSL 环境）

```
┌─────────────────────────────────────────────────────────┐
│ Windows 系统                                            │
│                                                         │
│  1. 用户双击 start_flashdock_wsl.bat                   │
│                                                         │
│  2. 批处理脚本执行:                                     │
│     wsl -d Ubuntu-24.04 bash -c "..."                  │
│                                                         │
│  3. WSL 子系统启动                                      │
└─────────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────────┐
│ WSL (Ubuntu-24.04) 子系统                               │
│                                                         │
│  1. 初始化 conda:                                       │
│     source ~/miniconda3/etc/profile.d/conda.sh         │
│                                                         │
│  2. 激活环境:                                           │
│     conda activate flash_dock_wsl                      │
│                                                         │
│  3. 设置环境变量:                                       │
│     export PYTHONPATH=/mnt/e/Qinchaojun/AIDD-TRAIN    │
│                                                         │
│  4. 进入项目目录:                                       │
│     cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main   │
│                                                         │
│  5. 运行 Streamlit:                                     │
│     streamlit run FlashDock.py --server.port 8501     │
│                                                         │
│  6. Streamlit 服务器启动，监听 0.0.0.0:8501            │
└─────────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────────┐
│ 网络层                                                  │
│                                                         │
│  WSL2 自动端口转发:                                     │
│  localhost:8501 (Windows) ←→ localhost:8501 (WSL)     │
│                                                         │
│  浏览器访问: http://localhost:8501                     │
└─────────────────────────────────────────────────────────┘
```

## 关键技术点

### 1. Conda 环境激活

在 WSL 中，conda 需要先初始化：

```bash
# 方法 1: 使用 conda.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate flash_dock_wsl

# 方法 2: 使用 conda init（需要先运行 conda init bash）
eval "$(conda shell.bash hook)"
conda activate flash_dock_wsl
```

### 2. 路径转换

项目提供了路径转换工具 `path_converter.py`：

```python
# Windows 路径 → WSL 路径
E:\Qinchaojun\AIDD-TRAIN → /mnt/e/Qinchaojun/AIDD-TRAIN

# WSL 路径 → Windows 路径
/mnt/e/Qinchaojun/AIDD-TRAIN → E:\Qinchaojun\AIDD-TRAIN
```

### 3. 环境变量传递

在 WSL 中设置环境变量：

```bash
# 在 WSL 中
export PYTHONPATH=/mnt/e/Qinchaojun/AIDD-TRAIN:$PYTHONPATH

# 从 Windows 调用时
wsl -d Ubuntu-24.04 bash -c "export PYTHONPATH=/mnt/e/Qinchaojun/AIDD-TRAIN && python script.py"
```

### 4. 端口转发

WSL2 自动转发端口，无需额外配置：

```
WSL 中监听: 0.0.0.0:8501
Windows 访问: http://localhost:8501
```

## 文件说明

### WSL 脚本（在 WSL 中运行）

1. **`setup_flashdock_wsl.sh`**
   - 创建 conda 环境
   - 安装所有依赖
   - 验证安装

2. **`start_flashdock_wsl.sh`**
   - 激活 conda 环境
   - 启动 FlashDock
   - 在 WSL 中运行

### Windows 脚本（从 Windows 调用 WSL）

1. **`start_flashdock_wsl.bat`**
   - 检查 WSL 可用性
   - 调用 WSL 中的启动脚本
   - 从 Windows 启动

### 工具脚本

1. **`path_converter.py`**
   - Windows/WSL 路径转换
   - 自动检测环境
   - 双向转换

## 常见问题

### Q1: 为什么需要路径转换？

**A**: Windows 和 WSL 使用不同的路径格式：
- Windows: `E:\path\to\file`
- WSL: `/mnt/e/path/to/file`

WSL 通过 `/mnt/` 挂载 Windows 驱动器。

### Q2: 环境变量如何传递？

**A**: 在 `wsl` 命令中设置：

```batch
wsl -d Ubuntu-24.04 bash -c "export VAR=value && command"
```

### Q3: 如何调试 WSL 中的问题？

**A**: 
1. 直接在 WSL 中运行脚本（方式 1）
2. 查看详细输出
3. 检查 conda 环境是否正确激活

### Q4: 性能影响？

**A**: 
- 直接在 WSL 中运行：性能最好
- 从 Windows 调用：有轻微开销，但通常可忽略

### Q5: 如何同时使用 Windows 和 WSL 环境？

**A**: 它们是独立的：
- Windows: `flash_dock` 环境
- WSL: `flash_dock_wsl` 环境

可以同时运行，使用不同端口即可。

## 最佳实践

1. **开发阶段**：直接在 WSL 中运行（方式 1）
2. **生产/演示**：从 Windows 启动（方式 2）
3. **CI/CD**：使用 WSL 脚本，便于自动化

## 总结

让项目使用 WSL 中的环境的核心是：

1. **路径映射**：通过 `/mnt/` 访问 Windows 文件
2. **命令调用**：使用 `wsl` 命令从 Windows 执行 WSL 命令
3. **环境隔离**：WSL 中的 conda 环境独立于 Windows
4. **自动转发**：WSL2 自动转发网络端口

这样既可以利用 Linux 环境的优势，又保持了与 Windows 的集成。

