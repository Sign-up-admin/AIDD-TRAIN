# 手动启动服务指南

如果自动启动脚本无法正常工作，请按照以下步骤手动启动服务：

## 方法1：使用独立的批处理文件（推荐）

### 步骤1：启动服务注册中心

双击运行 `start_registry_direct.bat`，或者打开命令行执行：

```batch
cd E:\Qinchaojun\AIDD-TRAIN
C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe services/registry/server.py --host 0.0.0.0 --port 8500
```

### 步骤2：启动 COMPASS 服务

打开**新的命令行窗口**，双击运行 `start_compass_direct.bat`，或者执行：

```batch
cd E:\Qinchaojun\AIDD-TRAIN
C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
```

## 方法2：使用 Anaconda Prompt

### 步骤1：打开 Anaconda Prompt

1. 点击开始菜单
2. 搜索 "Anaconda Prompt"
3. 以管理员身份运行（可选）

### 步骤2：激活环境并启动服务注册中心

```bash
conda activate AIDDTRAIN
cd E:\Qinchaojun\AIDD-TRAIN
python services/registry/server.py --host 0.0.0.0 --port 8500
```

### 步骤3：在新窗口启动 COMPASS 服务

打开**新的 Anaconda Prompt 窗口**：

```bash
conda activate AIDDTRAIN
cd E:\Qinchaojun\AIDD-TRAIN
python compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
```

## 方法3：使用 conda run（无需激活）

### 步骤1：启动服务注册中心

```batch
cd E:\Qinchaojun\AIDD-TRAIN
conda run -n AIDDTRAIN python services/registry/server.py --host 0.0.0.0 --port 8500
```

### 步骤2：启动 COMPASS 服务（新窗口）

```batch
cd E:\Qinchaojun\AIDD-TRAIN
conda run -n AIDDTRAIN python compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
```

## 验证服务

启动后，在浏览器中访问：

- 服务注册中心: http://localhost:8500/health
- COMPASS 服务: http://localhost:8080/health
- API 文档: http://localhost:8080/docs

## 常见问题排查

### 1. 端口被占用

检查端口占用：
```batch
netstat -ano | findstr ":8500"
netstat -ano | findstr ":8080"
```

如果端口被占用，可以：
- 修改启动命令中的端口号
- 或者停止占用端口的进程

### 2. 模块导入错误

如果提示找不到模块，请安装依赖：
```batch
C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe -m pip install -r requirements_service.txt
```

### 3. 查看详细错误信息

如果服务无法启动，查看命令行窗口中的错误信息，通常会有详细的错误提示。

## 停止服务

直接关闭对应的命令行窗口即可停止服务。

