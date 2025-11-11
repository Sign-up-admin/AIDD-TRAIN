# COMPASS微服务架构使用指南

## 概述

本项目实现了COMPASS的微服务化架构，将COMPASS从单一脚本转换为可分布式部署的服务，并通过服务注册中心实现服务发现和负载均衡。

## 架构组件

1. **服务注册中心** (`services/registry/`) - 提供服务注册、发现和健康检查
2. **COMPASS服务** (`compass/service/`) - 训练、推理、数据管理服务
3. **FLASH-DOCK客户端** (`FLASH_DOCK-main/services/`) - Web前端服务客户端

## 环境要求

- **COMPASS服务**: 运行在 `AIDDTRAIN` 环境 (`C:\ProgramData\Anaconda3\envs\AIDDTRAIN`)
- **FLASH-DOCK**: 运行在 `flash_dock` 环境
- **服务注册中心**: 可以运行在任一环境

## 安装依赖

### COMPASS服务环境 (AIDDTRAIN)

```bash
conda activate AIDDTRAIN
pip install -r requirements_service.txt
```

### FLASH-DOCK环境 (flash_dock)

```bash
conda activate flash_dock
pip install -r FLASH_DOCK-main/requirements_service.txt
```

### 服务注册中心

```bash
conda activate AIDDTRAIN  # 或任何环境
pip install -r services/registry/requirements.txt
```

## 启动服务

### 1. 启动服务注册中心

```bash
conda activate AIDDTRAIN
python services/registry/server.py --host 0.0.0.0 --port 8500
```

### 2. 启动COMPASS服务

```bash
conda activate AIDDTRAIN
python compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
```

可以启动多个COMPASS服务实例（在不同端口），实现负载均衡。

### 3. 启动FLASH-DOCK

```bash
conda activate flash_dock
streamlit run FLASH_DOCK-main/FlashDock.py
```

## 使用说明

### 预测亲和力

1. 在FLASH-DOCK中选择"预测亲和力"页面
2. 选择模型：COMPASS (ViSNet) 或 PLANET (传统)
3. 上传蛋白质PDB文件和配体SDF文件
4. 点击"开始预测"

如果COMPASS服务不可用，系统会自动回退到PLANET模型。

### 训练管理

1. 选择"训练管理"页面
2. 创建新的训练任务，配置参数
3. 查看任务列表和状态
4. 查看训练日志和指标

### 数据管理

1. 选择"数据管理"页面
2. 上传数据集（ZIP/TAR格式）
3. 查看和管理数据集列表

### 服务监控

1. 选择"服务监控"页面
2. 查看COMPASS服务状态
3. 查看可用模型列表
4. 监控推理服务状态

## API文档

COMPASS服务启动后，可以访问API文档：
- Swagger UI: http://localhost:8080/docs
- ReDoc: http://localhost:8080/redoc

## 配置

### COMPASS服务配置

环境变量：
- `COMPASS_HOST`: 服务主机地址（默认: 0.0.0.0）
- `COMPASS_PORT`: 服务端口（默认: 8080）
- `REGISTRY_URL`: 注册中心URL（默认: http://localhost:8500）
- `COMPASS_MAX_WORKERS`: 最大并发训练任务数（默认: 4）

### 服务注册中心配置

环境变量：
- `REGISTRY_CHECK_INTERVAL`: 健康检查间隔（秒，默认: 10）
- `REGISTRY_TIMEOUT`: 健康检查超时（秒，默认: 5）

## 故障排除

1. **无法连接到COMPASS服务**
   - 检查服务注册中心是否运行
   - 检查COMPASS服务是否注册
   - 查看服务监控页面

2. **预测失败**
   - 检查模型是否已训练
   - 查看COMPASS服务日志
   - 尝试使用PLANET作为备选

3. **训练任务失败**
   - 查看任务日志
   - 检查数据集是否有效
   - 确认GPU/CPU资源充足

## 开发

### 添加新API端点

在 `compass/service/routes/` 中添加新的路由文件，并在 `compass/service/server.py` 中注册。

### 扩展服务功能

在 `compass/service/services/` 中添加新的服务实现类。

## 许可证

见主项目LICENSE文件。












