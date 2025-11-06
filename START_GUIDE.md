# COMPASS 服务启动指南

## 快速启动

### 方法1：使用启动脚本（推荐）

双击运行 `start_services.bat` 文件，脚本会自动：
1. 启动服务注册中心（端口 8500）
2. 启动 COMPASS 服务（端口 8080）
3. 每个服务在独立窗口中运行

### 方法2：手动启动

#### 1. 启动服务注册中心

打开命令行，执行：

```bash
conda activate AIDDTRAIN
cd E:\Qinchaojun\AIDD-TRAIN
python services/registry/server.py --host 0.0.0.0 --port 8500
```

#### 2. 启动 COMPASS 服务

打开新的命令行窗口，执行：

```bash
conda activate AIDDTRAIN
cd E:\Qinchaojun\AIDD-TRAIN
python compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
```

## 验证服务状态

### 检查服务注册中心
打开浏览器访问：http://localhost:8500/health

### 检查 COMPASS 服务
打开浏览器访问：http://localhost:8080/health

### 查看 API 文档
打开浏览器访问：http://localhost:8080/docs

## 服务地址

- **服务注册中心**: http://localhost:8500
- **COMPASS 服务**: http://localhost:8080
- **API 文档**: http://localhost:8080/docs
- **健康检查**: http://localhost:8080/health
- **性能指标**: http://localhost:8080/metrics
- **告警信息**: http://localhost:8080/alerts
- **备份列表**: http://localhost:8080/backups

## 常见问题

### 1. 端口被占用
如果端口被占用，可以：
- 修改启动脚本中的端口号
- 或者停止占用端口的进程

### 2. Conda 环境未激活
确保在执行启动命令前先激活 conda 环境：
```bash
conda activate AIDDTRAIN
```

### 3. 依赖未安装
如果遇到模块导入错误，请安装依赖：
```bash
pip install -r requirements_service.txt
```

### 4. 服务无法连接
- 检查防火墙设置
- 确认服务已成功启动（查看启动窗口的输出）
- 检查日志文件（位于 `logs/` 目录）

## 停止服务

直接关闭对应的命令行窗口即可停止服务。

## 日志查看

- 服务注册中心日志：`logs/registry.log`
- COMPASS 服务日志：`logs/compass.log`

## 下一步

服务启动后，您可以：
1. 访问 API 文档查看所有可用接口
2. 使用健康检查端点监控服务状态
3. 查看性能指标了解服务运行情况
4. 开始使用训练、推理、数据管理等功能

