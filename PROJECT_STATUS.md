# COMPASS 项目启动状态

## 服务列表

### 1. 服务注册中心
- **端口**: 8500
- **状态**: 运行中
- **健康检查**: http://localhost:8500/health
- **服务列表**: http://localhost:8500/services

### 2. COMPASS 服务
- **端口**: 8080
- **状态**: 运行中
- **API 文档**: http://localhost:8080/docs
- **健康检查**: http://localhost:8080/health
- **性能指标**: http://localhost:8080/metrics
- **告警信息**: http://localhost:8080/alerts
- **备份列表**: http://localhost:8080/backups

### 3. FLASH-DOCK 前端
- **端口**: 8501
- **状态**: 运行中
- **访问地址**: http://localhost:8501

## 快速访问

### 主要界面
- **FLASH-DOCK 前端**: http://localhost:8501
- **COMPASS API 文档**: http://localhost:8080/docs

### 监控和管理
- **健康检查**: http://localhost:8080/health
- **性能指标**: http://localhost:8080/metrics
- **告警信息**: http://localhost:8080/alerts

## 启动命令

### 启动所有服务
```batch
.\start_all.bat
```

### 单独启动服务
```batch
# 服务注册中心
.\start_registry_direct.bat

# COMPASS 服务
.\start_compass_direct.bat

# FLASH-DOCK 前端
cd FLASH_DOCK-main
.\start_flashdock_fixed.bat
```

### 重启服务
```batch
.\restart_services.bat
```

## 停止服务

直接关闭对应的命令行窗口即可停止服务。

## 故障排除

如果服务无法启动：

1. 检查端口是否被占用：
   ```batch
   netstat -ano | findstr ":8500 :8080 :8501"
   ```

2. 检查 conda 环境：
   ```batch
   conda env list
   ```

3. 查看服务窗口中的错误信息

4. 检查日志文件（位于 `logs/` 目录）

