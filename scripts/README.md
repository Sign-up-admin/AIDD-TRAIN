# 启动脚本说明

本目录包含所有服务启动、停止和管理脚本。

## 推荐使用的脚本（已修复路径问题）

### 启动脚本
- `start_all_services.bat` - 一键启动所有服务（推荐）
- `start_registry_safe.bat` - 启动服务注册中心（安全版本）
- `start_compass_safe.bat` - 启动COMPASS服务（安全版本）

### 管理脚本
- `restart_services_safe.bat` - 重启所有服务（推荐）
- `stop_all.bat` - 停止所有服务

### Python工具
- `check_and_start_services.py` - 检查并自动启动服务
- `check_service_status.py` - 检查服务状态
- `check_ports.py` - 检查端口占用

## 已废弃的脚本（保留用于兼容）

以下脚本可能存在问题，建议使用 `*_safe.bat` 版本：
- `start_services.bat` - 旧版本（可能有路径问题）
- `start_registry.bat` - 旧版本
- `start_compass_direct.bat` - 旧版本
- `restart_services.bat` - 旧版本（可能有路径问题）

## 脚本说明

### 批处理脚本（.bat）

#### 启动脚本
- `start_all_services.bat` - 启动所有服务（注册中心、COMPASS、FLASH-DOCK）
- `start_registry_safe.bat` - 仅启动服务注册中心
- `start_compass_safe.bat` - 仅启动COMPASS服务

#### 管理脚本
- `restart_services_safe.bat` - 停止并重启所有服务
- `stop_all.bat` - 停止所有服务

### Python脚本

#### 检查工具
- `check_service_status.py` - 检查所有服务的运行状态
- `check_ports.py` - 检查端口8500和8080的占用情况

#### 自动启动
- `check_and_start_services.py` - 检查服务状态，自动启动缺失的服务

## 使用建议

1. **首次使用**: 使用 `start_all_services.bat` 一键启动
2. **日常使用**: 使用 `restart_services_safe.bat` 重启服务
3. **故障排查**: 使用 `check_service_status.py` 检查状态
4. **端口问题**: 使用 `check_ports.py` 检查端口占用

## 注意事项

- 所有 `*_safe.bat` 脚本已修复路径问题
- 使用安全脚本可以避免 `OSError: failed to make path absolute` 错误
- Python脚本会自动处理路径和编码问题

