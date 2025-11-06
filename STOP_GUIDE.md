# 停止服务指南

## 快速停止所有服务

### 方法1：使用停止脚本（推荐）

双击运行 `stop_all.bat` 文件，脚本会自动停止所有服务。

### 方法2：手动停止

#### 方法2a：使用任务管理器
1. 打开任务管理器（Ctrl + Shift + Esc）
2. 找到以下进程并结束：
   - `python.exe` (运行服务注册中心或COMPASS服务)
   - `python.exe` (运行FLASH-DOCK前端)
   - 或直接关闭对应的命令行窗口

#### 方法2b：使用命令行

**停止特定端口的服务：**
```batch
# 查找占用端口的进程
netstat -ano | findstr ":8500"
netstat -ano | findstr ":8080"
netstat -ano | findstr ":8501"

# 停止进程（替换 PID 为实际进程ID）
taskkill /F /PID <PID>
```

**或使用PowerShell：**
```powershell
# 停止所有相关进程
Get-Process | Where-Object {$_.Path -like "*AIDDTRAIN*" -or $_.Path -like "*flash_dock*"} | Stop-Process -Force
```

### 方法3：关闭窗口

最简单的办法是直接关闭运行服务的命令行窗口。

## 验证服务已停止

### 检查端口状态
```batch
netstat -ano | findstr ":8500 :8080 :8501"
```

如果没有输出，说明所有服务已停止。

### 检查健康检查端点
```batch
curl http://localhost:8500/health
curl http://localhost:8080/health
curl http://localhost:8501
```

如果无法连接，说明服务已停止。

## 重新启动服务

停止服务后，如需重新启动：

```batch
# 启动所有服务
start_all.bat

# 或重启服务
restart_services.bat
```

## 注意事项

1. **强制停止**：使用 `taskkill /F` 会强制终止进程，可能导致数据丢失
2. **优雅停止**：建议先尝试正常关闭命令行窗口，让服务有机会清理资源
3. **正在运行的任务**：如果有关键任务正在运行，停止前请确保已保存数据


