# Compass终端显示问题排查指南

## 问题：无法显示Compass终端

### 重要说明

**关于"websockets library not available"警告**：
- 这个警告**不影响**前端终端显示
- 前端终端使用**浏览器原生WebSocket API**，不依赖Python的websockets库
- 该警告只影响Python后端的WebSocket客户端（TaskStreamClient），前端不受影响

### 排查步骤

#### 1. 检查任务状态
终端只在以下状态时显示：
- ✅ `running` - 任务正在运行
- ✅ `initializing` - 任务正在初始化
- ❌ 其他状态（如 `pending`, `paused`, `completed` 等）不会显示终端

**解决方法**：
- 确保任务状态为 `running` 或 `initializing`
- 如果任务状态不对，点击"启动任务"按钮

#### 2. 检查COMPASS服务状态
确保COMPASS服务已启动并注册到服务注册中心：
- 检查服务注册中心是否运行（默认端口8500）
- 检查COMPASS服务是否运行（默认端口8080）
- 在"服务监控"页面查看服务状态

#### 3. 检查WebSocket端点
终端需要连接到：`ws://<compass_host>:<compass_port>/api/v1/training/tasks/<task_id>/stream`

**调试方法**：
- 在终端显示区域，展开"🔧 调试信息"查看WebSocket URL
- 检查URL是否正确
- 尝试在浏览器中直接访问COMPASS服务的健康检查端点

#### 4. 检查浏览器控制台
打开浏览器开发者工具（F12），查看：
- Console标签：是否有WebSocket连接错误
- Network标签：WebSocket连接是否建立（状态码101）

#### 5. 检查COMPASS服务日志
查看COMPASS服务的日志，确认：
- WebSocket端点是否已注册
- 是否有连接请求到达
- 是否有错误信息

### 常见问题

#### Q: 看到"websockets library not available"警告，终端不显示
**A**: 这个警告不影响前端终端。检查：
1. 任务状态是否为 `running` 或 `initializing`
2. COMPASS服务是否正常运行
3. 浏览器控制台是否有错误

#### Q: 任务状态是running，但终端不显示
**A**: 
1. 检查"启用实时终端转播"复选框是否已勾选
2. 查看"🔧 调试信息"中的WebSocket URL是否正确
3. 检查浏览器控制台的错误信息

#### Q: 终端显示但连接失败
**A**:
1. 检查COMPASS服务的WebSocket端点是否正常
2. 检查防火墙/网络设置
3. 查看COMPASS服务日志

### 技术细节

#### 前端终端实现
- 使用 **xterm.js** 终端模拟器
- 使用 **浏览器原生WebSocket API** 连接
- 不依赖Python的websockets库

#### 后端WebSocket端点
- 路径：`/api/v1/training/tasks/{task_id}/stream`
- 协议：WebSocket (ws://)
- 消息格式：JSON
  ```json
  {
    "type": "log" | "resources" | "connected" | "error",
    "data": "..."
  }
  ```

### 验证步骤

1. ✅ 确认COMPASS服务运行：`http://localhost:8080/health`
2. ✅ 确认服务已注册：查看"服务监控"页面
3. ✅ 创建并启动训练任务
4. ✅ 等待任务状态变为 `running`
5. ✅ 勾选"启用实时终端转播"
6. ✅ 查看终端是否显示
7. ✅ 检查浏览器控制台是否有错误

### 如果仍然无法解决

1. 查看COMPASS服务日志
2. 查看浏览器控制台错误
3. 检查网络连接
4. 验证WebSocket端点是否可访问











