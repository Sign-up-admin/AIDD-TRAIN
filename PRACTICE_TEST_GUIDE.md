# 实践测试指南

本指南将帮助您验证COMPASS服务的所有改进功能。

## 快速开始

### 1. 启动服务

**方法1：使用批处理文件（推荐）**
```bash
start_service_for_test.bat
```

**方法2：手动启动**
```bash
cd E:\Qinchaojun\AIDD-TRAIN
set PYTHONPATH=%CD%
python compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
```

**方法3：使用项目启动脚本**
```bash
start_compass_safe.bat
```

### 2. 运行测试脚本

在另一个终端窗口中运行：
```bash
python test_service_practice.py
```

## 测试内容

测试脚本将验证以下功能：

### 1. 健康检查
- 验证服务是否正常运行
- 检查服务状态和版本信息

### 2. CORS配置
- 验证CORS头是否正确设置
- 检查跨域请求是否被正确处理

### 3. 安全头
- 验证安全响应头：
  - `X-Content-Type-Options: nosniff`
  - `X-Frame-Options: DENY`
  - `Content-Security-Policy`

### 4. 输入清理
- 测试XSS攻击防护
- 验证危险输入是否被正确清理

### 5. 速率限制
- 测试速率限制功能
- 验证是否在超过限制时返回429状态码

### 6. 错误处理
- 测试无效输入的错误处理
- 验证错误消息格式是否正确

### 7. 指标端点
- 验证`/metrics`端点是否可访问
- 检查速率限制统计信息

### 8. API文档
- 验证Swagger文档是否可访问
- 检查API文档地址：http://localhost:8080/docs

## 手动测试

### 测试健康检查
```bash
curl http://localhost:8080/health
```

### 测试CORS
```bash
curl -H "Origin: http://localhost:8501" -v http://localhost:8080/health
```

### 测试安全头
```bash
curl -I http://localhost:8080/health
```

### 测试速率限制
```bash
# 快速发送多个请求
for i in {1..15}; do curl http://localhost:8080/api/v1/training/tasks; done
```

### 测试输入验证
```bash
# 测试无效的task_id
curl http://localhost:8080/api/v1/training/tasks/invalid-id
```

### 测试指标端点
```bash
curl http://localhost:8080/metrics
```

### 访问API文档
在浏览器中打开：
```
http://localhost:8080/docs
```

## 预期结果

### 成功标志
- ✅ 所有测试通过
- ✅ 服务正常响应请求
- ✅ 安全头正确设置
- ✅ CORS配置正确
- ✅ 速率限制正常工作
- ✅ 错误处理正确

### 常见问题

**问题1：连接被拒绝**
- 确保服务正在运行
- 检查端口8080是否被占用
- 验证防火墙设置

**问题2：CORS测试失败**
- 检查环境变量`CORS_ORIGINS`设置
- 验证CORS中间件是否正确加载

**问题3：速率限制未触发**
- 速率限制可能设置较高（默认100请求/分钟）
- 可以设置环境变量降低限制进行测试：
  ```bash
  set RATE_LIMIT_TRAINING=5
  set RATE_LIMIT_DEFAULT=10
  ```

**问题4：安全头缺失**
- 检查安全头中间件是否正确加载
- 验证中间件顺序

## 环境变量配置

### 测试环境变量
```bash
# 禁用认证（用于测试）
set AUTH_ENABLED=false
set FORCE_AUTH_CRITICAL=false

# 降低速率限制（用于测试）
set RATE_LIMIT_TRAINING=5
set RATE_LIMIT_DEFAULT=10

# CORS配置
set CORS_ORIGINS=http://localhost:8501,http://127.0.0.1:8501

# 日志级别
set LOG_LEVEL=INFO
```

## 验证清单

- [ ] 服务成功启动
- [ ] 健康检查返回200
- [ ] CORS头正确设置
- [ ] 安全头全部存在
- [ ] XSS防护正常工作
- [ ] 速率限制功能正常
- [ ] 错误处理正确
- [ ] 指标端点可访问
- [ ] API文档可访问

## 下一步

测试通过后，您可以：
1. 查看API文档：http://localhost:8080/docs
2. 测试实际功能（创建任务、上传数据等）
3. 检查日志文件了解服务运行情况
4. 监控指标端点了解服务性能

## 故障排除

如果遇到问题，请检查：
1. 服务日志文件
2. 端口占用情况
3. 环境变量配置
4. 依赖项是否正确安装

更多信息请参考：
- `README.md` - 项目主文档
- `ROBUSTNESS_IMPROVEMENTS.md` - 改进总结
- `FINAL_TEST_SUMMARY.md` - 测试总结


