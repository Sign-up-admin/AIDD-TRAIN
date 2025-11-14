# COMPASS 配置指南

本文档详细说明 COMPASS 服务的配置选项。

## 配置文件

COMPASS 服务使用环境变量进行配置。你可以：

1. 创建 `.env` 文件（开发环境）
2. 创建 `.env.production` 文件（生产环境）
3. 直接设置环境变量

参考 `.env.example` 和 `.env.production.example` 获取配置模板。

## 配置验证

使用配置验证工具检查配置：

```bash
python scripts/validate_config.py
```

## 配置分类

### 1. 基础配置

#### ENVIRONMENT
- **描述**: 运行环境
- **类型**: 字符串
- **可选值**: `development`, `production`
- **默认值**: `development`
- **示例**: `ENVIRONMENT=production`

#### COMPASS_HOST
- **描述**: 服务监听地址
- **类型**: 字符串
- **默认值**: `0.0.0.0`
- **示例**: `COMPASS_HOST=0.0.0.0`

#### COMPASS_PORT
- **描述**: 服务监听端口
- **类型**: 整数 (1-65535)
- **默认值**: `8080`
- **示例**: `COMPASS_PORT=8080`

### 2. 目录配置

#### COMPASS_DATA_DIR
- **描述**: 数据文件目录
- **类型**: 路径字符串
- **默认值**: `./data`
- **示例**: `COMPASS_DATA_DIR=/var/lib/compass/data`

#### COMPASS_CHECKPOINT_DIR
- **描述**: 模型检查点目录
- **类型**: 路径字符串
- **默认值**: `./checkpoints`
- **示例**: `COMPASS_CHECKPOINT_DIR=/var/lib/compass/checkpoints`

#### COMPASS_LOG_DIR
- **描述**: 日志文件目录
- **类型**: 路径字符串
- **默认值**: `./logs`
- **示例**: `COMPASS_LOG_DIR=/var/log/compass`

### 3. 安全配置

#### AUTH_ENABLED
- **描述**: 是否启用API密钥认证
- **类型**: 布尔值 (`true`/`false`)
- **默认值**: `false`
- **生产环境**: 必须设置为 `true`
- **示例**: `AUTH_ENABLED=true`

#### API_KEY
- **描述**: API密钥（单密钥）
- **类型**: 字符串
- **默认值**: 空
- **注意**: 与 `API_KEYS` 二选一
- **示例**: `API_KEY=your-secure-api-key-here`

#### API_KEYS
- **描述**: 多个API密钥（逗号分隔，用于密钥轮换）
- **类型**: 字符串
- **默认值**: 空
- **注意**: 与 `API_KEY` 二选一
- **示例**: `API_KEYS=key1,key2,key3`

#### FORCE_AUTH_CRITICAL
- **描述**: 生产环境是否强制关键端点认证
- **类型**: 布尔值
- **默认值**: `true`
- **示例**: `FORCE_AUTH_CRITICAL=true`

#### CORS_ORIGINS
- **描述**: 允许的CORS来源（逗号分隔）
- **类型**: 字符串
- **默认值**: `http://localhost:8501,http://127.0.0.1:8501,http://localhost:3000,http://127.0.0.1:3000`
- **注意**: 生产环境必须设置，不能使用通配符 `*`
- **示例**: `CORS_ORIGINS=https://yourdomain.com,https://api.yourdomain.com`

### 4. 速率限制配置

#### RATE_LIMIT_DEFAULT
- **描述**: 默认速率限制（每分钟请求数）
- **类型**: 整数
- **默认值**: `100`
- **示例**: `RATE_LIMIT_DEFAULT=100`

#### RATE_LIMIT_TRAINING
- **描述**: 训练端点速率限制（每分钟请求数）
- **类型**: 整数
- **默认值**: `10`
- **示例**: `RATE_LIMIT_TRAINING=10`

#### RATE_LIMIT_UPLOAD
- **描述**: 上传端点速率限制（每分钟请求数）
- **类型**: 整数
- **默认值**: `3`
- **示例**: `RATE_LIMIT_UPLOAD=3`

#### RATE_LIMIT_INFERENCE
- **描述**: 推理端点速率限制（每分钟请求数）
- **类型**: 整数
- **默认值**: `20`
- **示例**: `RATE_LIMIT_INFERENCE=20`

### 5. 资源限制配置

#### COMPASS_MAX_WORKERS
- **描述**: 最大工作线程数
- **类型**: 整数
- **默认值**: `4`
- **示例**: `COMPASS_MAX_WORKERS=8`

#### MAX_CONCURRENT_UPLOADS
- **描述**: 最大并发上传数
- **类型**: 整数
- **默认值**: `2`
- **示例**: `MAX_CONCURRENT_UPLOADS=2`

#### MAX_CONCURRENT_TASKS
- **描述**: 最大并发训练任务数
- **类型**: 整数
- **默认值**: `4`
- **示例**: `MAX_CONCURRENT_TASKS=4`

#### COMPASS_UPLOAD_MAX_SIZE
- **描述**: 最大上传文件大小（字节）
- **类型**: 整数
- **默认值**: `10737418240` (10GB)
- **示例**: `COMPASS_UPLOAD_MAX_SIZE=10737418240`

### 6. 数据库配置

#### DB_CONNECTION_TIMEOUT
- **描述**: 数据库连接超时（秒）
- **类型**: 浮点数
- **默认值**: `10.0`
- **示例**: `DB_CONNECTION_TIMEOUT=10.0`

#### DB_BUSY_TIMEOUT
- **描述**: 数据库忙碌超时（毫秒）
- **类型**: 整数
- **默认值**: `5000`
- **示例**: `DB_BUSY_TIMEOUT=5000`

#### DB_CACHE_SIZE
- **描述**: 数据库缓存大小（负数表示KB）
- **类型**: 整数
- **默认值**: `-2000` (2000KB)
- **示例**: `DB_CACHE_SIZE=-2000`

### 7. 服务注册配置

#### REGISTRY_URL
- **描述**: 服务注册中心URL
- **类型**: URL字符串
- **默认值**: `http://localhost:8500`
- **示例**: `REGISTRY_URL=http://registry:8500`

#### REGISTRY_DB_PATH
- **描述**: 注册中心数据库路径
- **类型**: 路径字符串
- **默认值**: `./registry.db`
- **示例**: `REGISTRY_DB_PATH=/var/lib/compass/registry.db`

### 8. 日志配置

#### LOG_LEVEL
- **描述**: 日志级别
- **类型**: 字符串
- **可选值**: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`
- **默认值**: `INFO`
- **示例**: `LOG_LEVEL=INFO`

## 配置示例

### 开发环境配置

```bash
ENVIRONMENT=development
COMPASS_HOST=0.0.0.0
COMPASS_PORT=8080
AUTH_ENABLED=false
CORS_ORIGINS=http://localhost:8501,http://127.0.0.1:8501
```

### 生产环境配置

```bash
ENVIRONMENT=production
COMPASS_HOST=0.0.0.0
COMPASS_PORT=8080
AUTH_ENABLED=true
API_KEY=your-secure-api-key-here
FORCE_AUTH_CRITICAL=true
CORS_ORIGINS=https://yourdomain.com
COMPASS_MAX_WORKERS=8
LOG_LEVEL=INFO
```

## 配置验证清单

### 开发环境
- [ ] `COMPASS_HOST` 和 `COMPASS_PORT` 已设置
- [ ] `AUTH_ENABLED=false`（可选）

### 生产环境
- [ ] `ENVIRONMENT=production`
- [ ] `AUTH_ENABLED=true`
- [ ] `API_KEY` 或 `API_KEYS` 已设置
- [ ] `FORCE_AUTH_CRITICAL=true`
- [ ] `CORS_ORIGINS` 已设置（不使用通配符）
- [ ] 所有目录路径已设置且可写
- [ ] 日志级别设置为 `INFO` 或更高

## 故障排查

### 配置未生效
1. 检查环境变量是否正确设置
2. 运行 `python scripts/validate_config.py` 验证配置
3. 检查服务日志中的配置相关错误

### 认证失败
1. 确认 `AUTH_ENABLED=true`
2. 确认 `API_KEY` 或 `API_KEYS` 已设置
3. 检查请求头中的 `X-API-Key` 或 `Authorization: Bearer <key>`

### CORS错误
1. 确认 `CORS_ORIGINS` 包含前端域名
2. 确认不使用通配符 `*`
3. 检查协议（http/https）是否匹配




