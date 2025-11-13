# COMPASS 监控系统设置指南

本文档说明如何设置和使用 COMPASS 服务的监控系统（Prometheus + Grafana）。

## 前置要求

- Docker 和 Docker Compose 已安装
- COMPASS 服务正在运行（端口 8080）

## 快速开始

1. **启动监控服务**：
   ```bash
   docker-compose -f docker-compose.monitoring.yml up -d
   ```

2. **访问 Grafana**：
   - URL: http://localhost:3000
   - 用户名: `admin`
   - 密码: `admin`

3. **访问 Prometheus**：
   - URL: http://localhost:9090

## 配置说明

### Prometheus 配置

配置文件位于 `monitoring/prometheus/prometheus.yml`：

- **scrape_interval**: 指标采集间隔（默认 15 秒）
- **targets**: COMPASS 服务地址（默认 `host.docker.internal:8080`）

### Grafana 配置

- **数据源**: 自动配置为 Prometheus（http://prometheus:9090）
- **仪表板**: 自动加载 `monitoring/grafana/dashboards/` 中的仪表板

## 指标说明

COMPASS 服务导出以下 Prometheus 指标：

- `compass_http_requests_total`: HTTP 请求总数（按方法、端点、状态码）
- `compass_http_request_duration_seconds`: HTTP 请求持续时间（直方图）
- `compass_http_errors_total`: HTTP 错误总数（按方法、端点、状态码）
- `compass_active_requests`: 当前活跃请求数

## 告警规则

告警规则定义在 `monitoring/prometheus/alerts.yml`：

- **HighErrorRate**: 错误率超过 10% 持续 5 分钟
- **HighResponseTime**: 95 百分位响应时间超过 5 秒持续 5 分钟
- **ServiceDown**: 服务不可用超过 1 分钟

## 访问指标

### JSON 格式（默认）
```bash
curl http://localhost:8080/metrics
```

### Prometheus 格式
```bash
curl http://localhost:8080/metrics?format=prometheus
```

## 故障排查

### Prometheus 无法采集指标

1. 检查 COMPASS 服务是否运行：
   ```bash
   curl http://localhost:8080/health
   ```

2. 检查 Prometheus 目标状态：
   - 访问 http://localhost:9090/targets
   - 确认 `compass-service` 目标状态为 UP

3. 检查网络连接：
   - 在 Docker 容器中，使用 `host.docker.internal` 访问宿主机服务
   - 如果使用 Linux，可能需要使用 `host.docker.internal:8080` 或实际 IP

### Grafana 无法显示数据

1. 检查数据源连接：
   - 登录 Grafana
   - 进入 Configuration > Data Sources
   - 测试 Prometheus 数据源连接

2. 检查仪表板配置：
   - 确认仪表板 JSON 文件在 `monitoring/grafana/dashboards/` 目录
   - 检查 Grafana 日志

## 停止监控服务

```bash
docker-compose -f docker-compose.monitoring.yml down
```

## 数据持久化

监控数据存储在 Docker volumes 中：
- `prometheus-data`: Prometheus 时序数据
- `grafana-data`: Grafana 配置和仪表板

要完全删除数据：
```bash
docker-compose -f docker-compose.monitoring.yml down -v
```


