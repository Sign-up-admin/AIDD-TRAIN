# FlashDock 测试快速开始指南

## 快速开始

### 1. 安装测试依赖

```bash
pip install -r requirements-dev.txt
```

### 2. 运行所有测试

```bash
pytest tests/flashdock/ -v
```

### 3. 运行特定类型测试

```bash
# 单元测试（最快，不需要服务运行）
pytest tests/flashdock/unit/ -v

# 集成测试（需要mock服务）
pytest tests/flashdock/integration/ -v

# 端到端测试（需要实际服务运行）
pytest tests/flashdock/e2e/ -v
```

## 测试示例

### 运行单个测试文件

```bash
# 测试负载均衡器
pytest tests/flashdock/unit/test_load_balancer.py -v

# 测试CompassClient
pytest tests/flashdock/unit/test_compass_client.py -v
```

### 运行单个测试用例

```bash
# 测试轮询策略
pytest tests/flashdock/unit/test_load_balancer.py::TestLoadBalancer::test_select_service_round_robin -v
```

## 查看测试覆盖率

```bash
# 生成HTML覆盖率报告
pytest tests/flashdock/ --cov=FLASH_DOCK-main/services --cov-report=html

# 查看报告
# 打开 htmlcov/index.html
```

## 常见命令

```bash
# 运行测试并显示输出
pytest tests/flashdock/ -v -s

# 运行测试并在第一个失败时停止
pytest tests/flashdock/ -x

# 运行上次失败的测试
pytest tests/flashdock/ --lf

# 运行标记为slow的测试
pytest tests/flashdock/ -m slow

# 跳过标记为slow的测试
pytest tests/flashdock/ -m "not slow"
```

## 需要帮助？

查看完整文档：`tests/flashdock/README.md`

