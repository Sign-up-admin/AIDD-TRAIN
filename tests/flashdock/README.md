# FlashDock 自动化测试套件

## 概述

本测试套件为FlashDock应用提供全面的自动化测试，包括单元测试、集成测试和端到端测试。

## 测试结构

```
tests/flashdock/
├── unit/                    # 单元测试
│   ├── test_compass_client.py
│   ├── test_registry_client.py
│   ├── test_service_manager.py
│   ├── test_load_balancer.py
│   ├── test_training_management.py
│   ├── test_data_management.py
│   └── test_service_monitor.py
├── integration/             # 集成测试
│   ├── test_compass_integration.py
│   ├── test_registry_integration.py
│   └── test_training_lifecycle.py
├── e2e/                     # 端到端测试
│   ├── test_flashdock_ui.py
│   └── conftest_e2e.py
├── fixtures/                # 测试fixtures和mock数据
│   ├── mock_services.py
│   └── test_data/
├── conftest.py              # 测试配置
└── README.md                # 本文档
```

## 安装依赖

### 基础测试依赖

```bash
pip install -r requirements-dev.txt
```

### 端到端测试依赖（可选）

如果需要运行浏览器自动化测试：

```bash
pip install selenium
# 或
pip install playwright
playwright install
```

## 运行测试

### 运行所有测试

```bash
# 运行所有FlashDock测试
pytest tests/flashdock/ -v

# 运行所有测试并生成覆盖率报告
pytest tests/flashdock/ --cov=FLASH_DOCK-main --cov-report=html --cov-report=term
```

### 运行特定类型的测试

```bash
# 只运行单元测试
pytest tests/flashdock/unit/ -v -m unit

# 只运行集成测试
pytest tests/flashdock/integration/ -v -m integration

# 只运行端到端测试
pytest tests/flashdock/e2e/ -v -m e2e
```

### 运行特定测试文件

```bash
# 运行CompassClient测试
pytest tests/flashdock/unit/test_compass_client.py -v

# 运行训练任务生命周期测试
pytest tests/flashdock/integration/test_training_lifecycle.py -v
```

### 运行特定测试用例

```bash
# 运行特定测试类
pytest tests/flashdock/unit/test_load_balancer.py::TestLoadBalancer -v

# 运行特定测试方法
pytest tests/flashdock/unit/test_load_balancer.py::TestLoadBalancer::test_select_service_round_robin -v
```

## 测试标记

测试使用pytest标记进行分类：

- `@pytest.mark.unit` - 单元测试
- `@pytest.mark.integration` - 集成测试
- `@pytest.mark.e2e` - 端到端测试
- `@pytest.mark.slow` - 运行较慢的测试
- `@pytest.mark.flashdock` - FlashDock特定测试

### 使用标记运行测试

```bash
# 运行所有单元测试
pytest -m unit

# 运行所有集成测试
pytest -m integration

# 运行所有端到端测试
pytest -m e2e

# 排除慢速测试
pytest -m "not slow"
```

## 测试覆盖范围

### 单元测试

#### 服务客户端模块
- **LoadBalancer**: 负载均衡策略（轮询、随机、最少连接）
- **FlashDockRegistryClient**: 服务发现和健康检查
- **ServiceManager**: 服务管理和缓存
- **CompassClient**: COMPASS服务API客户端
  - 训练任务管理（创建、启动、停止、查询）
  - 数据集管理（上传、列表、删除）
  - 模型管理（列表、查询）
  - 推理服务状态

#### 页面模块
- **training_management**: 训练管理页面辅助函数
- **data_management**: 数据管理页面
- **service_monitor**: 服务监控页面

### 集成测试

- **服务发现集成**: 测试服务注册中心和服务发现流程
- **训练任务生命周期**: 测试完整的训练任务流程
- **API集成**: 测试COMPASS服务API的集成

### 端到端测试

- **服务可用性**: 检查所有服务是否运行
- **UI测试**: 测试Streamlit应用的用户界面（需要浏览器自动化工具）

## 测试配置

### pytest配置

测试配置在 `pytest.ini` 中：

```ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
markers =
    unit: Unit tests
    integration: Integration tests
    slow: Slow running tests
    e2e: End-to-end tests
    flashdock: FlashDock specific tests
```

### 环境变量

测试可以使用以下环境变量：

- `COMPASS_HOST`: COMPASS服务主机（默认: 127.0.0.1）
- `COMPASS_PORT`: COMPASS服务端口（默认: 8080）
- `REGISTRY_URL`: 服务注册中心URL（默认: http://localhost:8500）
- `FLASHDOCK_URL`: FlashDock前端URL（默认: http://localhost:8501）

## 编写新测试

### 单元测试示例

```python
import pytest
from unittest.mock import Mock, patch

class TestMyComponent:
    def test_basic_functionality(self):
        # 测试基本功能
        assert True
    
    @patch('module.external_dependency')
    def test_with_mock(self, mock_dep):
        # 使用mock测试
        mock_dep.return_value = "mocked"
        result = my_function()
        assert result == "expected"
```

### 集成测试示例

```python
import pytest
import requests_mock

@pytest.mark.integration
class TestMyIntegration:
    def test_api_integration(self):
        with requests_mock.Mocker() as m:
            m.get("http://api.example.com/data", json={"result": "ok"})
            # 测试集成逻辑
            assert True
```

### 端到端测试示例

```python
import pytest
import requests

@pytest.mark.e2e
class TestMyE2E:
    def test_service_available(self, service_url):
        response = requests.get(service_url, timeout=5)
        assert response.status_code == 200
```

## 持续集成

### GitHub Actions示例

```yaml
name: FlashDock Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.8'
      - name: Install dependencies
        run: |
          pip install -r requirements-dev.txt
      - name: Run unit tests
        run: pytest tests/flashdock/unit/ -v
      - name: Run integration tests
        run: pytest tests/flashdock/integration/ -v
```

## 故障排除

### 常见问题

1. **导入错误**
   - 确保项目根目录在Python路径中
   - 检查 `sys.path` 设置

2. **服务连接错误**
   - 确保相关服务正在运行
   - 检查服务URL和端口配置

3. **Mock不工作**
   - 确保正确使用 `@patch` 装饰器
   - 检查patch路径是否正确

4. **测试超时**
   - 增加超时时间
   - 检查网络连接
   - 使用 `@pytest.mark.slow` 标记慢速测试

## 贡献指南

1. 为新功能编写测试
2. 确保所有测试通过
3. 保持测试代码简洁和可读
4. 添加适当的文档字符串
5. 使用有意义的测试名称

## 参考资源

- [pytest文档](https://docs.pytest.org/)
- [unittest.mock文档](https://docs.python.org/3/library/unittest.mock.html)
- [requests-mock文档](https://requests-mock.readthedocs.io/)

