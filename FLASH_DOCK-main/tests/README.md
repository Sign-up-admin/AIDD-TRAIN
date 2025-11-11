# FLASH_DOCK-main 单元测试

## 测试结构

```
tests/
├── __init__.py
├── test_compass_client.py      # CompassClient 测试
├── test_load_balancer.py       # LoadBalancer 测试
└── test_registry_client.py     # FlashDockRegistryClient 测试
```

## 运行测试

### 运行所有测试
```bash
python -m pytest FLASH_DOCK-main/tests -v
```

### 运行特定测试文件
```bash
python -m pytest FLASH_DOCK-main/tests/test_compass_client.py -v
python -m pytest FLASH_DOCK-main/tests/test_load_balancer.py -v
python -m pytest FLASH_DOCK-main/tests/test_registry_client.py -v
```

### 运行特定测试用例
```bash
python -m pytest FLASH_DOCK-main/tests/test_compass_client.py::TestCompassClient::test_list_training_tasks -v
```

### 生成覆盖率报告
```bash
python -m pytest FLASH_DOCK-main/tests --cov=FLASH_DOCK-main/services --cov-report=html
```

## 测试依赖

测试需要以下依赖（已在requirements_dev.txt中）：
- pytest
- pytest-cov
- unittest (标准库)

注意：某些测试可能需要运行的服务（如服务注册中心）。这些测试使用mock来避免实际依赖。

## 测试覆盖

当前测试覆盖：
- CompassClient 基本功能
- LoadBalancer 负载均衡策略
- FlashDockRegistryClient 服务发现
- TaskStreamClient 回调处理

## 添加新测试

当添加新功能时，请添加相应的测试：

1. 在 `tests/` 目录创建新的测试文件
2. 使用 `unittest.TestCase` 或 `pytest` 编写测试
3. 使用 `unittest.mock` 模拟外部依赖
4. 确保测试独立且可重复运行

## 示例

```python
import unittest
from unittest.mock import Mock, patch

class TestMyFeature(unittest.TestCase):
    def setUp(self):
        """设置测试环境"""
        pass
    
    def test_feature_behavior(self):
        """测试功能行为"""
        # 测试代码
        self.assertEqual(expected, actual)
```











