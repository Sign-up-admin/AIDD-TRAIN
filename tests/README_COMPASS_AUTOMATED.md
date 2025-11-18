# COMPASS自动化测试指南

本文档介绍COMPASS自动化测试套件的使用方法和测试结构。

## 目录结构

```
tests/
├── test_compass_automated.py          # 主测试套件
├── test_compass_service_integration.py # 服务集成测试
├── test_compass_training_unit.py      # 训练模块单元测试
├── test_compass_data_processing.py    # 数据处理测试
├── conftest.py                        # pytest配置和fixtures
└── README_COMPASS_AUTOMATED.md        # 本文档
```

## 快速开始

### 运行所有测试

```bash
# 使用pytest直接运行
pytest tests/ -v

# 使用测试运行器
python run_compass_tests.py

# 运行特定类型测试
python run_compass_tests.py --type unit
python run_compass_tests.py --type integration
python run_compass_tests.py --type e2e
```

### 运行特定测试文件

```bash
# 使用pytest
pytest tests/test_compass_automated.py -v

# 使用测试运行器
python run_compass_tests.py --test test_compass_automated.py
```

### 生成覆盖率报告

```bash
# 自动生成覆盖率报告
python run_compass_tests.py --type all

# 查看HTML覆盖率报告
# 打开 reports/test_reports/htmlcov/index.html
```

## 测试类型

### 1. 单元测试 (Unit Tests)

测试单个组件和函数的功能。

**标记**: `@pytest.mark.unit`

**测试范围**:
- 配置加载和验证
- 数据模型验证
- 工具函数
- 业务逻辑组件

**运行方式**:
```bash
pytest -m unit
python run_compass_tests.py --type unit
```

### 2. 集成测试 (Integration Tests)

测试多个组件之间的交互。

**标记**: `@pytest.mark.integration`

**测试范围**:
- API端点
- 服务间通信
- 数据库操作
- 文件系统操作

**运行方式**:
```bash
pytest -m integration
python run_compass_tests.py --type integration
```

### 3. 端到端测试 (E2E Tests)

测试完整的用户工作流程。

**标记**: `@pytest.mark.e2e` 和 `@pytest.mark.slow`

**测试范围**:
- 完整的训练流程
- 端到端推理流程
- 系统集成场景

**运行方式**:
```bash
pytest -m e2e
python run_compass_tests.py --type e2e
```

## 测试套件说明

### test_compass_automated.py

主测试套件，包含：
- `CompassTestSuite` 类：测试结果收集和报告生成
- 单元测试：配置、模型、错误代码等
- 集成测试：API端点测试
- 端到端测试：完整工作流测试

### test_compass_service_integration.py

服务集成测试，包括：
- 服务启动和健康检查
- API文档验证
- 训练任务CRUD操作
- 任务停止功能
- 数据集和模型端点
- 错误处理
- CORS和速率限制

### test_compass_training_unit.py

训练模块单元测试，包括：
- 进度跟踪器
- 训练配置准备
- 任务生命周期日志
- 训练异常处理
- 检查点管理
- 模型架构
- 训练辅助函数

### test_compass_data_processing.py

数据处理测试，包括：
- 数据集初始化
- 数据加载器
- 图转换和特征提取
- 数据服务
- 文件上传处理

## 测试配置

### 环境变量

测试会自动设置以下环境变量：
- `AUTH_ENABLED=false` - 禁用认证
- `FORCE_AUTH_CRITICAL=false` - 禁用强制认证
- `RATE_LIMIT_TRAINING=10000` - 提高速率限制
- `RATE_LIMIT_DEFAULT=10000` - 提高默认速率限制

### Fixtures

在 `conftest.py` 中定义的fixtures：
- `temp_dir` - 临时目录
- `mock_service_config` - Mock服务配置
- `sample_training_config` - 示例训练配置
- `sample_inference_request` - 示例推理请求
- `client` - FastAPI测试客户端

## 测试报告

### 控制台输出

测试运行时会显示：
- 测试执行进度
- 通过/失败的测试列表
- 错误信息和堆栈跟踪
- 覆盖率摘要

### JSON报告

测试结果会保存为JSON格式：
- 位置: `reports/test_reports/test_report_YYYYMMDD_HHMMSS.json`
- 包含: 测试摘要、详细结果、时间戳等

### HTML覆盖率报告

覆盖率报告以HTML格式生成：
- 位置: `reports/test_reports/htmlcov/index.html`
- 包含: 代码覆盖率、未覆盖行、覆盖率趋势等

### XML覆盖率报告

用于CI/CD集成：
- 位置: `reports/test_reports/coverage.xml`
- 格式: Cobertura XML

## 最佳实践

### 编写新测试

1. **遵循命名约定**:
   - 测试文件: `test_*.py`
   - 测试函数: `test_*`
   - 测试类: `Test*`

2. **使用适当的标记**:
   ```python
   @pytest.mark.unit
   def test_something():
       pass
   ```

3. **使用fixtures**:
   ```python
   def test_with_fixture(temp_dir, client):
       # 使用fixture
       pass
   ```

4. **Mock外部依赖**:
   ```python
   @patch('module.external_service')
   def test_with_mock(mock_service):
       mock_service.return_value = expected_value
       # 测试代码
   ```

5. **清理资源**:
   - 使用fixtures自动清理
   - 在teardown中清理临时文件

### 测试组织

- **单元测试**: 快速、独立、可重复
- **集成测试**: 测试组件交互
- **端到端测试**: 测试完整流程，可能较慢

### 性能考虑

- 单元测试应该很快（< 1秒）
- 集成测试可以稍慢（< 10秒）
- 端到端测试可能很慢（> 1分钟），使用 `@pytest.mark.slow` 标记

## CI/CD集成

### GitHub Actions示例

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
        with:
          python-version: '3.9'
      - run: pip install -r requirements.txt
      - run: pip install pytest pytest-cov pytest-json-report
      - run: python run_compass_tests.py --type unit --output json
      - uses: codecov/codecov-action@v2
        with:
          files: reports/test_reports/coverage.xml
```

## 故障排除

### 常见问题

1. **导入错误**
   - 确保项目根目录在Python路径中
   - 检查 `conftest.py` 中的路径设置

2. **测试失败但功能正常**
   - 检查Mock是否正确设置
   - 验证测试环境变量

3. **覆盖率报告不准确**
   - 确保使用 `--cov=compass` 参数
   - 检查 `.coveragerc` 配置

4. **测试运行缓慢**
   - 使用 `-m unit` 只运行单元测试
   - 检查是否有不必要的端到端测试

## 贡献指南

添加新测试时：
1. 遵循现有代码风格
2. 添加适当的文档字符串
3. 使用有意义的测试名称
4. 确保测试独立且可重复
5. 更新本文档（如需要）

## 参考资源

- [pytest文档](https://docs.pytest.org/)
- [FastAPI测试文档](https://fastapi.tiangolo.com/tutorial/testing/)
- [Python unittest.mock文档](https://docs.python.org/3/library/unittest.mock.html)

