# COMPASS自动化测试快速入门

## 快速开始

### 方式1: 使用Python测试运行器（推荐）

```bash
# 运行所有测试
python run_compass_tests.py

# 运行单元测试
python run_compass_tests.py --type unit

# 运行集成测试
python run_compass_tests.py --type integration

# 运行端到端测试
python run_compass_tests.py --type e2e

# 运行特定测试文件
python run_compass_tests.py --test test_compass_automated.py

# 不生成覆盖率报告（更快）
python run_compass_tests.py --no-coverage

# 查看测试摘要
python run_compass_tests.py --summary
```

### 方式2: 使用pytest直接运行

```bash
# 运行所有测试
pytest tests/ -v

# 运行单元测试
pytest tests/ -m unit -v

# 运行集成测试
pytest tests/ -m integration -v

# 运行端到端测试
pytest tests/ -m e2e -v

# 运行特定测试文件
pytest tests/test_compass_automated.py -v

# 生成覆盖率报告
pytest tests/ --cov=compass --cov-report=html
```

### 方式3: 使用脚本（Windows/Linux）

**Windows:**
```cmd
scripts\run_compass_tests.bat
scripts\run_compass_tests.bat unit
scripts\run_compass_tests.bat integration
```

**Linux/Mac:**
```bash
chmod +x scripts/run_compass_tests.sh
./scripts/run_compass_tests.sh
./scripts/run_compass_tests.sh unit
./scripts/run_compass_tests.sh integration
```

## 测试类型说明

### 单元测试 (Unit Tests)
- **速度**: 快速（< 1秒/测试）
- **范围**: 单个函数和组件
- **标记**: `@pytest.mark.unit`
- **示例**: 配置加载、数据验证、工具函数

### 集成测试 (Integration Tests)
- **速度**: 中等（< 10秒/测试）
- **范围**: 组件间交互
- **标记**: `@pytest.mark.integration`
- **示例**: API端点、服务通信、数据库操作

### 端到端测试 (E2E Tests)
- **速度**: 较慢（> 1分钟/测试）
- **范围**: 完整工作流程
- **标记**: `@pytest.mark.e2e` 和 `@pytest.mark.slow`
- **示例**: 完整训练流程、端到端推理

## 查看测试报告

### HTML覆盖率报告
```bash
# 生成报告后，打开浏览器访问：
reports/test_reports/htmlcov/index.html
```

### 控制台输出
测试运行时会实时显示：
- ✓ 通过的测试
- ✗ 失败的测试
- 错误信息和堆栈跟踪
- 覆盖率摘要

## 常用命令

```bash
# 只运行失败的测试
pytest --lf

# 运行上次失败的测试
pytest --ff

# 显示最慢的10个测试
pytest --durations=10

# 并行运行测试（需要pytest-xdist）
pytest -n auto

# 详细输出，显示print语句
pytest -v -s

# 只运行特定标记的测试
pytest -m "unit and not slow"
```

## 测试文件结构

```
tests/
├── test_compass_automated.py          # 主测试套件
├── test_compass_service_integration.py # 服务集成测试
├── test_compass_training_unit.py      # 训练单元测试
├── test_compass_data_processing.py    # 数据处理测试
└── conftest.py                        # pytest配置
```

## 环境要求

- Python 3.8+
- pytest
- pytest-cov (用于覆盖率)
- 项目依赖（见 requirements.txt）

安装测试依赖：
```bash
pip install pytest pytest-cov
```

## 故障排除

### 问题1: 导入错误
**解决方案**: 确保在项目根目录运行测试
```bash
cd /path/to/AIDD-TRAIN
python run_compass_tests.py
```

### 问题2: 测试失败但功能正常
**解决方案**: 检查环境变量和Mock设置
```bash
# 查看详细错误信息
pytest -v -s tests/test_compass_automated.py
```

### 问题3: 覆盖率报告不准确
**解决方案**: 确保使用正确的覆盖率参数
```bash
pytest --cov=compass --cov-report=html
```

## 下一步

- 查看 [完整测试文档](tests/README_COMPASS_AUTOMATED.md)
- 了解如何 [编写新测试](tests/README_COMPASS_AUTOMATED.md#编写新测试)
- 查看 [CI/CD集成指南](tests/README_COMPASS_AUTOMATED.md#cicd集成)

## 获取帮助

- 运行 `python run_compass_tests.py --help` 查看所有选项
- 查看 `tests/README_COMPASS_AUTOMATED.md` 获取详细文档
- 查看 `pytest.ini` 了解pytest配置

