# FlashDock 端到端测试

## 概述

端到端测试用于测试FlashDock应用的完整用户流程，包括UI交互和功能验证。

## 测试类型

### 1. 基础服务可用性测试

这些测试检查服务是否运行并可访问：
- FlashDock前端服务（端口8501）
- COMPASS服务（端口8080）
- 服务注册中心（端口8500）

### 2. API端点测试

测试COMPASS服务的API端点是否可访问和响应。

### 3. 浏览器自动化测试（可选）

使用Selenium进行完整的浏览器自动化测试。需要：
- 安装selenium: `pip install selenium`
- 安装ChromeDriver或使用其他浏览器驱动

## 运行测试

### 运行所有E2E测试

```bash
pytest tests/flashdock/e2e/ -v -m e2e
```

### 运行特定测试

```bash
# 只运行服务可用性测试
pytest tests/flashdock/e2e/test_flashdock_ui.py::TestFlashDockUI -v

# 运行Selenium测试（如果已安装）
pytest tests/flashdock/e2e/test_flashdock_ui.py::TestFlashDockUISelenium -v -m slow
```

## 前置条件

1. **启动所有服务**：
   - FlashDock前端服务
   - COMPASS服务
   - 服务注册中心

2. **安装测试依赖**（可选，用于浏览器自动化）：
   ```bash
   pip install selenium
   ```

## 注意事项

- E2E测试需要实际运行的服务，确保在运行测试前启动所有服务
- 浏览器自动化测试可能需要较长时间
- 某些测试可能会被跳过（如果服务未运行或依赖未安装）

