# 快速测试指南

## 三步快速测试

### 步骤1：启动服务
```bash
start_service_for_test.bat
```
或
```bash
python compass\service_main.py --host 0.0.0.0 --port 8080
```

### 步骤2：运行测试
在另一个终端：
```bash
python test_service_practice.py
```

### 步骤3：查看结果
测试脚本会显示：
- ✅ 通过的测试
- ❌ 失败的测试
- 📊 测试统计

## 快速验证命令

```bash
# 1. 健康检查
curl http://localhost:8080/health

# 2. 查看API文档（浏览器）
start http://localhost:8080/docs

# 3. 查看指标
curl http://localhost:8080/metrics
```

## 预期输出

测试通过时应该看到：
```
[SUCCESS] 所有测试通过！服务运行正常。
总计: 8/8 测试通过
```

## 需要帮助？

查看 `PRACTICE_TEST_GUIDE.md` 获取详细说明。




