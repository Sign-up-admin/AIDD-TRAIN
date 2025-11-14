# 前端代码检查报告目录

此目录用于存放前端代码质量检查报告。

## 生成报告

运行以下命令生成检查报告：

```bash
# 1. 提取前端代码
python scripts/extract_frontend_code.py FLASH_DOCK-main

# 2. 运行检查（需要先安装 Node.js 和 npm 依赖）
python scripts/check_frontend.py
```

## 前置要求

在运行检查之前，请确保：

1. **已安装 Node.js** (>= 14.0.0)
   - 下载地址: https://nodejs.org/
   - 验证: `node --version`

2. **已安装 npm 依赖**
   ```bash
   npm install
   ```

## 报告格式

检查完成后，将生成以下文件：

- `frontend_lint_report_YYYYMMDD_HHMMSS.md` - Markdown 格式报告
- `frontend_lint_report_YYYYMMDD_HHMMSS.json` - JSON 格式报告

报告包含：
- 检查摘要（总计、错误、警告、通过）
- 按工具分组的问题列表
- 问题位置和修复建议

## 更多信息

详细使用说明请参考：[Frontend Code Quality Guide](../docs/FRONTEND_CODE_QUALITY.md)



