# 前端代码质量检查工具使用指南

本文档介绍项目中配置的前端代码质量检查工具及其使用方法。

## 工具概览

项目配置了以下前端代码质量检查工具：

### HTML 检查

- **HTMLHint** - HTML 代码质量检查
  - 检查 HTML 语法错误
  - 验证标签闭合、属性格式等
  - 配置文件: `.htmlhintrc.json`

- **html5validator** - HTML5 验证（Python 包）
  - 验证 HTML5 标准合规性
  - 可通过 Python 脚本调用

### CSS 检查

- **Stylelint** - CSS/SCSS 代码检查
  - 检查 CSS 语法和最佳实践
  - 支持标准 CSS 规范
  - 配置文件: `.stylelintrc.json`

### JavaScript 检查

- **ESLint** - JavaScript 代码检查
  - 检查 JavaScript 语法错误
  - 代码风格和最佳实践
  - 配置文件: `.eslintrc.json`

### 通用格式化

- **Prettier** - 代码格式化工具
  - 自动格式化 HTML/CSS/JS/JSON
  - 统一代码风格
  - 配置文件: `.prettierrc.json`

## 安装和配置

### 方法一：使用 Docker（推荐，无需安装 Node.js）

**优点：**
- ✅ 无需在系统上安装 Node.js
- ✅ 环境隔离，不影响系统配置
- ✅ 版本统一，团队一致
- ✅ 跨平台支持

**缺点：**
- ⚠️ 需要先安装 Docker
- ⚠️ 首次构建镜像需要时间
- ⚠️ Windows 上文件挂载可能较慢

**使用步骤：**

1. **安装 Docker**
   - 下载地址: https://www.docker.com/get-started
   - Windows: 下载 Docker Desktop

2. **构建 Docker 镜像**
   ```bash
   docker build -f Dockerfile.frontend-lint -t aidd-frontend-linter:latest .
   ```

3. **运行检查**
   ```bash
   # 使用 Python 脚本（自动处理 Docker）
   python scripts/check_frontend_docker.py

   # 或直接使用 Docker
   docker run --rm -v ${PWD}:/app -w /app aidd-frontend-linter:latest npm run lint:all
   ```

### 方法二：本地安装 Node.js

**前置要求：**

1. **Node.js** (>= 14.0.0)
   - 下载地址: https://nodejs.org/
   - 验证安装: `node --version`

2. **npm** (随 Node.js 一起安装)
   - 验证安装: `npm --version`

**安装依赖：**

```bash
# 安装所有前端检查工具
npm install
```

这将安装以下工具：
- ESLint 及相关插件
- Stylelint 及相关配置
- Prettier
- HTMLHint

## 使用方法

### 1. 提取内嵌的前端代码

由于项目中的前端代码内嵌在 Python 文件中（如 Streamlit 页面），需要先提取：

```bash
# 从 FLASH_DOCK-main 目录提取代码
python scripts/extract_frontend_code.py FLASH_DOCK-main

# 指定输出目录
python scripts/extract_frontend_code.py FLASH_DOCK-main -o temp_frontend_code

# 保存文件映射到 JSON
python scripts/extract_frontend_code.py FLASH_DOCK-main -j file_map.json
```

提取的代码将保存在 `temp_frontend_code/` 目录中。

### 2. 运行代码检查

#### 方法一：使用 Docker（推荐，无需本地 Node.js）

```bash
# 使用 Docker 运行检查（自动构建镜像）
python scripts/check_frontend_docker.py

# 指定命令
python scripts/check_frontend_docker.py lint:js
python scripts/check_frontend_docker.py lint:css
python scripts/check_frontend_docker.py format

# 重新构建镜像
python scripts/check_frontend_docker.py --build
```

#### 方法二：使用 Python 脚本（需要本地 Node.js）

```bash
# 运行所有检查并生成报告
python scripts/check_frontend.py

# 指定提取的代码目录
python scripts/check_frontend.py -d temp_frontend_code

# 指定报告输出目录
python scripts/check_frontend.py -o lint_reports/frontend
```

#### 方法三：直接使用 npm 脚本（需要本地 Node.js）

```bash
# 检查 JavaScript
npm run lint:js

# 检查 CSS
npm run lint:css

# 检查 HTML
npm run lint:html

# 运行所有检查
npm run lint:all
```

### 3. 代码格式化

```bash
# 检查代码格式（不修改文件）
npm run format:check

# 自动格式化代码
npm run format
```

## 工具配置说明

### ESLint 配置

配置文件: `.eslintrc.json`

主要规则：
- 使用 ES2021 语法
- 支持浏览器和 Node.js 环境
- 禁止使用 `var`，推荐 `const`/`let`
- 要求使用严格相等 (`===`)
- 禁止使用 `eval()` 等危险函数
- 与 Prettier 集成

### Stylelint 配置

配置文件: `.stylelintrc.json`

主要规则：
- 遵循标准 CSS 规范
- 颜色值使用小写和短格式
- 禁止重复属性
- 与 Prettier 集成

### Prettier 配置

配置文件: `.prettierrc.json`

格式化规则：
- 使用双引号
- 行宽 100 字符
- 使用 2 空格缩进
- 行尾使用 LF
- 末尾添加分号

忽略文件: `.prettierignore`

### HTMLHint 配置

配置文件: `.htmlhintrc.json`

主要规则：
- 标签名小写
- 属性值使用双引号
- 要求标签配对
- ID 唯一性检查
- 图片需要 alt 属性

## 检查报告

检查完成后，报告将保存在 `lint_reports/frontend/` 目录：

- `frontend_lint_report_YYYYMMDD_HHMMSS.md` - Markdown 格式报告
- `frontend_lint_report_YYYYMMDD_HHMMSS.json` - JSON 格式报告

报告包含：
- 检查摘要（总计、错误、警告、通过）
- 按工具分组的问题列表
- 问题位置（文件名、行号、列号）
- 问题描述和修复建议

## 工作流程示例

### 完整检查流程

```bash
# 1. 提取前端代码
python scripts/extract_frontend_code.py FLASH_DOCK-main -o temp_frontend_code

# 2. 运行检查
python scripts/check_frontend.py -d temp_frontend_code

# 3. 查看报告
# 报告保存在 lint_reports/frontend/ 目录
```

### 开发时使用

```bash
# 格式化代码
npm run format

# 检查特定文件类型
npm run lint:js
npm run lint:css
```

## 常见问题

### Q: Node.js 未安装怎么办？

A: 有两种选择：

**选项 1：使用 Docker（推荐）**
```bash
# 安装 Docker Desktop
# 然后运行：
python scripts/check_frontend_docker.py
```

**选项 2：安装 Node.js**
- Windows: 从 https://nodejs.org/ 下载安装包
- macOS: `brew install node`
- Linux: `sudo apt install nodejs npm`

### Q: npm install 失败怎么办？

A: 尝试以下方法：
- 使用国内镜像: `npm install --registry=https://registry.npmmirror.com`
- 清除缓存: `npm cache clean --force`
- 删除 `node_modules` 和 `package-lock.json` 后重新安装

### Q: 如何只检查特定文件？

A: 可以手动指定文件：
```bash
npx eslint path/to/file.js
npx stylelint path/to/file.css
npx htmlhint path/to/file.html
```

### Q: 如何忽略某些检查规则？

A: 在代码中添加注释：
```javascript
// eslint-disable-next-line no-console
console.log('debug info');

/* stylelint-disable-next-line color-named */
color: red;

<!-- htmlhint attr-value-double-quotes:false -->
<div class='example'>
```

### Q: 如何集成到 CI/CD？

A: 在 CI 配置中添加：
```yaml
# GitHub Actions 示例
- name: Install dependencies
  run: npm install

- name: Extract frontend code
  run: python scripts/extract_frontend_code.py FLASH_DOCK-main

- name: Run linting
  run: python scripts/check_frontend.py
```

## 最佳实践

1. **定期检查**: 在提交代码前运行检查
2. **修复问题**: 优先修复错误，然后处理警告
3. **格式化代码**: 使用 Prettier 统一代码风格
4. **查看报告**: 定期查看检查报告，了解代码质量趋势
5. **更新配置**: 根据项目需求调整工具配置

## 相关资源

- [ESLint 官方文档](https://eslint.org/)
- [Stylelint 官方文档](https://stylelint.io/)
- [Prettier 官方文档](https://prettier.io/)
- [HTMLHint 官方文档](https://htmlhint.com/)

## 更新日志

- 2024-01-XX: 初始配置，添加 ESLint、Stylelint、Prettier、HTMLHint

