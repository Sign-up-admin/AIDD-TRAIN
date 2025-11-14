# Docker Node.js 快速开始

## 当前状态

Docker 已安装，但需要配置镜像加速器才能从 Docker Hub 拉取镜像。

## 解决方案

### 方案一：配置 Docker 镜像加速器（推荐）

1. **打开 Docker Desktop**
2. **点击设置（Settings）** → **Docker Engine**
3. **添加以下配置**：

```json
{
  "registry-mirrors": [
    "https://docker.mirrors.ustc.edu.cn",
    "https://hub-mirror.c.163.com",
    "https://mirror.baidubce.com"
  ]
}
```

4. **点击 "Apply & Restart"**

5. **验证配置**：
```bash
docker info | findstr "Registry Mirrors"
```

6. **重新构建镜像**：
```bash
python scripts/check_frontend_docker.py --build
```

### 方案二：使用国内镜像源构建

如果镜像加速器配置后仍有问题，可以使用国内镜像源：

```bash
docker build -f Dockerfile.frontend-lint.cn -t aidd-frontend-linter:latest .
```

### 方案三：手动拉取镜像

```bash
# 使用镜像加速器拉取
docker pull node:18-alpine

# 然后构建
docker build -f Dockerfile.frontend-lint -t aidd-frontend-linter:latest .
```

## 配置完成后使用

```bash
# 运行前端代码检查
python scripts/check_frontend_docker.py

# 或指定命令
python scripts/check_frontend_docker.py lint:js
python scripts/check_frontend_docker.py lint:css
python scripts/check_frontend_docker.py format
```

## 详细文档

- Docker 设置指南: `docs/DOCKER_SETUP.md`
- Docker vs 本地 Node.js 对比: `docs/DOCKER_NODEJS_COMPARISON.md`
- 前端代码质量指南: `docs/FRONTEND_CODE_QUALITY.md`



