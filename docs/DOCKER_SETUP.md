# Docker 设置指南

## Docker 镜像加速配置

如果遇到 Docker Hub 连接问题，需要配置国内镜像加速器。

### Windows Docker Desktop

1. 打开 Docker Desktop
2. 点击设置（Settings）
3. 选择 "Docker Engine"
4. 在 JSON 配置中添加镜像加速器：

```json
{
  "registry-mirrors": [
    "https://docker.mirrors.ustc.edu.cn",
    "https://hub-mirror.c.163.com",
    "https://mirror.baidubce.com"
  ]
}
```

5. 点击 "Apply & Restart"

### 使用国内镜像源构建

如果配置了镜像加速器后仍然无法访问，可以使用国内镜像源：

```bash
# 使用腾讯云镜像
docker build -f Dockerfile.frontend-lint.cn -t aidd-frontend-linter:latest .
```

### 手动拉取镜像

也可以先手动拉取镜像：

```bash
# 使用镜像加速器拉取
docker pull node:18-alpine

# 或使用国内镜像
docker pull ccr.ccs.tencentyun.com/library/node:18-alpine
docker tag ccr.ccs.tencentyun.com/library/node:18-alpine node:18-alpine
```

## 验证配置

```bash
# 检查 Docker 版本
docker --version

# 测试拉取镜像
docker pull hello-world

# 检查镜像加速器配置
docker info | grep -A 10 "Registry Mirrors"
```

## 常见问题

### Q: 连接超时怎么办？

A: 
1. 检查网络连接
2. 配置镜像加速器（见上方）
3. 使用代理（如果可用）

### Q: 镜像拉取失败？

A:
1. 尝试使用不同的镜像源
2. 检查防火墙设置
3. 使用 VPN 或代理

### Q: 如何验证镜像加速器是否生效？

A:
```bash
docker info | findstr "Registry Mirrors"
# 应该显示配置的镜像地址
```



