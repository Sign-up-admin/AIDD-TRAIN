# WSL2 网络问题分析报告

## 检查结果

### 1. 基本网络连接 ✅
- **Ping 测试**: 可以 ping 通 8.8.8.8（延迟 ~190ms）
- **HTTPS 连接**: curl 可以连接到 files.pythonhosted.org（返回 HTTP/2 200）
- **IP 地址**: 172.30.226.241/20（正常 WSL2 网络配置）
- **DNS 服务器**: 10.255.255.254（WSL 虚拟 DNS）

### 2. 问题诊断

#### 问题现象
- pip 下载包时频繁超时
- 连接速度很慢（延迟 190ms+）
- 下载大文件时容易中断

#### 可能原因

1. **网络速度慢**
   - 延迟较高（190ms）
   - 可能是网络环境或防火墙限制

2. **pip 默认超时时间过短**
   - pip 默认超时 15 秒
   - 对于慢速网络可能不够

3. **DNS 解析可能较慢**
   - 虽然能解析，但可能较慢

4. **Windows 防火墙或代理设置**
   - 可能影响 WSL2 的网络性能

## 解决方案

### 方案 1: 增加 pip 超时和重试（已实施）
```bash
pip install --timeout=100 --retries=10 <package>
```

### 方案 2: 配置 DNS（如果 DNS 慢）
```bash
# 在 WSL2 中
sudo bash -c 'cat > /etc/wsl.conf << EOF
[network]
generateResolvConf = false
EOF'

sudo bash -c 'cat > /etc/resolv.conf << EOF
nameserver 8.8.8.8
nameserver 8.8.4.4
EOF'
```

### 方案 3: 使用国内镜像源（推荐）
```bash
# 配置 pip 使用清华镜像
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
pip config set global.trusted-host pypi.tuna.tsinghua.edu.cn
```

### 方案 4: 配置代理（如果有代理）
```bash
# 在 WSL2 中设置代理
export http_proxy=http://127.0.0.1:10808
export https_proxy=http://127.0.0.1:10808
```

### 方案 5: 使用 conda 代替 pip（部分包）
```bash
# conda 通常更稳定
conda install -c conda-forge <package>
```

## 推荐操作

1. **立即尝试**: 使用国内镜像源（方案 3）
2. **如果仍慢**: 增加超时时间（方案 1）
3. **如果有代理**: 配置代理（方案 4）

## 当前状态

- ✅ 网络基本连通
- ⚠️ 速度较慢，需要优化
- ⚠️ pip 下载容易超时

