# Java安装指南 - WSL版本

## 📋 重要说明

**根据项目配置，FlashDock应用在WSL（Ubuntu-24.04）中运行，因此Java必须在WSL中安装！**

## 🔍 如何确认运行环境

### 检查Streamlit运行位置

如果您的FlashDock是通过以下方式启动的：
- `start_flashdock_wsl.bat`
- `start_flashdock_fixed.bat`
- `fix_and_start_services.py`

那么Streamlit**在WSL中运行**，Java也必须在WSL中安装。

### 检查方法

在WSL中运行：
```bash
wsl -d Ubuntu-24.04
ps aux | grep streamlit
```

如果看到streamlit进程，说明在WSL中运行。

---

## 🚀 在WSL中安装Java

### 方法一：使用apt安装（推荐，最简单）

#### 步骤1：进入WSL
```bash
wsl -d Ubuntu-24.04
```

#### 步骤2：更新包列表
```bash
sudo apt update
```

#### 步骤3：安装OpenJDK 17（推荐）
```bash
sudo apt install -y openjdk-17-jdk
```

或者安装OpenJDK 21：
```bash
sudo apt install -y openjdk-21-jdk
```

#### 步骤4：验证安装
```bash
java -version
```

应该显示类似：
```
openjdk version "17.0.x" 2024-xx-xx
OpenJDK Runtime Environment (build 17.0.x+xx-Ubuntu-...)
OpenJDK 64-Bit Server VM (build 17.0.x+xx-Ubuntu-..., mixed mode, sharing)
```

#### 步骤5：设置JAVA_HOME（可选但推荐）
```bash
# 查找Java安装路径
sudo update-alternatives --config java
# 或
readlink -f $(which java)
```

通常Java安装在：`/usr/lib/jvm/java-17-openjdk-amd64`

添加到环境变量：
```bash
# 编辑 ~/.bashrc
nano ~/.bashrc

# 添加以下行（根据实际路径调整）
export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH

# 保存后重新加载
source ~/.bashrc
```

---

### 方法二：使用Adoptium（手动安装）

#### 步骤1：进入WSL
```bash
wsl -d Ubuntu-24.04
```

#### 步骤2：下载Adoptium JDK
```bash
# 创建临时目录
cd /tmp

# 下载Java 17（根据系统架构选择）
# 对于x64系统：
wget https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.13%2B11/OpenJDK17U-jdk_x64_linux_hotspot_17.0.13_11.tar.gz

# 或者下载Java 21：
# wget https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.7%2B11/OpenJDK21U-jdk_x64_linux_hotspot_21.0.7_11.tar.gz
```

#### 步骤3：解压并安装
```bash
# 解压
tar -xzf OpenJDK17U-jdk_x64_linux_hotspot_*.tar.gz

# 移动到系统目录
sudo mv jdk-17.0.13+11 /opt/java-17

# 设置环境变量
echo 'export JAVA_HOME=/opt/java-17' >> ~/.bashrc
echo 'export PATH=$JAVA_HOME/bin:$PATH' >> ~/.bashrc

# 重新加载
source ~/.bashrc
```

#### 步骤4：验证安装
```bash
java -version
```

---

## 🔍 验证安装

### 在WSL中运行诊断脚本

```bash
# 进入WSL
wsl -d Ubuntu-24.04

# 进入项目目录
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main

# 运行诊断脚本（需要在WSL中运行）
python3 check_java_environment.py
```

或者使用我创建的WSL专用诊断脚本：
```bash
bash check_java_environment_wsl.sh
```

---

## ⚠️ 常见问题

### 问题1：在Windows中安装了Java，但WSL中找不到

**原因**：Windows和WSL是独立的系统，环境变量不共享。

**解决**：必须在WSL中单独安装Java。

### 问题2：安装后仍然找不到Java

**检查步骤**：
1. 确认在WSL中运行 `java -version`
2. 检查PATH环境变量：`echo $PATH`
3. 检查JAVA_HOME：`echo $JAVA_HOME`
4. 重新加载环境变量：`source ~/.bashrc`

### 问题3：Java版本不符合要求

**要求**：Java 17-23

**检查版本**：
```bash
java -version
```

**升级/降级**：
```bash
# 卸载旧版本
sudo apt remove openjdk-*-jdk

# 安装正确版本
sudo apt install -y openjdk-17-jdk
```

---

## 📝 完整安装示例

```bash
# 1. 进入WSL
wsl -d Ubuntu-24.04

# 2. 更新系统
sudo apt update

# 3. 安装Java 17
sudo apt install -y openjdk-17-jdk

# 4. 验证安装
java -version

# 5. 设置JAVA_HOME（可选）
echo 'export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64' >> ~/.bashrc
echo 'export PATH=$JAVA_HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

# 6. 进入项目目录
cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main

# 7. 运行诊断脚本
python3 check_java_environment.py

# 8. 如果所有检查通过，重启Streamlit应用
```

---

## ✅ 安装完成后的检查清单

- [ ] Java已安装在WSL中（`wsl -d Ubuntu-24.04` 然后 `java -version`）
- [ ] Java版本在17-23之间
- [ ] JAVA_HOME环境变量已设置（可选，但推荐）
- [ ] 诊断脚本在WSL中运行显示所有检查通过
- [ ] 已重启Streamlit应用（在WSL中运行的）
- [ ] "加载示例文件"功能可以正常使用

---

## 🔄 重启Streamlit应用

安装Java后，需要重启在WSL中运行的Streamlit应用：

1. **停止当前应用**：
   - 如果通过批处理脚本启动，关闭对应的窗口
   - 或在WSL中运行：`pkill -f streamlit`

2. **重新启动**：
   ```bash
   # 在Windows中
   start_flashdock_wsl.bat
   
   # 或在WSL中手动启动
   wsl -d Ubuntu-24.04
   source ~/miniconda3/etc/profile.d/conda.sh
   conda activate flash_dock
   export PYTHONPATH=/mnt/e/Qinchaojun/AIDD-TRAIN
   cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main
   streamlit run FlashDock.py --server.port 8501 --server.address 0.0.0.0
   ```

---

## 📞 需要帮助？

如果遇到问题，请提供：
1. WSL中运行 `java -version` 的输出
2. WSL中运行 `echo $JAVA_HOME` 的输出
3. WSL中运行诊断脚本的完整输出
4. 操作系统版本和WSL版本

