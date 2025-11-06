# 快速启动指南

## 一键启动（最简单）

双击运行 `start_all_services.bat`

等待三个窗口打开：
1. 服务注册中心
2. COMPASS服务  
3. FLASH-DOCK

## 手动启动步骤

### 步骤1：启动注册中心
```cmd
conda activate AIDDTRAIN
cd E:\Qinchaojun\AIDD-TRAIN
python services/registry/server.py
```

### 步骤2：启动COMPASS服务
```cmd
conda activate AIDDTRAIN
cd E:\Qinchaojun\AIDD-TRAIN
python compass/service_main.py
```

### 步骤3：启动FLASH-DOCK
```cmd
conda activate flash_dock
cd E:\Qinchaojun\AIDD-TRAIN\FLASH_DOCK-main
streamlit run FlashDock.py
```

## 访问地址

- FLASH-DOCK: http://localhost:8501
- COMPASS API文档: http://localhost:8080/docs
- 注册中心: http://localhost:8500

## 注意事项

1. **必须按顺序启动**：先注册中心，再COMPASS，最后FLASH-DOCK
2. **每个服务单独窗口**：不要关闭窗口
3. **等待启动完成**：看到"started"或"running"消息后再启动下一个

