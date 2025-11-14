# 服务环境配置说明

## 环境分离策略

本项目使用**环境分离**策略，不同服务使用不同的Conda环境，避免依赖冲突。

## 环境配置

### 1. AIDDTRAIN 环境
- **用途**: COMPASS服务和Registry服务
- **环境路径**: `D:\conda_envs\AIDDTRAIN` (推荐) 或 `C:\ProgramData\Anaconda3\envs\AIDDTRAIN`
- **Python版本**: 3.10.19
- **关键依赖**:
  - PyTorch 2.9.0+cpu
  - torch-geometric 2.7.0
  - rdkit-pypi
  - FastAPI, Uvicorn
  - 其他COMPASS服务依赖

### 2. flash_dock 环境
- **用途**: FLASH-DOCK前端服务
- **环境路径**: `C:\ProgramData\Anaconda3\envs\flash_dock`
- **Python版本**: 3.8.20
- **关键依赖**:
  - streamlit
  - streamlit-molstar
  - rdkit-pypi==2022.9.3
  - numpy<2.0.0 (RDKit要求)
  - py3dmol, stmol
  - 其他FLASH-DOCK依赖

## 服务启动配置

### 自动启动脚本 (`start_all_services.bat`)

脚本会自动检测并使用正确的环境：

1. **Registry服务** (端口8500)
   - 使用: AIDDTRAIN环境
   - Python: `D:\conda_envs\AIDDTRAIN\python.exe`

2. **COMPASS服务** (端口8080)
   - 使用: AIDDTRAIN环境
   - Python: `D:\conda_envs\AIDDTRAIN\python.exe`
   - 需要: rdkit模块

3. **FLASH-DOCK服务** (端口8501)
   - 使用: **flash_dock环境** (独立环境)
   - Python: `C:\ProgramData\Anaconda3\envs\flash_dock\python.exe`
   - 需要: streamlit-molstar模块

## 环境检测顺序

### AIDDTRAIN环境检测顺序:
1. `D:\conda_envs\AIDDTRAIN\python.exe` (推荐)
2. `C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe`
3. `%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe`

### flash_dock环境检测顺序:
1. `C:\ProgramData\Anaconda3\envs\flash_dock\python.exe`
2. `D:\conda_envs\flash_dock\python.exe`
3. `%USERPROFILE%\anaconda3\envs\flash_dock\python.exe`

## 验证环境

### 检查AIDDTRAIN环境:
```bash
D:\conda_envs\AIDDTRAIN\python.exe -c "import rdkit; print('rdkit OK')"
```

### 检查flash_dock环境:
```bash
C:\ProgramData\Anaconda3\envs\flash_dock\python.exe -c "import streamlit_molstar; print('streamlit_molstar OK')"
```

## 常见问题

### 问题1: FLASH-DOCK找不到streamlit_molstar
**原因**: 使用了错误的Python环境（使用了AIDDTRAIN而不是flash_dock）

**解决**: 确保启动脚本使用flash_dock环境
- 检查 `start_all_services.bat` 是否正确检测到flash_dock环境
- 手动指定: `C:\ProgramData\Anaconda3\envs\flash_dock\python.exe -m streamlit run FlashDock.py`

### 问题2: COMPASS服务找不到rdkit
**原因**: 使用了错误的Python环境或rdkit未安装

**解决**: 
- 确保使用AIDDTRAIN环境
- 安装rdkit: `pip install rdkit-pypi`

### 问题3: 环境路径不匹配
**解决**: 
- 检查环境实际位置: `conda env list`
- 更新 `start_all_services.bat` 中的路径检测逻辑

## 手动启动服务

### 启动Registry服务:
```bash
D:\conda_envs\AIDDTRAIN\python.exe services\registry\server.py --host 0.0.0.0 --port 8500
```

### 启动COMPASS服务:
```bash
D:\conda_envs\AIDDTRAIN\python.exe compass\service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500
```

### 启动FLASH-DOCK服务:
```bash
C:\ProgramData\Anaconda3\envs\flash_dock\python.exe -m streamlit run FLASH_DOCK-main\FlashDock.py --server.port 8501
```

## 注意事项

1. **环境隔离**: FLASH-DOCK必须使用自己的环境，不能与COMPASS共享
2. **依赖版本**: flash_dock环境需要numpy<2.0.0，这是RDKit的要求
3. **路径引用**: 启动脚本使用绝对路径引用Python可执行文件，确保路径正确
4. **PYTHONPATH**: 所有服务都需要设置PYTHONPATH指向项目根目录



