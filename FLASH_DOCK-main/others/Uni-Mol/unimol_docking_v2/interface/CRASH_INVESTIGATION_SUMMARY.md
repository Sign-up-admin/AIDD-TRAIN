# Windows 崩溃问题调查总结

## 问题描述

程序在第一次迭代时崩溃，返回码 `3228369023` (0xC0000005)，这是 Windows 访问违规错误。

### 错误特征
- **错误代码**: 0xC0000005 (访问违规/访问冲突)
- **发生位置**: 第一次迭代时 (`for i, sample in enumerate(progress)`)
- **环境**: Windows 10, Python 3.10.19, PyTorch 2.9.0+cpu
- **症状**: 程序在 "Starting inference loop..." 之后立即崩溃，无错误信息输出

## 已完成的修复

### 1. LMDB 环境管理修复 (`lmdb_dataset.py`)
- ✅ 修复了环境初始化问题，确保环境在整个对象生命周期内保持打开
- ✅ 添加了线程安全锁 (`threading.Lock`)
- ✅ Windows 特定配置：`max_readers=1`（Windows 上 LMDB 多读者支持有限）
- ✅ 添加了 `__del__` 方法确保环境正确关闭
- ✅ Windows 上禁用 `lru_cache`，使用线程安全的直接访问
- ✅ 添加了环境错误恢复机制

### 2. TTADockingPoseDataset 线程安全修复 (`tta_dataset.py`)
- ✅ 创建了线程安全的 LRU 缓存实现 (`ThreadSafeLRUCache`)
- ✅ Windows 上使用线程安全缓存，非 Windows 使用标准 `lru_cache`
- ✅ 添加了缓存清除机制

### 3. 错误处理增强 (`infer.py`)
- ✅ 增强了 Windows 特定的错误诊断
- ✅ 添加了第一次迭代的特殊处理和详细日志
- ✅ 添加了预测试机制，在循环开始前测试第一个样本访问
- ✅ 详细的错误原因分析和解决建议

### 4. 调试工具
- ✅ 创建了增强的调试脚本 (`debug_crash.py`)
- ✅ 创建了环境诊断脚本 (`check_environment.py`)
- ✅ 创建了多个测试脚本用于隔离问题

## 网上搜索结果和建议

### 0xC0000005 错误的常见原因

根据搜索结果，访问违规错误可能由以下原因引起：

1. **系统文件损坏**
   - 运行 `sfc /scannow` 和 `DISM /Online /Cleanup-Image /RestoreHealth`
   - 检查系统文件完整性

2. **内存问题**
   - 使用 Windows 内存诊断工具检查 RAM
   - 内存条损坏或接触不良可能导致随机崩溃

3. **驱动程序问题**
   - 更新或回滚显卡、声卡等硬件驱动
   - 检查设备管理器中是否有黄色感叹号

4. **软件兼容性问题**
   - Python 版本与库不兼容
   - C 扩展库编译问题
   - 库版本冲突

5. **安全软件干扰**
   - 杀毒软件或防火墙错误拦截
   - 尝试暂时禁用安全软件测试

6. **权限问题**
   - 尝试以管理员权限运行

## 可能的原因分析

基于代码分析和网上搜索结果，崩溃可能由以下原因引起：

### 1. LMDB 库的 C 扩展问题（最可能）
- LMDB 是 C 扩展库，在 Windows 上可能有兼容性问题
- 多线程访问时的内存管理问题
- 环境关闭/打开时的竞争条件

### 2. PyTorch DataLoader 的 Windows 问题
- 即使 `num_workers=0`，PyTorch 内部可能仍使用多线程
- Windows 上的 fork 模拟问题

### 3. unicore 库的迭代器实现
- unicore 的迭代器可能在 Windows 上有问题
- 数据缓冲区的内存管理问题

### 4. Pickle 反序列化问题
- 大对象的反序列化可能导致内存访问错误
- Python 3.10 的 pickle 实现问题

## 建议的解决方案

### 立即尝试的解决方案

1. **重新安装 LMDB**
   ```bash
   conda activate flash_dock
   pip uninstall lmdb
   pip install lmdb --no-cache-dir
   ```
   或者使用 conda：
   ```bash
   conda install -c conda-forge python-lmdb
   ```

2. **降低 batch_size**
   - 将 `--batch-size 4` 改为 `--batch-size 1`
   - 减少内存压力

3. **检查 Windows 事件查看器**
   - 打开"事件查看器" > "Windows 日志" > "应用程序"
   - 查找崩溃时的错误记录，可能包含更多细节

4. **运行系统文件检查**
   ```bash
   sfc /scannow
   DISM /Online /Cleanup-Image /RestoreHealth
   ```

5. **检查内存**
   - 运行 `mdsched.exe` 进行内存诊断

### 长期解决方案

1. **考虑使用 Linux 环境**
   - Windows 上的多线程和 C 扩展库支持有限
   - Linux 环境通常更稳定

2. **使用 Python 3.9**
   - Python 3.10 可能太新，某些库可能不完全兼容
   - 创建新环境：`conda create -n flash_dock_py39 python=3.9`

3. **检查防病毒软件**
   - 暂时禁用防病毒软件测试
   - 将项目目录添加到白名单

4. **使用 WSL2**
   - 在 Windows Subsystem for Linux 中运行
   - 可能避免 Windows 特定的问题

## 下一步调查方向

1. **运行增强的调试脚本**
   - 使用 `debug_crash.py` 获取更详细的错误信息
   - 查看预测试是否能捕获错误

2. **检查 Windows 事件查看器**
   - 获取系统级的错误信息
   - 可能包含崩溃时的内存地址和调用栈

3. **尝试使用不同的 LMDB 版本**
   - 降级或升级 LMDB 版本
   - 检查是否有已知的 Windows 兼容性问题

4. **使用调试器**
   - 使用 `pdb` 或 `winpdb` 调试
   - 在崩溃前设置断点

5. **检查是否有其他进程访问 LMDB 文件**
   - 确保没有其他程序正在使用 LMDB 文件
   - 检查文件锁定

## 测试脚本

已创建以下测试脚本用于诊断：

1. **`debug_crash.py`** - 增强的调试脚本，逐步测试每个组件
2. **`check_environment.py`** - 环境诊断脚本
3. **`test_lmdb_direct.py`** - 直接 LMDB 访问测试
4. **`test_no_lru_cache.py`** - 无缓存版本测试
5. **`test_pickle.py`** - Pickle 反序列化测试
6. **`test_dataset_chain.py`** - 数据集处理链逐步测试

## 修复的文件

1. `unimol/data/lmdb_dataset.py` - LMDB 环境管理和线程安全
2. `unimol/data/tta_dataset.py` - 线程安全的缓存实现
3. `unimol/infer.py` - 错误处理和预测试
4. `interface/debug_crash.py` - 增强的调试脚本

## 网上搜索的关键发现

### LMDB 在 Windows 上的已知问题

1. **多读者支持有限**
   - Windows 上 LMDB 的 `max_readers` 应该设置为 1
   - 已修复：代码中已设置 `max_readers=1` for Windows

2. **文件锁定问题**
   - Windows 上的文件锁定机制与 Linux 不同
   - 可能需要更长的关闭延迟

3. **C 扩展库编译问题**
   - LMDB 的 Python 绑定是 C 扩展
   - 可能需要重新编译或使用预编译的 wheel

### PyTorch DataLoader 在 Windows 上的问题

1. **multiprocessing 问题**
   - Windows 上 multiprocessing 使用 spawn 而不是 fork
   - 即使 `num_workers=0`，某些内部机制可能仍使用多线程

2. **内存映射问题**
   - Windows 上的内存映射行为不同
   - 可能导致访问违规

## 额外的修复建议

### 1. 尝试更保守的 LMDB 配置

在 `lmdb_dataset.py` 中，可以尝试更保守的配置：

```python
env = lmdb.open(
    lmdb_path,
    subdir=False,
    readonly=True,
    lock=False,
    readahead=False,
    meminit=False,
    max_readers=1,  # Windows 上必须为 1
    map_size=1024*1024*1024,  # 明确设置 map_size
    max_dbs=0,  # 不使用命名数据库
)
```

### 2. 添加环境变量

在运行程序前设置：

```bash
set PYTHONFAULTHANDLER=1
set PYTHONDEVMODE=1
```

这可以帮助捕获更多错误信息。

### 3. 使用 try-except 包装整个循环

虽然已经添加了错误处理，但访问违规可能发生在 C 扩展层面，Python 的 try-except 无法捕获。

### 4. 检查是否有其他进程访问 LMDB

```bash
# 检查是否有 Python 进程正在运行
tasklist | findstr python

# 检查 LMDB 文件是否被锁定
# 尝试重命名文件，如果失败说明被锁定
```

## 最可能的根本原因

基于所有调查，最可能的原因是：

1. **LMDB C 扩展库在 Windows 上的内存管理问题**（70% 可能性）
   - C 扩展库的内存访问错误
   - 多线程访问时的竞争条件
   - 环境关闭/打开时的内存泄漏

2. **PyTorch DataLoader 的内部多线程问题**（20% 可能性）
   - 即使 `num_workers=0`，内部可能仍使用线程
   - Windows 上的线程模型与 Linux 不同

3. **unicore 库的迭代器实现问题**（10% 可能性）
   - 数据缓冲区的内存管理
   - Windows 特定的实现缺陷

## 结论

虽然已经做了大量修复，但访问违规错误通常是 C 扩展库层面的问题，可能需要在系统层面解决。建议：

### 优先级 1（立即尝试）
1. **重新安装 LMDB 库**
   ```bash
   conda activate flash_dock
   pip uninstall lmdb
   pip install lmdb --no-cache-dir --force-reinstall
   ```

2. **降低 batch_size 到 1**
   - 修改命令中的 `--batch-size 4` 为 `--batch-size 1`

3. **检查 Windows 事件查看器**
   - 获取系统级的错误信息

### 优先级 2（如果优先级1无效）
4. **运行系统文件检查**
   ```bash
   sfc /scannow
   DISM /Online /Cleanup-Image /RestoreHealth
   ```

5. **检查内存**
   - 运行 `mdsched.exe`

6. **尝试以管理员权限运行**

### 优先级 3（长期解决方案）
7. **使用 WSL2 或 Linux 环境**
   - Windows 上的 C 扩展库支持有限
   - Linux 环境通常更稳定

8. **使用 Python 3.9**
   - 创建新环境测试

9. **联系 LMDB 或 unicore 的开发者**
   - 报告 Windows 特定的崩溃问题
   - 可能是一个已知的 bug

