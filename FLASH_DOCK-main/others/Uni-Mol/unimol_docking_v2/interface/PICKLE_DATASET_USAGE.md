# PickleDataset 使用说明

## 为什么需要 PickleDataset？

LMDB 在 Windows 上存在兼容性问题，可能导致访问违规错误（0xC0000005）。`PickleDataset` 使用标准的 pickle 文件替代 LMDB，完全避免 LMDB 的 C 扩展库问题。

## 工作原理

`PickleDataset` 将 LMDB 数据转换为 pickle 文件格式：
- **LMDB 格式**: `data.lmdb` (单个文件)
- **Pickle 格式**: `data.pickle_dir/` (目录)
  - `0.pkl`
  - `1.pkl`
  - `2.pkl`
  - ...
  - `_metadata.pkl` (元数据)

## 自动转换

`PickleDataset` 会自动检测：
1. 如果 pickle 目录存在，直接使用
2. 如果只有 LMDB 文件，自动转换为 pickle 格式
3. 转换只需一次，之后直接使用 pickle 格式

## 使用方法

### 方法 1: 自动使用（推荐）

在 Windows 上，程序会自动使用 `PickleDataset`，无需任何配置。

### 方法 2: 手动启用

设置环境变量强制使用 pickle 格式：

```bash
# Windows PowerShell
$env:UNIMOL_USE_PICKLE="1"

# Windows CMD
set UNIMOL_USE_PICKLE=1

# Linux/Mac
export UNIMOL_USE_PICKLE=1
```

### 方法 3: 手动转换

如果需要提前转换，可以使用 Python 脚本：

```python
from unimol.data.pickle_dataset import convert_lmdb_to_pickle

# 转换单个文件
convert_lmdb_to_pickle('path/to/data.lmdb')

# 转换后会在同一目录生成 data.pickle_dir/
```

## 性能对比

| 特性 | LMDB | PickleDataset |
|------|------|---------------|
| Windows 兼容性 | ❌ 有问题 | ✅ 完全兼容 |
| 读取速度 | 快 | 稍慢（但可接受） |
| 文件大小 | 小（单个文件） | 稍大（多个文件） |
| 内存使用 | 低 | 稍高 |
| 线程安全 | 需要特殊处理 | 原生支持 |

## 注意事项

1. **首次使用**: 第一次使用时会自动转换，可能需要一些时间
2. **磁盘空间**: Pickle 格式可能占用稍多磁盘空间
3. **文件数量**: 每个数据项对应一个文件，如果数据量很大，文件数量会很多
4. **删除 LMDB**: 转换后可以删除原始 LMDB 文件以节省空间

## 故障排除

### 问题：转换失败

**原因**: LMDB 库未安装或损坏

**解决**: 
```bash
pip install lmdb
# 或者使用 conda
conda install -c conda-forge python-lmdb
```

### 问题：权限错误

**原因**: 没有写入权限

**解决**: 确保对数据目录有写入权限

### 问题：磁盘空间不足

**原因**: Pickle 格式占用更多空间

**解决**: 清理不需要的文件，或使用更大的磁盘

## 代码示例

```python
from unimol.data.pickle_dataset import PickleDataset

# 自动检测和转换
dataset = PickleDataset('path/to/data.lmdb')

# 使用方式与 LMDBDataset 完全相同
print(f"Dataset size: {len(dataset)}")
item = dataset[0]
print(f"First item keys: {item.keys()}")
```

## 优势

✅ **完全避免 LMDB 的 Windows 兼容性问题**
✅ **使用标准 Python 库，无需 C 扩展**
✅ **自动转换，无需手动操作**
✅ **线程安全，支持多线程访问**
✅ **向后兼容，API 与 LMDBDataset 相同**

## 总结

如果遇到 Windows 上的 LMDB 崩溃问题，使用 `PickleDataset` 是最简单有效的解决方案。程序会自动处理转换，你只需要正常运行即可。

