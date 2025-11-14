#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
环境诊断脚本：检查可能导致崩溃的环境问题
"""
import sys
import platform
import traceback
import io

# 设置输出编码为UTF-8（如果可能）
if sys.platform == 'win32':
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except:
        pass

print("=" * 60)
print("环境诊断报告")
print("=" * 60)

# 1. Python 版本
print("\n1. Python 环境:")
print(f"   Python 版本: {sys.version}")
print(f"   Python 可执行文件: {sys.executable}")
print(f"   平台: {platform.system()} {platform.release()}")
print(f"   架构: {platform.machine()}")

# 2. 关键库版本
print("\n2. 关键库版本:")
try:
    import torch
    print(f"   PyTorch: {torch.__version__}")
    print(f"   CUDA 可用: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"   CUDA 版本: {torch.version.cuda}")
except Exception as e:
    print(f"   PyTorch: 导入失败 - {e}")

try:
    import lmdb
    print(f"   LMDB: {lmdb.__version__ if hasattr(lmdb, '__version__') else '版本未知'}")
except Exception as e:
    print(f"   LMDB: 导入失败 - {e}")

try:
    import numpy as np
    print(f"   NumPy: {np.__version__}")
except Exception as e:
    print(f"   NumPy: 导入失败 - {e}")

try:
    import pickle
    print(f"   Pickle: 协议版本 {pickle.HIGHEST_PROTOCOL}")
except Exception as e:
    print(f"   Pickle: 导入失败 - {e}")

# 3. LMDB 测试
print("\n3. LMDB 功能测试:")
try:
    import lmdb
    import tempfile
    import os
    
    # 创建临时LMDB数据库测试
    with tempfile.TemporaryDirectory() as tmpdir:
        test_db = os.path.join(tmpdir, "test.db")
        env = lmdb.open(
            test_db,
            subdir=False,
            readonly=False,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=1 if platform.system() == 'Windows' else 256,
            map_size=1024*1024,  # 1MB
        )
        
        # 写入测试
        with env.begin(write=True) as txn:
            txn.put(b"0", b"test_data")
        
        # 读取测试
        with env.begin() as txn:
            data = txn.get(b"0")
            if data == b"test_data":
                print("   [OK] LMDB 基本功能正常")
            else:
                print("   [FAIL] LMDB 读取数据不匹配")
        
        # 确保环境完全关闭
        env.close()
        import time
        time.sleep(0.1)  # 给系统时间释放文件锁
        
        print("   [OK] LMDB 环境关闭正常")
        
except Exception as e:
    print(f"   [FAIL] LMDB 测试失败: {e}")
    print(f"   错误类型: {type(e).__name__}")
    if "PermissionError" in str(type(e).__name__) or "WinError 32" in str(e):
        print("   [WARNING] LMDB 文件被锁定 - 这可能是崩溃的原因！")
        print("   建议：关闭所有可能使用LMDB文件的程序，然后重试")
    traceback.print_exc()

# 4. 线程测试
print("\n4. 线程安全测试:")
try:
    import threading
    import time
    
    def test_thread():
        time.sleep(0.1)
        return True
    
    threads = []
    for i in range(5):
        t = threading.Thread(target=test_thread)
        threads.append(t)
        t.start()
    
    for t in threads:
        t.join()
    
    print("   [OK] 基本线程功能正常")
except Exception as e:
    print(f"   [FAIL] 线程测试失败: {e}")

# 5. Pickle 测试
print("\n5. Pickle 序列化测试:")
try:
    import pickle
    import numpy as np
    
    # 测试序列化numpy数组
    test_data = {
        'atoms': np.array([1, 2, 3]),
        'coordinates': np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        'smi': 'CCO',
    }
    
    pickled = pickle.dumps(test_data, protocol=pickle.HIGHEST_PROTOCOL)
    unpickled = pickle.loads(pickled)
    
    if np.array_equal(test_data['atoms'], unpickled['atoms']):
        print("   [OK] Pickle 序列化/反序列化正常")
    else:
        print("   [FAIL] Pickle 数据不匹配")
        
except Exception as e:
    print(f"   [FAIL] Pickle 测试失败: {e}")
    traceback.print_exc()

# 6. 内存测试
print("\n6. 内存访问测试:")
try:
    import numpy as np
    
    # 测试大数组分配
    large_array = np.zeros((1000, 1000), dtype=np.float32)
    large_array[0, 0] = 1.0
    del large_array
    
    print("   [OK] 内存分配/释放正常")
except Exception as e:
    print(f"   [FAIL] 内存测试失败: {e}")

# 7. 检查 unicore 库
print("\n7. Unicore 库检查:")
try:
    import unicore
    print(f"   [OK] Unicore 导入成功")
    
    try:
        from unicore import tasks, options
        print(f"   [OK] Unicore 核心模块导入成功")
    except Exception as e:
        print(f"   [FAIL] Unicore 核心模块导入失败: {e}")
        
except Exception as e:
    print(f"   [FAIL] Unicore 导入失败: {e}")

# 8. 检查 unimol 模块
print("\n8. Unimol 模块检查:")
try:
    import sys
    import os
    
    unimol_path = os.path.join(os.path.dirname(__file__), '..', 'unimol')
    if os.path.exists(unimol_path):
        sys.path.insert(0, unimol_path)
        from unimol.data import LMDBDataset
        print(f"   [OK] Unimol 模块导入成功")
    else:
        print(f"   [FAIL] Unimol 路径不存在: {unimol_path}")
        
except Exception as e:
    print(f"   [FAIL] Unimol 模块导入失败: {e}")
    traceback.print_exc()

# 9. Windows 特定检查
if platform.system() == 'Windows':
    print("\n9. Windows 特定检查:")
    try:
        import ctypes
        
        # 检查是否可以使用 Windows API
        kernel32 = ctypes.windll.kernel32
        print(f"   [OK] Windows API 可访问")
        
        # 检查路径长度限制
        max_path = 260  # Windows 默认最大路径长度
        print(f"   最大路径长度: {max_path} 字符")
        
        # 检查当前工作目录路径长度
        cwd = os.getcwd()
        if len(cwd) > max_path:
            print(f"   [WARNING] 当前工作目录路径过长 ({len(cwd)} 字符)")
        else:
            print(f"   [OK] 当前工作目录路径长度正常 ({len(cwd)} 字符)")
            
    except Exception as e:
        print(f"   [FAIL] Windows 特定检查失败: {e}")

# 10. 建议
print("\n" + "=" * 60)
print("诊断建议:")
print("=" * 60)

issues = []

if platform.system() == 'Windows':
    print("\nWindows 特定建议:")
    print("1. 确保使用 Python 3.8-3.11 (Python 3.12 可能有兼容性问题)")
    print("2. 尝试使用 conda 安装 LMDB: conda install -c conda-forge python-lmdb")
    print("3. 检查是否有防病毒软件干扰")
    print("4. 尝试以管理员权限运行")
    print("5. 检查 Windows 事件查看器中的应用程序错误")

print("\n通用建议:")
print("1. 尝试降低 batch_size (例如 --batch-size 1)")
print("2. 检查 LMDB 文件是否损坏")
print("3. 尝试重新创建 LMDB 文件")
print("4. 检查是否有其他进程正在访问 LMDB 文件")

print("\n" + "=" * 60)
print("诊断完成")
print("=" * 60)

