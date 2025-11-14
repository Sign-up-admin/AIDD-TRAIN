#!/usr/bin/env python3
"""
测试pickle反序列化，验证数据完整性和兼容性
专门测试Windows上的pickle问题
"""
import os
import sys
import platform
import logging
import traceback
import pickle
import lmdb
import numpy as np

# 设置日志
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    stream=sys.stdout,
)
logger = logging.getLogger("test_pickle")

def test_pickle_protocol_compatibility(lmdb_path):
    """测试pickle协议版本兼容性"""
    logger.info("=" * 60)
    logger.info("Test 1: Pickle Protocol Compatibility")
    logger.info("=" * 60)
    
    env = None
    try:
        env = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=256,
        )
        
        with env.begin() as txn:
            keys = list(txn.cursor().iternext(values=False))
        
        if len(keys) == 0:
            logger.warning("No keys found, skipping test")
            return
        
        # 读取第一个数据
        first_key = keys[0]
        with env.begin() as txn:
            original_data = txn.get(first_key)
        
        if not original_data:
            logger.warning("No data found for first key")
            return
        
        # 反序列化原始数据
        logger.info("Deserializing original data...")
        try:
            original_obj = pickle.loads(original_data)
            logger.info(f"[OK] Original deserialization successful, type: {type(original_obj)}")
        except Exception as e:
            logger.error(f"[FAIL] Original deserialization failed: {e}")
            raise
        
        # 测试不同协议版本
        logger.info("Testing different pickle protocol versions...")
        supported_protocols = []
        for protocol in range(6):  # 协议版本 0-5
            try:
                # 重新序列化
                repickled = pickle.dumps(original_obj, protocol=protocol)
                # 重新反序列化
                reunpickled = pickle.loads(repickled)
                supported_protocols.append(protocol)
                logger.info(f"  Protocol {protocol}: OK")
            except Exception as e:
                logger.warning(f"  Protocol {protocol}: FAILED - {e}")
        
        logger.info(f"[OK] Supported protocols: {supported_protocols}")
        
        # 检查原始数据的协议版本
        try:
            # 尝试检测原始数据的协议版本
            if original_data[0] == 0x80:  # 二进制协议
                protocol_version = original_data[1]
                logger.info(f"Original data appears to use protocol version: {protocol_version}")
        except:
            pass
        
    except Exception as e:
        logger.error(f"[FAIL] Protocol compatibility test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise
    finally:
        if env is not None:
            try:
                env.close()
            except:
                pass

def test_pickle_data_integrity(lmdb_path):
    """测试pickle数据完整性"""
    logger.info("\n" + "=" * 60)
    logger.info("Test 2: Pickle Data Integrity")
    logger.info("=" * 60)
    
    env = None
    try:
        env = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=256,
        )
        
        with env.begin() as txn:
            keys = list(txn.cursor().iternext(values=False))
        
        logger.info(f"Testing {min(10, len(keys))} keys...")
        
        success_count = 0
        error_count = 0
        error_details = []
        
        for i, key in enumerate(keys[:10]):
            try:
                with env.begin() as txn:
                    data = txn.get(key)
                
                if not data:
                    logger.warning(f"  Key {i+1}: No data")
                    continue
                
                # 尝试反序列化
                obj = pickle.loads(data)
                success_count += 1
                
                # 检查数据类型
                if isinstance(obj, dict):
                    logger.debug(f"  Key {i+1}: OK, dict with {len(obj)} keys")
                    # 检查常见的数据类型
                    for k, v in list(obj.items())[:3]:
                        if isinstance(v, np.ndarray):
                            logger.debug(f"    {k}: numpy array, shape={v.shape}, dtype={v.dtype}")
                        elif isinstance(v, (list, tuple)):
                            logger.debug(f"    {k}: {type(v).__name__}, len={len(v)}")
                        else:
                            logger.debug(f"    {k}: {type(v).__name__}")
                else:
                    logger.debug(f"  Key {i+1}: OK, type={type(obj)}")
                
            except Exception as e:
                error_count += 1
                error_details.append((i+1, key, str(e), type(e).__name__))
                logger.error(f"  Key {i+1} ({key}): FAILED - {e}")
                logger.error(f"    Error type: {type(e).__name__}")
        
        logger.info(f"[OK] Success: {success_count}, Errors: {error_count}")
        
        if error_count > 0:
            logger.warning("Error details:")
            for idx, key, msg, err_type in error_details:
                logger.warning(f"  Key {idx} ({key}): {err_type} - {msg}")
        
    except Exception as e:
        logger.error(f"[FAIL] Data integrity test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise
    finally:
        if env is not None:
            try:
                env.close()
            except:
                pass

def test_pickle_windows_specific(lmdb_path):
    """测试Windows特定的pickle问题"""
    logger.info("\n" + "=" * 60)
    logger.info("Test 3: Windows-Specific Pickle Issues")
    logger.info("=" * 60)
    
    if platform.system() != 'Windows':
        logger.info("Not Windows, skipping Windows-specific tests")
        return
    
    env = None
    try:
        env = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=256,
        )
        
        with env.begin() as txn:
            keys = list(txn.cursor().iternext(values=False))
        
        if len(keys) == 0:
            logger.warning("No keys found")
            return
        
        # 测试路径相关的问题
        logger.info("Test 3.1: Testing path-related issues...")
        first_key = keys[0]
        with env.begin() as txn:
            data = txn.get(first_key)
        
        if data:
            try:
                obj = pickle.loads(data)
                logger.info("[OK] No path-related issues detected")
                
                # 检查对象中是否包含路径
                if isinstance(obj, dict):
                    for key, value in obj.items():
                        if isinstance(value, str) and (os.sep in value or ':' in value):
                            logger.debug(f"  Found path-like string in key '{key}': {value[:50]}...")
            except Exception as e:
                if 'path' in str(e).lower() or 'file' in str(e).lower():
                    logger.error(f"[FAIL] Path-related error: {e}")
                else:
                    logger.warning(f"Error (may not be path-related): {e}")
        
        # 测试编码问题
        logger.info("Test 3.2: Testing encoding issues...")
        try:
            obj = pickle.loads(data)
            logger.info("[OK] No encoding issues detected")
        except UnicodeDecodeError as e:
            logger.error(f"[FAIL] Unicode decode error: {e}")
            logger.error("  This might be a Windows encoding issue")
        except Exception as e:
            if 'encoding' in str(e).lower() or 'decode' in str(e).lower():
                logger.error(f"[FAIL] Encoding-related error: {e}")
        
        # 测试大对象
        logger.info("Test 3.3: Testing large objects...")
        large_keys = []
        for key in keys[:10]:
            with env.begin() as txn:
                data = txn.get(key)
            if data and len(data) > 1024 * 1024:  # > 1MB
                large_keys.append((key, len(data)))
        
        if large_keys:
            logger.info(f"Found {len(large_keys)} large objects (>1MB)")
            for key, size in large_keys:
                try:
                    with env.begin() as txn:
                        data = txn.get(key)
                    obj = pickle.loads(data)
                    logger.info(f"  Key {key}: OK ({size / 1024 / 1024:.2f}MB)")
                except Exception as e:
                    logger.error(f"  Key {key}: FAILED - {e}")
        else:
            logger.info("No large objects found")
        
    except Exception as e:
        logger.error(f"[FAIL] Windows-specific test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise
    finally:
        if env is not None:
            try:
                env.close()
            except:
                pass

def test_pickle_error_recovery(lmdb_path):
    """测试pickle错误恢复"""
    logger.info("\n" + "=" * 60)
    logger.info("Test 4: Pickle Error Recovery")
    logger.info("=" * 60)
    
    env = None
    try:
        env = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=256,
        )
        
        with env.begin() as txn:
            keys = list(txn.cursor().iternext(values=False))
        
        logger.info("Testing error recovery strategies...")
        
        # 测试策略1: 使用不同的pickle协议
        logger.info("Strategy 1: Try different pickle protocols...")
        first_key = keys[0]
        with env.begin() as txn:
            data = txn.get(first_key)
        
        if data:
            # 先尝试正常反序列化
            try:
                obj = pickle.loads(data)
                logger.info("[OK] Normal deserialization works")
            except Exception as e:
                logger.warning(f"Normal deserialization failed: {e}")
                # 尝试使用不同的协议重新序列化（如果可能）
                logger.info("  Attempting recovery...")
                # 这里无法恢复，因为原始数据已经损坏
        
        # 测试策略2: 检查数据完整性
        logger.info("Strategy 2: Check data integrity...")
        corrupted_count = 0
        for i, key in enumerate(keys[:10]):
            with env.begin() as txn:
                data = txn.get(key)
            
            if data:
                # 检查基本完整性
                if len(data) < 10:  # 太小的数据可能有问题
                    logger.warning(f"  Key {i+1}: Data too small ({len(data)} bytes)")
                    corrupted_count += 1
                elif not data.startswith(b'\x80') and not data.startswith(b'('):  # 不是pickle格式
                    logger.warning(f"  Key {i+1}: Data doesn't look like pickle format")
                    corrupted_count += 1
        
        if corrupted_count == 0:
            logger.info("[OK] All data appears to be valid pickle format")
        else:
            logger.warning(f"Found {corrupted_count} potentially corrupted entries")
        
    except Exception as e:
        logger.error(f"[FAIL] Error recovery test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise
    finally:
        if env is not None:
            try:
                env.close()
            except:
                pass

def main():
    logger.info("=" * 60)
    logger.info("Pickle Deserialization Test")
    logger.info("=" * 60)
    logger.info(f"Platform: {platform.system()} {platform.release()}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Pickle protocol: {pickle.HIGHEST_PROTOCOL}")
    
    # 测试路径
    test_output_dir = "./test_output"
    lmdb_path = os.path.join(test_output_dir, "ligand_predict.lmdb")
    
    if not os.path.exists(lmdb_path):
        logger.error(f"LMDB file not found: {lmdb_path}")
        return
    
    try:
        test_pickle_protocol_compatibility(lmdb_path)
        test_pickle_data_integrity(lmdb_path)
        test_pickle_windows_specific(lmdb_path)
        test_pickle_error_recovery(lmdb_path)
        
        logger.info("\n" + "=" * 60)
        logger.info("All pickle tests completed!")
        logger.info("=" * 60)
        
    except Exception as e:
        logger.error(f"\nTest suite failed: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

