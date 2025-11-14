#!/usr/bin/env python3
"""
测试直接LMDB访问，绕过数据集包装器
验证LMDB环境管理和Windows兼容性
"""
import os
import sys
import platform
import logging
import traceback
import pickle
import lmdb
import gc

# 设置日志
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    stream=sys.stdout,
)
logger = logging.getLogger("test_lmdb")

def test_lmdb_environment_management(lmdb_path):
    """测试LMDB环境管理"""
    logger.info("=" * 60)
    logger.info("Test 1: LMDB Environment Management")
    logger.info("=" * 60)
    
    envs = []
    try:
        # 测试1: 基本打开和关闭
        logger.info("Test 1.1: Basic open and close")
        env1 = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=256,
        )
        envs.append(env1)
        logger.info(f"[OK] LMDB environment opened: {env1}")
        
        # 获取环境信息
        info = env1.info()
        logger.info(f"LMDB info: map_size={info['map_size']}, max_readers={env1.max_readers()}")
        
        # 测试2: 多次打开（模拟多线程场景）
        logger.info("Test 1.2: Multiple opens (simulating multi-thread scenario)")
        for i in range(3):
            env = lmdb.open(
                lmdb_path,
                subdir=False,
                readonly=True,
                lock=False,
                readahead=False,
                meminit=False,
                max_readers=256,
            )
            envs.append(env)
            logger.info(f"[OK] Opened environment {i+1}: {env}")
        
        # 测试3: 正确关闭所有环境
        logger.info("Test 1.3: Closing all environments")
        for i, env in enumerate(envs):
            try:
                env.close()
                logger.info(f"[OK] Closed environment {i+1}")
            except Exception as e:
                logger.error(f"[FAIL] Failed to close environment {i+1}: {e}")
        
        envs.clear()
        logger.info("[OK] All environments closed successfully")
        
    except Exception as e:
        logger.error(f"[FAIL] LMDB environment management test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        # 清理
        for env in envs:
            try:
                env.close()
            except:
                pass
        raise

def test_lmdb_read_operations(lmdb_path):
    """测试LMDB读取操作"""
    logger.info("\n" + "=" * 60)
    logger.info("Test 2: LMDB Read Operations")
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
        
        # 测试2.1: 读取所有键
        logger.info("Test 2.1: Reading all keys")
        with env.begin() as txn:
            keys = list(txn.cursor().iternext(values=False))
            logger.info(f"[OK] Found {len(keys)} keys in LMDB")
        
        if len(keys) == 0:
            logger.warning("No keys found in LMDB, skipping further tests")
            return
        
        # 测试2.2: 读取前几个键的数据
        logger.info("Test 2.2: Reading first few keys")
        for i, key in enumerate(keys[:5]):
            try:
                with env.begin() as txn:
                    data = txn.get(key)
                    logger.info(f"[OK] Key {i+1} ({key}): {len(data) if data else 0} bytes")
            except Exception as e:
                logger.error(f"[FAIL] Failed to read key {i+1} ({key}): {e}")
                raise
        
        # 测试2.3: 多次读取同一个键（测试缓存行为）
        logger.info("Test 2.3: Multiple reads of same key (testing cache behavior)")
        first_key = keys[0]
        for i in range(5):
            with env.begin() as txn:
                data = txn.get(first_key)
                logger.debug(f"Read {i+1}: {len(data) if data else 0} bytes")
        logger.info("[OK] Multiple reads successful")
        
    except Exception as e:
        logger.error(f"[FAIL] LMDB read operations test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise
    finally:
        if env is not None:
            try:
                env.close()
            except:
                pass

def test_pickle_deserialization(lmdb_path):
    """测试pickle反序列化"""
    logger.info("\n" + "=" * 60)
    logger.info("Test 3: Pickle Deserialization")
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
            logger.warning("No keys found, skipping pickle test")
            return
        
        # 测试3.1: 反序列化第一个键
        logger.info("Test 3.1: Deserializing first key")
        first_key = keys[0]
        with env.begin() as txn:
            data = txn.get(first_key)
        
        if data:
            try:
                unpickled = pickle.loads(data)
                logger.info(f"[OK] Pickle deserialization successful")
                logger.info(f"  Type: {type(unpickled)}")
                if isinstance(unpickled, dict):
                    logger.info(f"  Dict keys: {list(unpickled.keys())[:10]}...")
                    # 检查数据类型
                    for key, value in list(unpickled.items())[:5]:
                        logger.debug(f"    {key}: {type(value)}")
            except Exception as e:
                logger.error(f"[FAIL] Pickle deserialization failed: {e}")
                logger.error(f"  Error type: {type(e).__name__}")
                logger.error(f"  Traceback:\n{traceback.format_exc()}")
                raise
        
        # 测试3.2: 反序列化多个键
        logger.info("Test 3.2: Deserializing multiple keys")
        success_count = 0
        for i, key in enumerate(keys[:10]):
            try:
                with env.begin() as txn:
                    data = txn.get(key)
                if data:
                    unpickled = pickle.loads(data)
                    success_count += 1
                    logger.debug(f"  Key {i+1}: OK")
            except Exception as e:
                logger.error(f"[FAIL] Failed to deserialize key {i+1} ({key}): {e}")
                raise
        
        logger.info(f"[OK] Successfully deserialized {success_count} keys")
        
        # 测试3.3: 测试不同pickle协议版本
        logger.info("Test 3.3: Testing pickle protocol compatibility")
        with env.begin() as txn:
            data = txn.get(first_key)
        
        if data:
            # 尝试使用不同的协议版本重新序列化和反序列化
            try:
                obj = pickle.loads(data)
                # 测试不同协议版本
                for protocol in [0, 1, 2, 3, 4, 5]:
                    try:
                        repickled = pickle.dumps(obj, protocol=protocol)
                        reunpickled = pickle.loads(repickled)
                        logger.debug(f"  Protocol {protocol}: OK")
                    except Exception as e:
                        logger.warning(f"  Protocol {protocol} failed: {e}")
                logger.info("[OK] Pickle protocol compatibility test passed")
            except Exception as e:
                logger.error(f"[FAIL] Pickle protocol test failed: {e}")
        
    except Exception as e:
        logger.error(f"[FAIL] Pickle deserialization test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise
    finally:
        if env is not None:
            try:
                env.close()
            except:
                pass

def test_lmdb_without_wrapper(lmdb_path):
    """测试直接使用LMDBDataset（不通过任务包装器）"""
    logger.info("\n" + "=" * 60)
    logger.info("Test 4: Direct LMDBDataset Usage (without task wrapper)")
    logger.info("=" * 60)
    
    try:
        # 导入LMDBDataset
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'unimol'))
        from unimol.data import LMDBDataset
        
        logger.info("Test 4.1: Creating LMDBDataset instance")
        dataset = LMDBDataset(lmdb_path)
        logger.info(f"[OK] LMDBDataset created, size: {len(dataset)}")
        
        # 测试4.2: 直接访问数据项
        logger.info("Test 4.2: Direct dataset item access")
        if len(dataset) > 0:
            try:
                first_item = dataset[0]
                logger.info(f"[OK] Direct access successful, type: {type(first_item)}")
                if isinstance(first_item, dict):
                    logger.info(f"  Keys: {list(first_item.keys())[:10]}...")
            except Exception as e:
                logger.error(f"[FAIL] Direct access failed: {e}")
                logger.error(f"  Error type: {type(e).__name__}")
                logger.error(f"  Traceback:\n{traceback.format_exc()}")
                raise
        
        # 测试4.3: 多次访问（测试lru_cache）
        logger.info("Test 4.3: Multiple accesses (testing lru_cache)")
        for i in range(5):
            try:
                item = dataset[i % len(dataset)]
                logger.debug(f"  Access {i+1}: OK")
            except Exception as e:
                logger.error(f"[FAIL] Access {i+1} failed: {e}")
                raise
        logger.info("[OK] Multiple accesses successful")
        
        # 测试4.4: 检查LMDB环境状态
        logger.info("Test 4.4: Checking LMDB environment state")
        if hasattr(dataset, 'env'):
            logger.info(f"  Environment exists: {dataset.env is not None}")
            if dataset.env is not None:
                try:
                    info = dataset.env.info()
                    logger.info(f"  Environment info: map_size={info['map_size']}")
                except Exception as e:
                    logger.warning(f"  Failed to get environment info: {e}")
        else:
            logger.warning("  Dataset does not have 'env' attribute")
        
        # 清理
        if hasattr(dataset, 'env') and dataset.env is not None:
            try:
                dataset.env.close()
                logger.info("[OK] LMDB environment closed")
            except Exception as e:
                logger.warning(f"Failed to close environment: {e}")
        
    except Exception as e:
        logger.error(f"[FAIL] Direct LMDBDataset test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise

def main():
    logger.info("=" * 60)
    logger.info("Direct LMDB Access Test (Windows Compatibility)")
    logger.info("=" * 60)
    logger.info(f"Platform: {platform.system()} {platform.release()}")
    logger.info(f"Python version: {sys.version}")
    
    # 测试路径
    test_output_dir = "./test_output"
    lmdb_path = os.path.join(test_output_dir, "ligand_predict.lmdb")
    
    if not os.path.exists(lmdb_path):
        logger.error(f"LMDB file not found: {lmdb_path}")
        return
    
    try:
        # 运行所有测试
        test_lmdb_environment_management(lmdb_path)
        test_lmdb_read_operations(lmdb_path)
        test_pickle_deserialization(lmdb_path)
        test_lmdb_without_wrapper(lmdb_path)
        
        logger.info("\n" + "=" * 60)
        logger.info("All tests passed!")
        logger.info("=" * 60)
        
    except Exception as e:
        logger.error(f"\nTest suite failed: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

