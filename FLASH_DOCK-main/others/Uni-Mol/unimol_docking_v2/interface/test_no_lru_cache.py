#!/usr/bin/env python3
"""
测试禁用lru_cache的版本
检查是否是缓存导致的问题
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
logger = logging.getLogger("test_no_lru")

# 创建禁用lru_cache的LMDBDataset版本
class LMDBDatasetNoCache:
    """LMDBDataset without lru_cache"""
    def __init__(self, db_path):
        self.db_path = db_path
        assert os.path.isfile(self.db_path), "{} not found".format(self.db_path)
        env = self.connect_db(self.db_path)
        with env.begin() as txn:
            self._keys = list(txn.cursor().iternext(values=False))

    def connect_db(self, lmdb_path, save_to_self=False):
        env = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=256,
        )
        if not save_to_self:
            return env
        else:
            self.env = env

    def __len__(self):
        return len(self._keys)

    # 注意：这里移除了@lru_cache装饰器
    def __getitem__(self, idx):
        if not hasattr(self, "env"):
            self.connect_db(self.db_path, save_to_self=True)
        datapoint_pickled = self.env.begin().get(f"{idx}".encode("ascii"))
        data = pickle.loads(datapoint_pickled)
        return data

# 创建禁用lru_cache的TTADockingPoseDataset版本
class TTADockingPoseDatasetNoCache:
    """TTADockingPoseDataset without lru_cache"""
    def __init__(
        self,
        dataset,
        atoms,
        coordinates,
        pocket_atoms,
        pocket_coordinates,
        holo_coordinates,
        holo_pocket_coordinates,
        is_train=True,
        conf_size=10,
    ):
        self.dataset = dataset
        self.atoms = atoms
        self.coordinates = coordinates
        self.pocket_atoms = pocket_atoms
        self.pocket_coordinates = pocket_coordinates
        self.holo_coordinates = holo_coordinates
        self.holo_pocket_coordinates = holo_pocket_coordinates
        self.is_train = is_train
        self.conf_size = conf_size
        self.set_epoch(None)

    def set_epoch(self, epoch, **unused):
        self.epoch = epoch

    def __len__(self):
        return len(self.dataset) * self.conf_size

    # 注意：这里移除了@lru_cache装饰器
    def __getitem__(self, index: int):
        smi_idx = index // self.conf_size
        coord_idx = index % self.conf_size
        atoms = np.array(self.dataset[smi_idx][self.atoms])
        coordinates = np.array(self.dataset[smi_idx][self.coordinates][coord_idx])
        pocket_atoms = np.array(
            [item[0] for item in self.dataset[smi_idx][self.pocket_atoms]]
        )
        pocket_coordinates = np.array(self.dataset[smi_idx][self.pocket_coordinates][0])
        if self.is_train:
            holo_coordinates = np.array(self.dataset[smi_idx][self.holo_coordinates][0])
            holo_pocket_coordinates = np.array(
                self.dataset[smi_idx][self.holo_pocket_coordinates][0]
            )
        else:
            holo_coordinates = coordinates
            holo_pocket_coordinates = pocket_coordinates

        smi = self.dataset[smi_idx]["smi"]
        pocket = self.dataset[smi_idx]["pocket"]

        return {
            "atoms": atoms,
            "coordinates": coordinates.astype(np.float32),
            "pocket_atoms": pocket_atoms,
            "pocket_coordinates": pocket_coordinates.astype(np.float32),
            "holo_coordinates": holo_coordinates.astype(np.float32),
            "holo_pocket_coordinates": holo_pocket_coordinates.astype(np.float32),
            "smi": smi,
            "pocket": pocket,
        }

def test_lmdb_without_cache(lmdb_path):
    """测试不使用lru_cache的LMDBDataset"""
    logger.info("=" * 60)
    logger.info("Test 1: LMDBDataset without lru_cache")
    logger.info("=" * 60)
    
    try:
        logger.info("Creating LMDBDatasetNoCache instance...")
        dataset = LMDBDatasetNoCache(lmdb_path)
        logger.info(f"[OK] Dataset created, size: {len(dataset)}")
        
        if len(dataset) > 0:
            # 测试多次访问
            logger.info("Testing multiple accesses...")
            for i in range(10):
                try:
                    item = dataset[i % len(dataset)]
                    logger.debug(f"  Access {i+1}: OK, type: {type(item)}")
                except Exception as e:
                    logger.error(f"[FAIL] Access {i+1} failed: {e}")
                    logger.error(f"  Error type: {type(e).__name__}")
                    logger.error(f"  Traceback:\n{traceback.format_exc()}")
                    raise
            
            logger.info("[OK] Multiple accesses successful without cache")
        
        # 清理
        if hasattr(dataset, 'env') and dataset.env is not None:
            try:
                dataset.env.close()
                logger.info("[OK] LMDB environment closed")
            except Exception as e:
                logger.warning(f"Failed to close environment: {e}")
        
    except Exception as e:
        logger.error(f"[FAIL] Test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise

def test_tta_without_cache(lmdb_path):
    """测试不使用lru_cache的TTADockingPoseDataset"""
    logger.info("\n" + "=" * 60)
    logger.info("Test 2: TTADockingPoseDataset without lru_cache")
    logger.info("=" * 60)
    
    try:
        # 首先创建基础数据集
        logger.info("Creating base LMDBDataset...")
        base_dataset = LMDBDatasetNoCache(lmdb_path)
        logger.info(f"[OK] Base dataset created, size: {len(base_dataset)}")
        
        if len(base_dataset) == 0:
            logger.warning("Base dataset is empty, skipping TTA test")
            return
        
        # 创建TTA数据集
        logger.info("Creating TTADockingPoseDatasetNoCache...")
        tta_dataset = TTADockingPoseDatasetNoCache(
            base_dataset,
            atoms='atoms',
            coordinates='coordinates',
            pocket_atoms='pocket_atoms',
            pocket_coordinates='pocket_coordinates',
            holo_coordinates='holo_coordinates',
            holo_pocket_coordinates='holo_pocket_coordinates',
            is_train=False,
            conf_size=10,
        )
        logger.info(f"[OK] TTA dataset created, size: {len(tta_dataset)}")
        
        # 测试访问
        if len(tta_dataset) > 0:
            logger.info("Testing TTA dataset access...")
            try:
                first_item = tta_dataset[0]
                logger.info(f"[OK] First item accessed, type: {type(first_item)}")
                if isinstance(first_item, dict):
                    logger.info(f"  Keys: {list(first_item.keys())}")
            except Exception as e:
                logger.error(f"[FAIL] Failed to access first item: {e}")
                logger.error(f"  Error type: {type(e).__name__}")
                logger.error(f"  Traceback:\n{traceback.format_exc()}")
                raise
            
            # 测试多次访问
            logger.info("Testing multiple accesses...")
            for i in range(5):
                try:
                    item = tta_dataset[i % len(tta_dataset)]
                    logger.debug(f"  Access {i+1}: OK")
                except Exception as e:
                    logger.error(f"[FAIL] Access {i+1} failed: {e}")
                    raise
            
            logger.info("[OK] Multiple accesses successful without cache")
        
        # 清理
        if hasattr(base_dataset, 'env') and base_dataset.env is not None:
            try:
                base_dataset.env.close()
            except:
                pass
        
    except Exception as e:
        logger.error(f"[FAIL] Test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise

def test_comparison_with_cache(lmdb_path):
    """对比有缓存和无缓存的性能和行为"""
    logger.info("\n" + "=" * 60)
    logger.info("Test 3: Comparison with and without cache")
    logger.info("=" * 60)
    
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'unimol'))
        from unimol.data import LMDBDataset
        
        # 测试有缓存的版本
        logger.info("Testing with cache...")
        try:
            dataset_with_cache = LMDBDataset(lmdb_path)
            if len(dataset_with_cache) > 0:
                import time
                start = time.time()
                for i in range(10):
                    item = dataset_with_cache[i % len(dataset_with_cache)]
                time_with_cache = time.time() - start
                logger.info(f"[OK] With cache: {time_with_cache:.4f}s for 10 accesses")
                
                if hasattr(dataset_with_cache, 'env') and dataset_with_cache.env is not None:
                    dataset_with_cache.env.close()
        except Exception as e:
            logger.error(f"[FAIL] With cache failed: {e}")
            logger.error(f"  This might indicate a cache-related issue!")
        
        # 测试无缓存的版本
        logger.info("Testing without cache...")
        try:
            dataset_no_cache = LMDBDatasetNoCache(lmdb_path)
            if len(dataset_no_cache) > 0:
                import time
                start = time.time()
                for i in range(10):
                    item = dataset_no_cache[i % len(dataset_no_cache)]
                time_no_cache = time.time() - start
                logger.info(f"[OK] Without cache: {time_no_cache:.4f}s for 10 accesses")
                
                if hasattr(dataset_no_cache, 'env') and dataset_no_cache.env is not None:
                    dataset_no_cache.env.close()
        except Exception as e:
            logger.error(f"[FAIL] Without cache failed: {e}")
            logger.error(f"  This might indicate a non-cache-related issue!")
        
    except Exception as e:
        logger.error(f"[FAIL] Comparison test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")

def main():
    logger.info("=" * 60)
    logger.info("Test: Disabled lru_cache Version")
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
        test_lmdb_without_cache(lmdb_path)
        test_tta_without_cache(lmdb_path)
        test_comparison_with_cache(lmdb_path)
        
        logger.info("\n" + "=" * 60)
        logger.info("All tests completed!")
        logger.info("=" * 60)
        
    except Exception as e:
        logger.error(f"\nTest suite failed: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

