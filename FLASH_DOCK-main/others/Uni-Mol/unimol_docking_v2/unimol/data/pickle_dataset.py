# Copyright (c) DP Technology.
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

"""
PickleDataset: 使用pickle文件替代LMDB的Dataset实现
用于避免Windows上LMDB的兼容性问题
"""

import os
import pickle
import logging
import platform
import threading
from pathlib import Path

logger = logging.getLogger(__name__)

_is_windows = platform.system() == 'Windows'


class PickleDataset:
    """
    使用pickle文件替代LMDB的Dataset实现
    
    数据存储格式：
    - 目录结构: {db_path}.pickle_dir/
      - 0.pkl
      - 1.pkl
      - 2.pkl
      - ...
      - _metadata.pkl (存储键列表)
    """
    
    def __init__(self, db_path):
        """
        Args:
            db_path: 原始LMDB文件路径，会自动转换为pickle目录路径
        """
        # 如果传入的是.lmdb文件，转换为.pickle_dir目录
        if db_path.endswith('.lmdb'):
            self.pickle_dir = db_path[:-5] + '.pickle_dir'
        else:
            self.pickle_dir = db_path + '.pickle_dir'
        
        self.db_path = db_path  # 保留原始路径用于兼容性
        self._lock = threading.Lock()  # 用于线程安全
        
        # 检查pickle目录是否存在
        if not os.path.exists(self.pickle_dir):
            # 如果pickle目录不存在，尝试从LMDB转换
            if os.path.isfile(db_path):
                logger.warning(f"Pickle directory not found: {self.pickle_dir}")
                logger.warning(f"LMDB file exists: {db_path}")
                logger.warning("Attempting to convert LMDB to pickle format...")
                self._convert_from_lmdb(db_path, self.pickle_dir)
            else:
                raise FileNotFoundError(f"Neither pickle directory nor LMDB file found: {db_path}")
        
        # 加载元数据
        metadata_path = os.path.join(self.pickle_dir, '_metadata.pkl')
        if os.path.exists(metadata_path):
            with open(metadata_path, 'rb') as f:
                metadata = pickle.load(f)
                self._keys = metadata.get('keys', [])
                self._count = metadata.get('count', len(self._keys))
        else:
            # 如果没有元数据，扫描目录
            self._keys = []
            self._count = 0
            for i in range(1000000):  # 最大支持100万个文件
                pkl_path = os.path.join(self.pickle_dir, f"{i}.pkl")
                if os.path.exists(pkl_path):
                    self._keys.append(i)
                    self._count += 1
                else:
                    break
            # 保存元数据
            self._save_metadata()
        
        logger.info(f"PickleDataset initialized: {self.pickle_dir}, count={self._count}")

    def _convert_from_lmdb(self, lmdb_path, pickle_dir):
        """从LMDB文件转换为pickle格式"""
        try:
            import lmdb
            logger.info(f"Converting LMDB {lmdb_path} to pickle format...")
            
            # 创建pickle目录
            os.makedirs(pickle_dir, exist_ok=True)
            
            # 打开LMDB
            env = lmdb.open(
                lmdb_path,
                subdir=False,
                readonly=True,
                lock=False,
                readahead=False,
                meminit=False,
                max_readers=1 if _is_windows else 256,
            )
            
            keys = []
            with env.begin() as txn:
                # 先获取所有键
                cursor = txn.cursor()
                all_keys = list(cursor.iternext(values=False))
                logger.info(f"Found {len(all_keys)} keys in LMDB")
                
                # 然后逐个转换
                for key in all_keys:
                    try:
                        # 尝试将key转换为整数索引
                        idx = int(key.decode('ascii'))
                        keys.append(idx)
                        
                        # 获取对应的值
                        value = txn.get(key)
                        if value is None:
                            logger.warning(f"Key {key} has no value, skipping")
                            continue
                        
                        # 保存为pickle文件
                        pkl_path = os.path.join(pickle_dir, f"{idx}.pkl")
                        with open(pkl_path, 'wb') as f:
                            f.write(value)  # value已经是pickle数据
                    except Exception as e:
                        logger.warning(f"Failed to convert key {key}: {e}")
                        continue
            
            env.close()
            
            # 保存元数据
            self._keys = sorted(keys)
            self._count = len(keys)
            self._save_metadata()
            
            logger.info(f"Conversion complete: {self._count} items converted")
        except ImportError:
            raise ImportError("LMDB library not available. Cannot convert from LMDB format.")
        except Exception as e:
            logger.error(f"Failed to convert LMDB to pickle: {e}")
            raise

    def _save_metadata(self):
        """保存元数据"""
        metadata_path = os.path.join(self.pickle_dir, '_metadata.pkl')
        metadata = {
            'keys': self._keys,
            'count': self._count
        }
        with open(metadata_path, 'wb') as f:
            pickle.dump(metadata, f)

    def __len__(self):
        return self._count

    def __getitem__(self, idx):
        """获取数据项"""
        if idx < 0 or idx >= self._count:
            raise IndexError(f"Index {idx} out of range [0, {self._count})")
        
        # 使用锁确保线程安全（Windows上特别重要）
        with self._lock:
            pkl_path = os.path.join(self.pickle_dir, f"{idx}.pkl")
            
            if not os.path.exists(pkl_path):
                raise KeyError(f"Pickle file not found: {pkl_path}")
            
            try:
                with open(pkl_path, 'rb') as f:
                    data = pickle.load(f)
                return data
            except Exception as e:
                logger.error(f"Error reading pickle file {pkl_path}: {e}")
                raise

    def connect_db(self, lmdb_path, save_to_self=False):
        """保持向后兼容的方法（返回self）"""
        return self


def convert_lmdb_to_pickle(lmdb_path, pickle_dir=None):
    """
    将LMDB文件转换为pickle格式的辅助函数
    
    Args:
        lmdb_path: LMDB文件路径
        pickle_dir: 输出pickle目录路径（如果为None，自动生成）
    
    Returns:
        PickleDataset实例
    """
    if pickle_dir is None:
        if lmdb_path.endswith('.lmdb'):
            pickle_dir = lmdb_path[:-5] + '.pickle_dir'
        else:
            pickle_dir = lmdb_path + '.pickle_dir'
    
    # 创建临时dataset来触发转换
    dataset = PickleDataset(lmdb_path)
    return dataset

