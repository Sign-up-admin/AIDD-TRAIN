# Copyright (c) DP Technology.
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.


import lmdb
import os
import pickle
from functools import lru_cache
import logging
import platform
import threading

logger = logging.getLogger(__name__)

# Windows特定的LMDB配置
_is_windows = platform.system() == 'Windows'

class LMDBDataset:
    def __init__(self, db_path):
        self.db_path = db_path
        assert os.path.isfile(self.db_path), "{} not found".format(self.db_path)
        self._env = None
        self._lock = threading.Lock()  # 用于线程安全
        
        # 在初始化时打开环境并读取键（保持环境打开供后续使用）
        self._open_env()
        try:
            with self._env.begin() as txn:
                self._keys = list(txn.cursor().iternext(values=False))
        except Exception as e:
            # 如果读取键失败，关闭环境
            if self._env is not None:
                try:
                    self._env.close()
                except:
                    pass
                self._env = None
            raise

    def _open_env(self):
        """打开LMDB环境（线程安全）"""
        with self._lock:
            if self._env is None:
                try:
                    # Windows特定的配置
                    max_readers = 1 if _is_windows else 256
                    
                    # 获取文件大小以设置合适的map_size
                    try:
                        file_size = os.path.getsize(self.db_path)
                        # map_size应该至少是文件大小的2倍，最小1GB
                        map_size = max(file_size * 2, 1024 * 1024 * 1024)  # 至少1GB
                    except:
                        # 如果无法获取文件大小，使用默认值
                        map_size = 1024 * 1024 * 1024  # 1GB
                    
                    self._env = lmdb.open(
                        self.db_path,
                        subdir=False,
                        readonly=True,
                        lock=False,
                        readahead=False,
                        meminit=False,
                        max_readers=max_readers,
                        map_size=map_size,  # 明确设置map_size
                        max_dbs=0,  # 不使用命名数据库，更安全
                    )
                    logger.debug(f"LMDB environment opened: {self.db_path}, map_size={map_size}, max_readers={max_readers}")
                except Exception as e:
                    logger.error(f"Failed to open LMDB environment: {e}")
                    logger.error(f"Error type: {type(e).__name__}")
                    raise
            return self._env

    def connect_db(self, lmdb_path, save_to_self=False):
        """保持向后兼容的方法"""
        env = self._open_env()
        if not save_to_self:
            return env
        else:
            # 如果已经保存到self.env，直接返回
            if not hasattr(self, "env") or self.env is None:
                self.env = env
            return env

    def __len__(self):
        return len(self._keys)

    def __getitem__(self, idx):
        # 确保环境已打开
        if self._env is None:
            self._open_env()
        
        # Windows上不使用lru_cache，直接访问以避免线程安全问题
        # 非Windows上可以使用lru_cache以提高性能
        if _is_windows:
            # Windows: 直接访问，不使用缓存
            # 注意：不在锁内访问环境，因为LMDB的begin()本身是线程安全的
            # 但我们需要确保环境不会被关闭
            try:
                # 检查环境是否仍然有效
                if self._env is None:
                    self._open_env()
                
                # LMDB的begin()在只读模式下是线程安全的
                with self._env.begin() as txn:
                    datapoint_pickled = txn.get(f"{idx}".encode("ascii"))
                
                if datapoint_pickled is None:
                    raise KeyError(f"Key {idx} not found in LMDB")
                
                data = pickle.loads(datapoint_pickled)
                return data
            except lmdb.Error as e:
                logger.error(f"LMDB error reading item {idx}: {e}")
                # 如果是环境错误，尝试重新打开
                if "closed" in str(e).lower() or "deleted" in str(e).lower():
                    logger.warning("LMDB environment appears closed, attempting to reopen...")
                    with self._lock:
                        if self._env is not None:
                            try:
                                self._env.close()
                            except:
                                pass
                        self._env = None
                    self._open_env()
                    # 重试一次
                    with self._env.begin() as txn:
                        datapoint_pickled = txn.get(f"{idx}".encode("ascii"))
                    if datapoint_pickled is None:
                        raise KeyError(f"Key {idx} not found in LMDB")
                    data = pickle.loads(datapoint_pickled)
                    return data
                raise
            except Exception as e:
                logger.error(f"Error reading item {idx} from LMDB: {e}")
                logger.error(f"Error type: {type(e).__name__}")
                raise
        else:
            # 非Windows: 使用lru_cache
            return self._getitem_cached(idx)
    
    @lru_cache(maxsize=16)
    def _getitem_cached(self, idx):
        """带缓存的getitem（仅用于非Windows）"""
        try:
            with self._env.begin() as txn:
                datapoint_pickled = txn.get(f"{idx}".encode("ascii"))
            
            if datapoint_pickled is None:
                raise KeyError(f"Key {idx} not found in LMDB")
            
            data = pickle.loads(datapoint_pickled)
            return data
        except Exception as e:
            logger.error(f"Error reading item {idx} from LMDB: {e}")
            raise

    def __del__(self):
        """确保环境被正确关闭"""
        if hasattr(self, '_env') and self._env is not None:
            try:
                self._env.close()
                logger.debug("LMDB environment closed in __del__")
            except Exception as e:
                logger.warning(f"Error closing LMDB environment in __del__: {e}")
