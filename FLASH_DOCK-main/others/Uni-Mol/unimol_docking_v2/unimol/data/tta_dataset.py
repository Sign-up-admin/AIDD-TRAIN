# Copyright (c) DP Technology.
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import numpy as np
from functools import lru_cache
from unicore.data import BaseWrapperDataset
import platform
import threading
from collections import OrderedDict

# Windows特定的缓存实现
_is_windows = platform.system() == 'Windows'

class ThreadSafeLRUCache:
    """线程安全的LRU缓存实现（用于Windows）"""
    def __init__(self, maxsize=128):
        self.maxsize = maxsize
        self.cache = OrderedDict()
        self.lock = threading.Lock()
    
    def get(self, key):
        with self.lock:
            if key in self.cache:
                # 移动到末尾（最近使用）
                value = self.cache.pop(key)
                self.cache[key] = value
                return value
            return None
    
    def put(self, key, value):
        with self.lock:
            if key in self.cache:
                # 更新现有项
                self.cache.pop(key)
            elif len(self.cache) >= self.maxsize:
                # 移除最旧的项
                self.cache.popitem(last=False)
            self.cache[key] = value
    
    def clear(self):
        with self.lock:
            self.cache.clear()

def thread_safe_lru_cache(maxsize=128):
    """线程安全的LRU缓存装饰器"""
    if _is_windows:
        # Windows上使用线程安全的缓存
        def decorator(func):
            cache = ThreadSafeLRUCache(maxsize=maxsize)
            def wrapper(self, *args):
                key = args
                value = cache.get(key)
                if value is None:
                    value = func(self, *args)
                    cache.put(key, value)
                return value
            wrapper.cache_clear = cache.clear
            return wrapper
        return decorator
    else:
        # 非Windows上使用标准的lru_cache
        return lru_cache(maxsize=maxsize)


class TTADockingPoseDataset(BaseWrapperDataset):
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
        
        # Windows上使用线程安全的缓存
        if _is_windows:
            self._cache = ThreadSafeLRUCache(maxsize=16)

    def set_epoch(self, epoch, **unused):
        super().set_epoch(epoch)
        self.epoch = epoch
        # 清除缓存（因为epoch改变了）
        if _is_windows and hasattr(self, '_cache'):
            self._cache.clear()
        elif hasattr(self, '__cached_item__'):
            # 清除lru_cache
            if hasattr(self.__cached_item__, 'cache_clear'):
                self.__cached_item__.cache_clear()

    def __len__(self):
        return len(self.dataset) * self.conf_size

    def _get_item_impl(self, index: int, epoch: int):
        """实际的数据获取实现（不包含缓存）"""
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
        # id = self.dataset[smi_idx]['id']

        return {
            "atoms": atoms,
            "coordinates": coordinates.astype(np.float32),
            "pocket_atoms": pocket_atoms,
            "pocket_coordinates": pocket_coordinates.astype(np.float32),
            "holo_coordinates": holo_coordinates.astype(np.float32),
            "holo_pocket_coordinates": holo_pocket_coordinates.astype(np.float32),
            "smi": smi,
            "pocket": pocket,
            # 'id': id,
        }

    if _is_windows:
        # Windows上使用线程安全的缓存
        def __cached_item__(self, index: int, epoch: int):
            key = (index, epoch)
            value = self._cache.get(key)
            if value is None:
                value = self._get_item_impl(index, epoch)
                self._cache.put(key, value)
            return value
    else:
        # 非Windows上使用标准的lru_cache
        @lru_cache(maxsize=16)
        def __cached_item__(self, index: int, epoch: int):
            return self._get_item_impl(index, epoch)

    def __getitem__(self, index: int):
        return self.__cached_item__(index, self.epoch)
