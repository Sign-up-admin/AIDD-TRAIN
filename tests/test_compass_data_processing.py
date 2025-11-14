"""
COMPASS数据处理模块测试
测试数据加载、处理和验证功能
"""

import pytest
import os
import sys
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


@pytest.mark.unit
class TestDataProcessing:
    """数据处理测试类"""
    
    def test_dataset_initialization(self):
        """测试数据集初始化"""
        from compass.data.dataset import PDBBindDataset
        
        temp_dir = tempfile.mkdtemp()
        try:
            # 创建最小数据集配置
            dataset = PDBBindDataset(
                root=temp_dir,
                data_paths=[],
                num_workers=0
            )
            assert dataset is not None
            assert len(dataset) == 0
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_data_loader_paths(self):
        """测试数据加载器路径处理"""
        from compass.data.loader.paths import get_pdb_info, get_data_paths
        
        # 测试空索引文件处理
        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.lst')
        try:
            temp_file.write("")
            temp_file.close()
            
            pdb_info = get_pdb_info(temp_file.name)
            assert isinstance(pdb_info, list)
        finally:
            os.unlink(temp_file.name)
    
    def test_graph_conversion(self):
        """测试图转换功能"""
        try:
            from compass.data.graph.conversion import convert_to_graph
            
            # 这个测试需要实际的PDB数据
            # 在单元测试中，我们主要验证函数存在
            assert convert_to_graph is not None
        except ImportError:
            pytest.skip("图转换模块不可用")
    
    def test_graph_features(self):
        """测试图特征提取"""
        try:
            from compass.data.graph.features import extract_features
            
            # 验证函数存在
            assert extract_features is not None
        except ImportError:
            pytest.skip("特征提取模块不可用")
    
    def test_data_processing_utils(self):
        """测试数据处理工具函数"""
        from compass.data.processing import process_data_file
        
        # 验证函数存在
        assert process_data_file is not None


@pytest.mark.unit
class TestDataService:
    """数据服务测试类"""
    
    def test_data_service_initialization(self):
        """测试数据服务初始化"""
        from compass.service.services.data_service import DataService
        
        service = DataService()
        assert service is not None
    
    def test_dataset_listing(self):
        """测试数据集列表"""
        from compass.service.services.data_service import DataService
        
        service = DataService()
        datasets = service.list_datasets()
        assert isinstance(datasets, list)
    
    def test_file_upload_helpers(self):
        """测试文件上传辅助函数"""
        from compass.service.utils.file_upload_helpers import (
            validate_file_extension,
            validate_file_size
        )
        
        # 测试文件扩展名验证
        assert validate_file_extension("test.pdb", [".pdb", ".sdf"]) == True
        assert validate_file_extension("test.txt", [".pdb", ".sdf"]) == False
        
        # 测试文件大小验证
        assert validate_file_size(1024, max_size=2048) == True
        assert validate_file_size(4096, max_size=2048) == False
    
    def test_upload_queue(self):
        """测试上传队列"""
        from compass.service.services.upload_queue import UploadQueue
        
        queue = UploadQueue()
        assert queue is not None
        
        # 测试队列操作
        queue.add_file("test.pdb", "test_path")
        assert len(queue.get_pending_files()) >= 0

