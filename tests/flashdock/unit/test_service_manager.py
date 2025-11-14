"""
单元测试：ServiceManager服务管理器
"""
import pytest
import sys
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "services"))

from service_manager import ServiceManager
from registry_client import FlashDockRegistryClient


class TestServiceManager:
    """ServiceManager测试类"""
    
    def test_init_default(self):
        """测试默认初始化"""
        with patch('service_manager.FlashDockRegistryClient') as mock_registry_class:
            mock_registry = Mock()
            mock_registry.discover_compass_services.return_value = []
            mock_registry_class.return_value = mock_registry
            
            manager = ServiceManager()
            
            assert manager.registry_client is not None
            assert manager._services_cache == []
    
    def test_init_custom(self):
        """测试自定义初始化"""
        with patch('service_manager.FlashDockRegistryClient') as mock_registry_class:
            mock_registry = Mock()
            mock_registry.discover_compass_services.return_value = []
            mock_registry_class.return_value = mock_registry
            
            manager = ServiceManager(
                registry_url="http://test:9000",
                timeout=10.0
            )
            
            mock_registry_class.assert_called_once_with(
                registry_url="http://test:9000",
                timeout=10.0
            )
    
    def test_refresh_services_success(self):
        """测试成功刷新服务"""
        mock_service1 = Mock()
        mock_service1.service_id = "compass-1"
        mock_service2 = Mock()
        mock_service2.service_id = "compass-2"
        
        with patch('service_manager.FlashDockRegistryClient') as mock_registry_class:
            mock_registry = Mock()
            mock_registry.discover_compass_services.return_value = [mock_service1, mock_service2]
            mock_registry_class.return_value = mock_registry
            
            manager = ServiceManager()
            manager.refresh_services()
            
            assert len(manager._services_cache) == 2
            assert manager._services_cache[0].service_id == "compass-1"
            assert manager._services_cache[1].service_id == "compass-2"
            mock_registry.discover_compass_services.assert_called_with(healthy_only=False)
    
    def test_refresh_services_empty_list(self):
        """测试刷新服务返回空列表"""
        with patch('service_manager.FlashDockRegistryClient') as mock_registry_class:
            mock_registry = Mock()
            mock_registry.discover_compass_services.return_value = []
            mock_registry_class.return_value = mock_registry
            
            manager = ServiceManager()
            manager.refresh_services()
            
            assert manager._services_cache == []
    
    def test_refresh_services_exception(self):
        """测试刷新服务时发生异常"""
        with patch('service_manager.FlashDockRegistryClient') as mock_registry_class:
            mock_registry = Mock()
            mock_registry.discover_compass_services.side_effect = Exception("Error")
            mock_registry_class.return_value = mock_registry
            
            manager = ServiceManager()
            # 初始缓存为空
            manager._services_cache = []
            manager.refresh_services()
            
            # 异常时应该保持空列表
            assert manager._services_cache == []
    
    def test_refresh_services_exception_preserves_cache(self):
        """测试刷新服务异常时保留现有缓存"""
        existing_service = Mock()
        existing_service.service_id = "compass-existing"
        
        with patch('service_manager.FlashDockRegistryClient') as mock_registry_class:
            mock_registry = Mock()
            mock_registry.discover_compass_services.side_effect = Exception("Error")
            mock_registry_class.return_value = mock_registry
            
            manager = ServiceManager()
            # 设置现有缓存
            manager._services_cache = [existing_service]
            manager.refresh_services()
            
            # 异常时应该保留现有缓存
            assert len(manager._services_cache) == 1
            assert manager._services_cache[0].service_id == "compass-existing"

