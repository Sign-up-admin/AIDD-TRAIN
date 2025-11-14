"""
集成测试：服务注册中心集成
"""
import pytest
import sys
from pathlib import Path
from unittest.mock import Mock, patch

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "services"))

from registry_client import FlashDockRegistryClient
from service_manager import ServiceManager
from services.common.service_protocol import ServiceInfo


@pytest.mark.integration
class TestRegistryIntegration:
    """服务注册中心集成测试类"""
    
    @patch('registry_client.RegistryClient')
    def test_service_discovery_flow(self, mock_registry_client_class):
        """测试服务发现流程"""
        # 创建模拟服务
        mock_service1 = Mock(spec=ServiceInfo)
        mock_service1.service_id = "compass-1"
        mock_service1.host = "localhost"
        mock_service1.port = 8080
        mock_service1.status = Mock(value="healthy")
        
        mock_service2 = Mock(spec=ServiceInfo)
        mock_service2.service_id = "compass-2"
        mock_service2.host = "localhost"
        mock_service2.port = 8081
        mock_service2.status = Mock(value="unhealthy")
        
        # 设置mock
        mock_registry_client = Mock()
        mock_registry_client.discover_services.return_value = [mock_service1, mock_service2]
        mock_registry_client.check_registry_health.return_value = True
        mock_registry_client_class.return_value = mock_registry_client
        
        # 测试FlashDockRegistryClient
        client = FlashDockRegistryClient()
        
        # 测试发现所有服务
        all_services = client.discover_compass_services(healthy_only=False)
        assert len(all_services) == 2
        
        # 测试只发现健康服务
        healthy_services = client.discover_compass_services(healthy_only=True)
        assert len(healthy_services) == 1
        assert healthy_services[0].service_id == "compass-1"
        
        # 测试健康检查
        assert client.check_registry_available() is True
    
    @patch('registry_client.RegistryClient')
    def test_service_manager_integration(self, mock_registry_client_class):
        """测试ServiceManager与服务注册中心集成"""
        mock_service = Mock(spec=ServiceInfo)
        mock_service.service_id = "compass-1"
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        
        mock_registry_client = Mock()
        mock_registry_client.discover_services.return_value = [mock_service]
        mock_registry_client_class.return_value = mock_registry_client
        
        # 测试ServiceManager
        manager = ServiceManager()
        
        # 刷新服务
        manager.refresh_services()
        assert len(manager._services_cache) == 1
        assert manager._services_cache[0].service_id == "compass-1"
        
        # 再次刷新
        manager.refresh_services()
        assert len(manager._services_cache) == 1

