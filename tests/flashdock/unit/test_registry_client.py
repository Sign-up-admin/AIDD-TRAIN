"""
单元测试：FlashDockRegistryClient服务注册中心客户端
"""
import pytest
import sys
import requests
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "services"))

from registry_client import FlashDockRegistryClient
from services.common.service_protocol import ServiceInfo, ServiceStatus


class TestFlashDockRegistryClient:
    """FlashDockRegistryClient测试类"""
    
    def test_init_default(self):
        """测试默认初始化"""
        client = FlashDockRegistryClient()
        assert client.registry_url == "http://localhost:8500"
        assert client.timeout == 5.0
        assert client.client is not None
    
    def test_init_custom(self):
        """测试自定义初始化"""
        client = FlashDockRegistryClient(
            registry_url="http://test:9000",
            timeout=10.0
        )
        assert client.registry_url == "http://test:9000"
        assert client.timeout == 10.0
    
    @patch('registry_client.RegistryClient')
    def test_discover_compass_services_success(self, mock_registry_client_class):
        """测试成功发现服务"""
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
        mock_service2.status = Mock(value="healthy")
        
        # 设置mock
        mock_registry_client = Mock()
        mock_registry_client.discover_services.return_value = [mock_service1, mock_service2]
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        services = client.discover_compass_services(healthy_only=True)
        
        assert len(services) == 2
        assert services[0].service_id == "compass-1"
        assert services[1].service_id == "compass-2"
        mock_registry_client.discover_services.assert_called_once_with(
            service_name="compass-service",
            status_filter="healthy"
        )
    
    @patch('registry_client.RegistryClient')
    def test_discover_compass_services_healthy_only_false(self, mock_registry_client_class):
        """测试发现所有服务（包括不健康的）"""
        mock_registry_client = Mock()
        mock_registry_client.discover_services.return_value = []
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        services = client.discover_compass_services(healthy_only=False)
        
        assert services == []
        mock_registry_client.discover_services.assert_called_once_with(
            service_name="compass-service",
            status_filter=None
        )
    
    @patch('registry_client.RegistryClient')
    def test_discover_compass_services_timeout(self, mock_registry_client_class):
        """测试超时错误处理"""
        mock_registry_client = Mock()
        mock_registry_client.discover_services.side_effect = requests.exceptions.Timeout("Timeout")
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        services = client.discover_compass_services()
        
        assert services == []
    
    @patch('registry_client.RegistryClient')
    def test_discover_compass_services_connection_error(self, mock_registry_client_class):
        """测试连接错误处理"""
        mock_registry_client = Mock()
        mock_registry_client.discover_services.side_effect = requests.exceptions.ConnectionError("Connection error")
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        services = client.discover_compass_services()
        
        assert services == []
    
    @patch('registry_client.RegistryClient')
    def test_discover_compass_services_request_exception(self, mock_registry_client_class):
        """测试请求异常处理"""
        mock_registry_client = Mock()
        mock_registry_client.discover_services.side_effect = requests.exceptions.RequestException("Request error")
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        services = client.discover_compass_services()
        
        assert services == []
    
    @patch('registry_client.RegistryClient')
    def test_discover_compass_services_unexpected_error(self, mock_registry_client_class):
        """测试意外错误处理"""
        mock_registry_client = Mock()
        mock_registry_client.discover_services.side_effect = ValueError("Unexpected error")
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        services = client.discover_compass_services()
        
        assert services == []
    
    @patch('registry_client.RegistryClient')
    def test_check_registry_available_success(self, mock_registry_client_class):
        """测试检查注册中心可用性（成功）"""
        mock_registry_client = Mock()
        mock_registry_client.check_registry_health.return_value = True
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        result = client.check_registry_available()
        
        assert result is True
        mock_registry_client.check_registry_health.assert_called_once()
    
    @patch('registry_client.RegistryClient')
    def test_check_registry_available_failure(self, mock_registry_client_class):
        """测试检查注册中心可用性（失败）"""
        mock_registry_client = Mock()
        mock_registry_client.check_registry_health.return_value = False
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        result = client.check_registry_available()
        
        assert result is False
    
    @patch('registry_client.RegistryClient')
    def test_check_registry_available_exception(self, mock_registry_client_class):
        """测试检查注册中心可用性（异常）"""
        mock_registry_client = Mock()
        mock_registry_client.check_registry_health.side_effect = Exception("Error")
        mock_registry_client_class.return_value = mock_registry_client
        
        client = FlashDockRegistryClient()
        result = client.check_registry_available()
        
        assert result is False

