"""
单元测试：LoadBalancer负载均衡器
"""
import pytest
import sys
from pathlib import Path
from unittest.mock import Mock

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "services"))

from load_balancer import LoadBalancer, LoadBalanceStrategy
from services.common.service_protocol import ServiceInfo, ServiceStatus


@pytest.fixture
def mock_services():
    """创建模拟服务列表"""
    services = []
    for i in range(3):
        service = Mock(spec=ServiceInfo)
        service.service_id = f"compass-{i+1}"
        service.host = "localhost"
        service.port = 8080 + i
        service.status = Mock(value="healthy")
        services.append(service)
    return services


@pytest.fixture
def mock_unhealthy_services():
    """创建包含不健康服务的列表"""
    services = []
    for i in range(3):
        service = Mock(spec=ServiceInfo)
        service.service_id = f"compass-{i+1}"
        service.host = "localhost"
        service.port = 8080 + i
        if i == 0:
            service.status = Mock(value="unhealthy")
        else:
            service.status = Mock(value="healthy")
        services.append(service)
    return services


class TestLoadBalancer:
    """LoadBalancer测试类"""
    
    def test_init_default_strategy(self):
        """测试默认初始化"""
        lb = LoadBalancer()
        assert lb.strategy == LoadBalanceStrategy.ROUND_ROBIN
        assert lb.round_robin_index == 0
        assert lb.connection_counts == {}
    
    def test_init_custom_strategy(self):
        """测试自定义策略初始化"""
        lb = LoadBalancer(strategy=LoadBalanceStrategy.RANDOM)
        assert lb.strategy == LoadBalanceStrategy.RANDOM
    
    def test_select_service_empty_list(self):
        """测试空服务列表"""
        lb = LoadBalancer()
        assert lb.select_service([]) is None
    
    def test_select_service_round_robin(self, mock_services):
        """测试轮询策略"""
        lb = LoadBalancer(strategy=LoadBalanceStrategy.ROUND_ROBIN)
        
        # 第一次选择
        service1 = lb.select_service(mock_services)
        assert service1.service_id == "compass-1"
        
        # 第二次选择
        service2 = lb.select_service(mock_services)
        assert service2.service_id == "compass-2"
        
        # 第三次选择
        service3 = lb.select_service(mock_services)
        assert service3.service_id == "compass-3"
        
        # 第四次选择，应该回到第一个
        service4 = lb.select_service(mock_services)
        assert service4.service_id == "compass-1"
    
    def test_select_service_random(self, mock_services):
        """测试随机策略"""
        lb = LoadBalancer(strategy=LoadBalanceStrategy.RANDOM)
        
        # 多次选择，应该从服务列表中选择
        selected_services = set()
        for _ in range(10):
            service = lb.select_service(mock_services)
            assert service in mock_services
            selected_services.add(service.service_id)
        
        # 应该至少选择过不同的服务
        assert len(selected_services) > 0
    
    def test_select_service_least_connections(self, mock_services):
        """测试最少连接策略"""
        lb = LoadBalancer(strategy=LoadBalanceStrategy.LEAST_CONNECTIONS)
        
        # 初始状态，所有服务连接数为0，应该选择第一个
        service1 = lb.select_service(mock_services)
        assert service1.service_id == "compass-1"
        
        # 增加第一个服务的连接数
        lb.increment_connections("compass-1")
        lb.increment_connections("compass-1")
        
        # 现在应该选择连接数最少的服务
        service2 = lb.select_service(mock_services)
        assert service2.service_id in ["compass-2", "compass-3"]
    
    def test_select_service_filters_unhealthy(self, mock_unhealthy_services):
        """测试过滤不健康服务"""
        lb = LoadBalancer()
        
        # 应该只选择健康的服务
        service = lb.select_service(mock_unhealthy_services)
        assert service.status.value == "healthy"
        assert service.service_id != "compass-1"  # compass-1是不健康的
    
    def test_select_service_fallback_to_unhealthy(self):
        """测试当没有健康服务时回退到不健康服务"""
        services = []
        for i in range(2):
            service = Mock(spec=ServiceInfo)
            service.service_id = f"compass-{i+1}"
            service.status = Mock(value="unhealthy")
            services.append(service)
        
        lb = LoadBalancer()
        service = lb.select_service(services)
        # 应该仍然返回一个服务（虽然不健康）
        assert service is not None
        assert service in services
    
    def test_increment_connections(self):
        """测试增加连接数"""
        lb = LoadBalancer()
        lb.increment_connections("compass-1")
        assert lb.connection_counts["compass-1"] == 1
        
        lb.increment_connections("compass-1")
        assert lb.connection_counts["compass-1"] == 2
    
    def test_decrement_connections(self):
        """测试减少连接数"""
        lb = LoadBalancer()
        lb.connection_counts["compass-1"] = 3
        
        lb.decrement_connections("compass-1")
        assert lb.connection_counts["compass-1"] == 2
        
        lb.decrement_connections("compass-1")
        assert lb.connection_counts["compass-1"] == 1
    
    def test_decrement_connections_below_zero(self):
        """测试连接数不会低于0"""
        lb = LoadBalancer()
        lb.connection_counts["compass-1"] = 1
        
        lb.decrement_connections("compass-1")
        assert lb.connection_counts["compass-1"] == 0
        
        # 再次减少，应该保持为0
        lb.decrement_connections("compass-1")
        assert lb.connection_counts["compass-1"] == 0
    
    def test_decrement_nonexistent_service(self):
        """测试减少不存在的服务的连接数"""
        lb = LoadBalancer()
        # 不应该抛出异常
        lb.decrement_connections("nonexistent")
        assert "nonexistent" not in lb.connection_counts

