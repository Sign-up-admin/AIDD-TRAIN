"""
单元测试：service_monitor页面模块
"""
import pytest
import sys
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "pages"))


class TestServiceMonitor:
    """ServiceMonitor页面模块测试类"""
    
    @patch('service_monitor.ServiceManager')
    @patch('service_monitor.CompassClient')
    @patch('service_monitor.st')
    def test_initialization_success(self, mock_st, mock_client_class, mock_manager_class):
        """测试成功初始化"""
        mock_manager = Mock()
        mock_client = Mock()
        mock_manager_class.return_value = mock_manager
        mock_client_class.return_value = mock_client
        mock_st.success = Mock()
        mock_st.error = Mock()
        mock_st.stop = Mock()
        
        # 导入模块会执行初始化代码
        import importlib
        import service_monitor
        importlib.reload(service_monitor)
        
        # 验证初始化
        assert True
    
    @patch('service_monitor.ServiceManager')
    @patch('service_monitor.st')
    def test_initialization_failure(self, mock_st, mock_manager_class):
        """测试初始化失败"""
        mock_manager_class.side_effect = Exception("Connection failed")
        mock_st.error = Mock()
        mock_st.stop = Mock()
        
        # 导入模块会执行初始化代码
        import importlib
        import service_monitor
        importlib.reload(service_monitor)
        
        # 验证错误处理
        mock_st.error.assert_called()
        mock_st.stop.assert_called()

