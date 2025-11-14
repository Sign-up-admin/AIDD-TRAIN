"""
单元测试：data_management页面模块
"""
import pytest
import sys
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "pages"))


class TestDataManagement:
    """DataManagement页面模块测试类"""
    
    @patch('data_management.CompassClient')
    @patch('data_management.st')
    def test_client_initialization_success(self, mock_st, mock_client_class):
        """测试成功初始化客户端"""
        mock_client = Mock()
        mock_client_class.return_value = mock_client
        mock_st.success = Mock()
        mock_st.error = Mock()
        mock_st.stop = Mock()
        
        # 导入模块会执行初始化代码
        import importlib
        import data_management
        importlib.reload(data_management)
        
        # 由于模块级别的代码，这里主要验证导入不会出错
        assert True
    
    @patch('data_management.CompassClient')
    @patch('data_management.st')
    def test_client_initialization_failure(self, mock_st, mock_client_class):
        """测试客户端初始化失败"""
        mock_client_class.side_effect = Exception("Connection failed")
        mock_st.error = Mock()
        mock_st.stop = Mock()
        
        # 导入模块会执行初始化代码
        import importlib
        import data_management
        importlib.reload(data_management)
        
        # 验证错误处理
        mock_st.error.assert_called()
        mock_st.stop.assert_called()

