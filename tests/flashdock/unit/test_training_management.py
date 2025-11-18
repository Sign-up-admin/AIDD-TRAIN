"""
单元测试：training_management页面模块
"""
import pytest
import sys
import json
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "pages"))

# Import functions from training_management
import importlib.util
training_management_path = project_root / "FLASH_DOCK-main" / "pages" / "training_management.py"
spec = importlib.util.spec_from_file_location("training_management", training_management_path)
training_management = importlib.util.module_from_spec(spec)
spec.loader.exec_module(training_management)

from compass_client import CompassError


class TestTrainingManagement:
    """TrainingManagement页面模块测试类"""
    
    def test_format_error_message_compass_error(self):
        """测试格式化CompassError错误消息"""
        error = CompassError(
            message="Test error",
            status_code=400,
            error_code="ERR_1000",
            detail={"field": "config", "message": "Invalid config"}
        )
        
        formatted = training_management._format_error_message(error)
        
        assert "Test error" in formatted
        assert "ERR_1000" in formatted
        assert "400" in formatted
        assert "config" in formatted
    
    def test_format_error_message_compass_error_minimal(self):
        """测试格式化最小CompassError错误消息"""
        error = CompassError(message="Simple error")
        
        formatted = training_management._format_error_message(error)
        
        assert "Simple error" in formatted
    
    def test_format_error_message_generic_exception(self):
        """测试格式化通用异常消息"""
        error = ValueError("Generic error")
        
        formatted = training_management._format_error_message(error)
        
        assert "Generic error" in formatted
    
    @patch('streamlit.warning')
    @patch('streamlit.rerun')
    def test_handle_stop_task_not_running(self, mock_rerun, mock_warning):
        """测试停止非运行状态的任务"""
        mock_client = Mock()
        
        # 直接调用函数
        training_management._handle_stop_task(mock_client, "task-123", "completed")
        
        # 验证warning被调用（因为状态不是running或initializing）
        # 注意：如果状态不在["running", "initializing"]中，应该调用warning
        mock_warning.assert_called_once()
        mock_client.stop_training_task.assert_not_called()
    
    @patch('training_management.st')
    @patch('training_management.time')
    def test_handle_stop_task_success(self, mock_time, mock_st):
        """测试成功停止任务"""
        mock_client = Mock()
        mock_client.stop_training_task.return_value = {
            "message": "Task stopped successfully",
            "task_id": "task-123"
        }
        mock_client.get_training_task.return_value = {
            "task_id": "task-123",
            "status": "cancelled"
        }
        
        mock_st.spinner = Mock()
        mock_st.spinner.return_value.__enter__ = Mock()
        mock_st.spinner.return_value.__exit__ = Mock(return_value=None)
        mock_st.success = Mock()
        mock_st.rerun = Mock()
        mock_time.sleep = Mock()
        
        training_management._handle_stop_task(mock_client, "task-123", "running")
        
        mock_client.stop_training_task.assert_called_once_with("task-123")
        mock_client.get_training_task.assert_called()
    
    @patch('streamlit.spinner')
    @patch('streamlit.error')
    @patch('streamlit.markdown')
    @patch('streamlit.info')
    @patch('training_management.time.sleep')
    def test_handle_stop_task_compass_error(self, mock_sleep, mock_info, mock_markdown, mock_error, mock_spinner):
        """测试停止任务时发生CompassError"""
        mock_client = Mock()
        error = CompassError(
            message="Cannot stop task",
            status_code=400,
            error_code="ERR_1000"
        )
        mock_client.stop_training_task.side_effect = error
        mock_client.get_training_task = Mock()  # 添加这个mock，因为错误处理中可能会调用
        
        # 设置spinner为context manager
        mock_spinner_context = Mock()
        mock_spinner_context.__enter__ = Mock(return_value=None)
        mock_spinner_context.__exit__ = Mock(return_value=None)
        mock_spinner.return_value = mock_spinner_context
        
        # 调用函数，应该捕获CompassError并显示错误
        training_management._handle_stop_task(mock_client, "task-123", "running")
        
        # 验证错误处理被调用 - CompassError应该在except块中被捕获
        # 检查error或markdown是否被调用（因为_format_error_message会调用markdown）
        assert mock_error.called or mock_markdown.called

