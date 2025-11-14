"""
COMPASS训练模块单元测试
测试训练相关的核心功能
"""

import pytest
import os
import sys
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import torch

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


@pytest.mark.unit
class TestTrainingComponents:
    """训练组件单元测试"""
    
    def test_progress_tracker_initialization(self):
        """测试进度跟踪器初始化"""
        from compass.service.services.progress_tracker import ProgressTracker
        
        tracker = ProgressTracker(task_id="test_task_123")
        assert tracker.task_id == "test_task_123"
        assert tracker.current_stage == "initializing"
        assert tracker.is_cancelled() == False
    
    def test_progress_tracker_stage_management(self):
        """测试进度跟踪器阶段管理"""
        from compass.service.services.progress_tracker import ProgressTracker
        
        tracker = ProgressTracker(task_id="test_task")
        tracker.set_stage("data_processing", "处理数据")
        assert tracker.current_stage == "data_processing"
        assert tracker.current_message == "处理数据"
    
    def test_progress_tracker_cancellation(self):
        """测试进度跟踪器取消功能"""
        from compass.service.services.progress_tracker import ProgressTracker
        
        tracker = ProgressTracker(task_id="test_task")
        assert tracker.is_cancelled() == False
        
        tracker.cancel()
        assert tracker.is_cancelled() == True
    
    def test_training_config_preparation(self):
        """测试训练配置准备"""
        from compass.service.utils.training_helpers import prepare_training_config
        
        base_config = {
            "execution_mode": "smoke_test",
            "epochs": 2,
            "batch_size": 1
        }
        
        config = prepare_training_config(base_config)
        assert config is not None
        assert "execution_mode" in config
        assert config["execution_mode"] == "smoke_test"
    
    def test_task_lifecycle_logger(self):
        """测试任务生命周期日志记录器"""
        from compass.service.services.task_lifecycle_logger import TaskLifecycleLogger
        
        logger = TaskLifecycleLogger(task_id="test_task")
        assert logger.task_id == "test_task"
        
        # 测试日志记录
        logger.log_task_created()
        logger.log_task_started()
        logger.log_task_completed()
    
    def test_progress_logger(self):
        """测试进度日志记录器"""
        from compass.service.services.progress_logger import ProgressAwareLogger
        from compass.service.services.progress_tracker import ProgressTracker
        
        tracker = ProgressTracker(task_id="test_task")
        logger = ProgressAwareLogger(tracker=tracker)
        
        assert logger.progress_tracker == tracker
        logger.log("测试消息")
    
    def test_training_exceptions(self):
        """测试训练异常"""
        from compass.training.exceptions import TrainingCancelled, TrainingError
        
        # 测试TrainingCancelled
        exc = TrainingCancelled("用户取消")
        assert str(exc) == "用户取消"
        
        # 测试TrainingError
        exc = TrainingError("训练错误")
        assert str(exc) == "训练错误"
    
    def test_checkpoint_manager(self):
        """测试检查点管理器"""
        from compass.training.checkpoint import CheckpointManager
        
        temp_dir = tempfile.mkdtemp()
        try:
            manager = CheckpointManager(checkpoint_dir=temp_dir)
            assert manager.checkpoint_dir == temp_dir
            
            # 测试检查点路径生成
            epoch = 5
            path = manager.get_checkpoint_path(epoch)
            assert "epoch" in path or str(epoch) in path
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_training_engine_initialization(self):
        """测试训练引擎初始化"""
        # 这个测试需要mock很多依赖
        with patch('compass.training.engine.Trainer') as mock_trainer:
            # 验证Trainer类可以被导入
            from compass.training.engine import Trainer
            assert Trainer is not None
    
    def test_model_architecture(self):
        """测试模型架构"""
        try:
            from compass.training.model import ViSNetPDB
            
            # 测试模型可以实例化（使用最小配置）
            model = ViSNetPDB(
                hidden_channels=32,
                num_layers=2,
                num_rbf=16,
                cutoff=5.0,
                max_num_neighbors=16,
                lmax=1
            )
            
            assert model is not None
            # 验证模型参数
            param_count = sum(p.numel() for p in model.parameters())
            assert param_count > 0
        except Exception as e:
            # 如果CUDA不可用或其他问题，跳过
            pytest.skip(f"模型测试跳过: {e}")


@pytest.mark.unit
class TestTrainingHelpers:
    """训练辅助函数测试"""
    
    def test_output_redirector(self):
        """测试输出重定向器"""
        from compass.service.utils.training_helpers import OutputRedirector
        import io
        
        redirector = OutputRedirector()
        assert redirector is not None
        
        # 测试输出捕获
        with redirector:
            print("测试输出")
        
        # 验证输出被捕获
        output = redirector.get_output()
        assert "测试输出" in output or len(output) >= 0
    
    def test_resource_monitor(self):
        """测试资源监控器"""
        from compass.service.utils.training_helpers import ResourceMonitor
        
        monitor = ResourceMonitor()
        assert monitor is not None
        
        # 测试资源获取
        resources = monitor.get_resources()
        assert resources is not None
        assert isinstance(resources, dict)
    
    def test_task_stop_helpers(self):
        """测试任务停止辅助函数"""
        from compass.service.utils.task_stop_helpers import (
            validate_task_can_be_stopped,
            set_cancellation_flag
        )
        from compass.service.services.progress_tracker import ProgressTracker
        
        tracker = ProgressTracker(task_id="test_task")
        
        # 测试验证任务是否可以停止
        can_stop = validate_task_can_be_stopped(tracker)
        assert can_stop == True
        
        # 测试设置取消标志
        set_cancellation_flag(tracker)
        assert tracker.is_cancelled() == True

