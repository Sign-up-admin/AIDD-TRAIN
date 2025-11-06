"""
Tests for upload queue manager.
"""
import pytest
import time
import threading
from compass.service.services.upload_queue import UploadQueueManager, UploadStatus


def test_upload_queue_manager_init():
    """Test UploadQueueManager initialization."""
    manager = UploadQueueManager(max_concurrent=3)
    assert manager.max_concurrent == 3
    assert manager.get_active_count() == 0
    assert manager.get_queue_size() == 0


def test_upload_queue_manager_submit():
    """Test submitting upload task."""
    manager = UploadQueueManager(max_concurrent=2)
    
    completed = []
    
    def upload_func(dataset_id: str):
        time.sleep(0.1)
        completed.append(dataset_id)
        return True
    
    task = manager.submit_upload("task-1", "dataset-1", upload_func)
    
    assert task.task_id == "task-1"
    assert task.dataset_id == "dataset-1"
    assert task.status == UploadStatus.PENDING
    
    # Wait for completion
    time.sleep(0.2)
    
    assert "dataset-1" in completed


def test_upload_queue_manager_concurrency():
    """Test concurrency limit."""
    manager = UploadQueueManager(max_concurrent=2)
    
    active_count = []
    
    def upload_func(dataset_id: str):
        active_count.append(manager.get_active_count())
        time.sleep(0.2)
        return True
    
    # Submit 3 tasks
    manager.submit_upload("task-1", "dataset-1", upload_func)
    manager.submit_upload("task-2", "dataset-2", upload_func)
    manager.submit_upload("task-3", "dataset-3", upload_func)
    
    # Wait a bit
    time.sleep(0.1)
    
    # Should not exceed max_concurrent
    assert manager.get_active_count() <= 2


def test_upload_queue_manager_task_status():
    """Test getting task status."""
    manager = UploadQueueManager(max_concurrent=2)
    
    def upload_func(dataset_id: str):
        time.sleep(0.1)
        return True
    
    task = manager.submit_upload("task-1", "dataset-1", upload_func)
    
    # Check status
    status = manager.get_task_status("task-1")
    assert status is not None
    assert status.task_id == "task-1"
    
    # Wait for completion
    time.sleep(0.2)
    
    status = manager.get_task_status("task-1")
    assert status is not None
    assert status.status == UploadStatus.COMPLETED

