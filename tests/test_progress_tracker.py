"""
Tests for progress tracker.
"""
import pytest
from compass.service.services.progress_tracker import ProgressTracker


def test_progress_tracker_init():
    """Test ProgressTracker initialization."""
    tracker = ProgressTracker("test-task")
    assert tracker.task_id == "test-task"
    assert tracker.current_stage == "initializing"
    assert tracker.cancelled is False


def test_progress_tracker_set_stage():
    """Test setting stage."""
    tracker = ProgressTracker("test-task")
    tracker.set_stage("training", "Training in progress")
    
    assert tracker.current_stage == "training"
    assert tracker.stage_message == "Training in progress"


def test_progress_tracker_update_training():
    """Test updating training progress."""
    tracker = ProgressTracker("test-task")
    tracker.update_training(
        epoch=5,
        total_epochs=10,
        batch=50,
        total_batches=100,
        train_loss=0.5,
        val_loss=0.4
    )
    
    assert tracker.current_epoch == 5
    assert tracker.total_epochs == 10
    assert tracker.current_batch == 50
    assert tracker.total_batches == 100
    assert tracker.train_loss == 0.5
    assert tracker.val_loss == 0.4


def test_progress_tracker_cancel():
    """Test cancelling tracker."""
    tracker = ProgressTracker("test-task")
    assert tracker.is_cancelled() is False
    
    tracker.cancel()
    assert tracker.is_cancelled() is True
    assert "cancelled" in tracker.stage_message.lower()


def test_progress_tracker_get_progress():
    """Test getting progress information."""
    tracker = ProgressTracker("test-task")
    tracker.update_training(
        epoch=3,
        total_epochs=10,
        batch=30,
        total_batches=100,
        train_loss=0.5
    )
    
    progress = tracker.get_progress()
    assert progress['stage'] == "training"
    assert progress['training']['current_epoch'] == 3
    assert progress['training']['total_epochs'] == 10
    assert progress['cancelled'] is False


