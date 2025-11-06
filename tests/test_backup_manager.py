"""
Tests for backup manager.
"""
import pytest
import tempfile
import shutil
from pathlib import Path
from compass.service.utils.backup import BackupManager, BackupConfig


@pytest.fixture
def backup_config():
    """Backup configuration for testing."""
    return BackupConfig(
        enabled=True,
        backup_dir="./test_backups",
        max_backups=3,
        compress=False  # Disable compression for faster tests
    )


@pytest.fixture
def temp_dirs():
    """Create temporary directories for testing."""
    data_dir = tempfile.mkdtemp()
    checkpoint_dir = tempfile.mkdtemp()
    yield data_dir, checkpoint_dir
    shutil.rmtree(data_dir, ignore_errors=True)
    shutil.rmtree(checkpoint_dir, ignore_errors=True)


def test_backup_manager_init(backup_config):
    """Test BackupManager initialization."""
    manager = BackupManager(backup_config)
    assert manager.config.enabled is True
    assert manager.backup_dir.exists()


def test_backup_manager_create_backup(backup_config, temp_dirs):
    """Test creating a backup."""
    data_dir, checkpoint_dir = temp_dirs
    
    # Create test files
    (Path(data_dir) / "test.txt").write_text("test data")
    (Path(checkpoint_dir) / "model.pth").write_text("test model")
    
    manager = BackupManager(backup_config)
    backup_path = manager.create_backup(
        datasets_dir=data_dir,
        models_dir=checkpoint_dir
    )
    
    assert backup_path is not None
    assert Path(backup_path).exists()
    
    # Cleanup
    shutil.rmtree(backup_config.backup_dir, ignore_errors=True)


def test_backup_manager_list_backups(backup_config):
    """Test listing backups."""
    manager = BackupManager(backup_config)
    
    # Create a backup
    backup_path = manager.create_backup()
    
    # List backups
    backups = manager.list_backups()
    assert len(backups) > 0
    
    # Cleanup
    shutil.rmtree(backup_config.backup_dir, ignore_errors=True)


def test_backup_manager_cleanup_old_backups(backup_config):
    """Test cleanup of old backups."""
    manager = BackupManager(backup_config)
    
    # Create multiple backups
    for i in range(5):
        manager.create_backup()
    
    # Should only keep max_backups
    backups = manager.list_backups()
    assert len(backups) <= backup_config.max_backups
    
    # Cleanup
    shutil.rmtree(backup_config.backup_dir, ignore_errors=True)


def test_backup_manager_disabled(backup_config):
    """Test backup manager when disabled."""
    backup_config.enabled = False
    manager = BackupManager(backup_config)
    
    backup_path = manager.create_backup()
    assert backup_path is None


