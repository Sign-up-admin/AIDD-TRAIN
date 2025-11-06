"""
Tests for configuration manager.
"""
import pytest
import os
import tempfile
import json
from pathlib import Path
from compass.service.config_manager import ConfigManager, ServiceConfig


def test_config_manager_defaults():
    """Test ConfigManager with default values."""
    manager = ConfigManager()
    config = manager.config
    
    assert config.service_name == "compass-service"
    assert config.port == 8080
    assert config.max_workers == 4


def test_config_manager_from_env(monkeypatch):
    """Test ConfigManager loading from environment variables."""
    monkeypatch.setenv('COMPASS_PORT', '9090')
    monkeypatch.setenv('COMPASS_MAX_WORKERS', '8')
    
    manager = ConfigManager()
    assert manager.config.port == 9090
    assert manager.config.max_workers == 8
    assert manager.get_source('port') == manager.ConfigSource.ENV


def test_config_manager_from_file():
    """Test ConfigManager loading from JSON file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        config_data = {
            'port': 9090,
            'max_workers': 8,
            'service_name': 'test-service'
        }
        json.dump(config_data, f)
        config_file = f.name
    
    try:
        manager = ConfigManager(config_file=config_file)
        assert manager.config.port == 9090
        assert manager.config.max_workers == 8
        assert manager.config.service_name == 'test-service'
    finally:
        os.unlink(config_file)


def test_config_manager_validation():
    """Test configuration validation."""
    manager = ConfigManager()
    
    # Valid config
    errors = manager.validate()
    assert len(errors) == 0
    
    # Invalid port
    manager.config.port = 70000
    errors = manager.validate()
    assert len(errors) > 0
    assert any('port' in error.lower() for error in errors)


def test_config_manager_get_dict():
    """Test getting configuration as dictionary."""
    manager = ConfigManager()
    config_dict = manager.get_dict()
    
    assert isinstance(config_dict, dict)
    assert 'service_name' in config_dict
    assert 'port' in config_dict
    assert 'max_workers' in config_dict


def test_config_manager_save_to_file():
    """Test saving configuration to file."""
    manager = ConfigManager()
    
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        config_file = f.name
    
    try:
        manager.save_to_file(config_file, format='json')
        assert Path(config_file).exists()
        
        # Verify file contents
        with open(config_file, 'r') as f:
            saved_config = json.load(f)
        assert saved_config['service_name'] == manager.config.service_name
    finally:
        os.unlink(config_file)

