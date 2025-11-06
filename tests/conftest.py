"""
Pytest configuration and fixtures.
"""
import pytest
import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    temp_path = tempfile.mkdtemp()
    yield temp_path
    shutil.rmtree(temp_path, ignore_errors=True)


@pytest.fixture
def mock_service_config(monkeypatch):
    """Mock service configuration for testing."""
    monkeypatch.setenv('COMPASS_HOST', '127.0.0.1')
    monkeypatch.setenv('COMPASS_PORT', '8080')
    monkeypatch.setenv('COMPASS_DATA_DIR', './test_data')
    monkeypatch.setenv('COMPASS_CHECKPOINT_DIR', './test_checkpoints')
    monkeypatch.setenv('COMPASS_LOG_DIR', './test_logs')
    monkeypatch.setenv('REGISTRY_URL', 'http://localhost:8500')


@pytest.fixture
def sample_training_config():
    """Sample training configuration for testing."""
    return {
        'execution_mode': 'smoke_test',
        'epochs': 2,
        'batch_size': 1,
        'learning_rate': 0.001,
        'optimizer': 'adam'
    }


@pytest.fixture
def sample_inference_request():
    """Sample inference request for testing."""
    return {
        'protein_path': 'test_protein.pdb',
        'ligand_path': 'test_ligand.sdf',
        'model_id': None
    }

