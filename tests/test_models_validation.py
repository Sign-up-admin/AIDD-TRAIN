"""
Tests for model validation.
"""
import pytest
from pydantic import ValidationError as PydanticValidationError
from compass.service.models.task import TrainingTaskCreate
from compass.service.models.model import InferenceRequest, BatchInferenceRequest


def test_training_task_create_valid():
    """Test valid training task creation."""
    task = TrainingTaskCreate(
        config={
            'execution_mode': 'validation_tuned',
            'epochs': 100,
            'batch_size': 32,
            'learning_rate': 0.001,
            'optimizer': 'adam'
        },
        description="Test task"
    )
    assert task.config['epochs'] == 100


def test_training_task_create_invalid_epochs():
    """Test invalid epochs validation."""
    with pytest.raises(PydanticValidationError):
        TrainingTaskCreate(
            config={
                'epochs': 20000,  # Exceeds max
            }
        )


def test_training_task_create_invalid_batch_size():
    """Test invalid batch size validation."""
    with pytest.raises(PydanticValidationError):
        TrainingTaskCreate(
            config={
                'batch_size': 200,  # Exceeds max
            }
        )


def test_training_task_create_invalid_learning_rate():
    """Test invalid learning rate validation."""
    with pytest.raises(PydanticValidationError):
        TrainingTaskCreate(
            config={
                'learning_rate': 2.0,  # Exceeds max
            }
        )


def test_training_task_create_invalid_optimizer():
    """Test invalid optimizer validation."""
    with pytest.raises(PydanticValidationError):
        TrainingTaskCreate(
            config={
                'optimizer': 'invalid_optimizer',
            }
        )


def test_inference_request_valid():
    """Test valid inference request."""
    request = InferenceRequest(
        protein_path="test.pdb",
        ligand_path="test.sdf"
    )
    assert request.protein_path == "test.pdb"
    assert request.ligand_path == "test.sdf"


def test_inference_request_invalid_protein():
    """Test invalid protein path."""
    with pytest.raises(PydanticValidationError):
        InferenceRequest(
            protein_path="test.txt",  # Invalid extension
            ligand_path="test.sdf"
        )


def test_inference_request_invalid_ligand():
    """Test invalid ligand path."""
    with pytest.raises(PydanticValidationError):
        InferenceRequest(
            protein_path="test.pdb",
            ligand_path="test.txt"  # Invalid extension
        )


def test_batch_inference_request_valid():
    """Test valid batch inference request."""
    request = BatchInferenceRequest(
        protein_ligand_pairs=[
            {'protein': 'test1.pdb', 'ligand': 'test1.sdf'},
            {'protein': 'test2.pdb', 'ligand': 'test2.sdf'}
        ]
    )
    assert len(request.protein_ligand_pairs) == 2


def test_batch_inference_request_empty():
    """Test empty batch inference request."""
    with pytest.raises(PydanticValidationError):
        BatchInferenceRequest(
            protein_ligand_pairs=[]
        )


def test_batch_inference_request_too_large():
    """Test batch inference request exceeding limit."""
    with pytest.raises(PydanticValidationError):
        pairs = [
            {'protein': f'test{i}.pdb', 'ligand': f'test{i}.sdf'}
            for i in range(1001)  # Exceeds max
        ]
        BatchInferenceRequest(protein_ligand_pairs=pairs)


