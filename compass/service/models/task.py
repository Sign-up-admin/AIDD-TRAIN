"""
Task models for training service.
"""
from pydantic import BaseModel, Field, validator
from typing import Dict, Optional, List, Any
from datetime import datetime
from enum import Enum


class TaskStatus(str, Enum):
    """Training task status."""
    CREATING = "creating"  # Task is being created and validated
    INITIALIZING = "initializing"  # Task is initializing (resources, threads, etc.)
    PENDING = "pending"  # Task is created but not started
    RUNNING = "running"  # Task is running
    PAUSED = "paused"  # Task is paused
    COMPLETED = "completed"  # Task completed successfully
    FAILED = "failed"  # Task failed
    CANCELLED = "cancelled"  # Task was cancelled


class TrainingTaskCreate(BaseModel):
    """Request model for creating a training task."""
    config: Dict = Field(..., description="Training configuration")
    dataset_id: Optional[str] = Field(None, description="Dataset ID to use")
    description: Optional[str] = Field(None, description="Task description")
    
    @validator('config')
    def validate_config(cls, v):
        """Validate training configuration parameters."""
        if not isinstance(v, dict):
            raise ValueError("config must be a dictionary")
        
        # Validate execution_mode
        valid_modes = ['validation_tuned', 'validation', 'prototyping', 'smoke_test', 'production']
        if 'execution_mode' in v and v['execution_mode'] not in valid_modes:
            raise ValueError(f"execution_mode must be one of {valid_modes}")
        
        # Validate epochs
        if 'epochs' in v:
            epochs = v['epochs']
            if not isinstance(epochs, int) or epochs < 1 or epochs > 10000:
                raise ValueError("epochs must be an integer between 1 and 10000")
        
        # Validate batch_size
        if 'batch_size' in v:
            batch_size = v['batch_size']
            if not isinstance(batch_size, int) or batch_size < 1 or batch_size > 128:
                raise ValueError("batch_size must be an integer between 1 and 128")
        
        # Validate learning_rate
        if 'learning_rate' in v:
            lr = v['learning_rate']
            if not isinstance(lr, (int, float)) or lr <= 0 or lr > 1.0:
                raise ValueError("learning_rate must be a positive number between 0 and 1.0")
        
        # Validate optimizer
        valid_optimizers = ['adam', 'adamw', 'sgd', 'rmsprop']
        if 'optimizer' in v and v['optimizer'].lower() not in valid_optimizers:
            raise ValueError(f"optimizer must be one of {valid_optimizers}")
        
        return v
    
    @validator('dataset_id')
    def validate_dataset_id(cls, v):
        """Validate dataset ID format."""
        if v is not None and not isinstance(v, str):
            raise ValueError("dataset_id must be a string")
        if v is not None and len(v) == 0:
            raise ValueError("dataset_id cannot be empty")
        return v


class TrainingTaskResponse(BaseModel):
    """Response model for training task."""
    task_id: str
    status: TaskStatus
    config: Dict
    created_at: datetime
    updated_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    progress: Dict = Field(default_factory=dict)
    error: Optional[str] = None
    description: Optional[str] = None


class TaskListResponse(BaseModel):
    """Response model for task list."""
    tasks: List[TrainingTaskResponse]
    count: int


class TaskLogResponse(BaseModel):
    """Response model for task logs."""
    logs: List[str]
    total_lines: int


class TaskMetricsResponse(BaseModel):
    """Response model for task metrics."""
    metrics: Dict
    timestamp: datetime = Field(default_factory=datetime.now)

