"""
Model service implementation.
"""
import os
import uuid
import logging
from typing import Dict, Optional, List
from datetime import datetime
from pathlib import Path

from compass.service.models.model import ModelResponse
from compass.service.config import SERVICE_CONFIG
from compass.training.checkpoint import load_latest_checkpoint_file

logger = logging.getLogger(__name__)


class ModelService:
    """Service for managing models."""
    
    def __init__(self):
        """Initialize model service."""
        self.models: Dict[str, ModelResponse] = {}
        self.checkpoint_dir = Path(SERVICE_CONFIG['checkpoint_dir'])
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self._scan_existing_models()
    
    def _scan_existing_models(self):
        """Scan checkpoint directory for existing models."""
        if not self.checkpoint_dir.exists():
            return
        
        for checkpoint_dir in self.checkpoint_dir.iterdir():
            if checkpoint_dir.is_dir():
                checkpoint_files = list(checkpoint_dir.glob("*.pth.tar"))
                if checkpoint_files:
                    # Use latest checkpoint
                    latest_checkpoint = max(checkpoint_files, key=lambda p: p.stat().st_mtime)
                    self._register_model_from_checkpoint(checkpoint_dir, latest_checkpoint)
    
    def _register_model_from_checkpoint(self, checkpoint_dir: Path, checkpoint_file: Path):
        """Register a model from checkpoint file."""
        model_id = str(uuid.uuid4())
        model_name = checkpoint_dir.name
        
        model = ModelResponse(
            model_id=model_id,
            name=model_name,
            version="1.0.0",
            checkpoint_path=str(checkpoint_file),
            file_size=checkpoint_file.stat().st_size,
            created_at=datetime.fromtimestamp(checkpoint_file.stat().st_mtime),
            metadata={
                'checkpoint_dir': str(checkpoint_dir),
                'checkpoint_file': checkpoint_file.name
            }
        )
        
        self.models[model_id] = model
        logger.info(f"Registered existing model: {model_id} ({model_name})")
    
    def register_model(self, checkpoint_path: str, task_id: Optional[str] = None, metadata: Optional[Dict] = None) -> str:
        """
        Register a new model.
        
        Args:
            checkpoint_path: Path to checkpoint file
            task_id: Optional task ID that created this model
            metadata: Optional metadata
            
        Returns:
            str: Model ID
        """
        checkpoint_file = Path(checkpoint_path)
        if not checkpoint_file.exists():
            raise ValueError(f"Checkpoint file not found: {checkpoint_path}")
        
        model_id = str(uuid.uuid4())
        model_name = checkpoint_file.parent.name
        
        model = ModelResponse(
            model_id=model_id,
            name=model_name,
            version="1.0.0",
            task_id=task_id,
            checkpoint_path=str(checkpoint_file),
            file_size=checkpoint_file.stat().st_size,
            created_at=datetime.now(),
            metadata=metadata or {}
        )
        
        self.models[model_id] = model
        logger.info(f"Registered model: {model_id} ({model_name})")
        return model_id
    
    def get_model(self, model_id: str) -> Optional[ModelResponse]:
        """
        Get model by ID.
        
        Args:
            model_id: Model ID
            
        Returns:
            Optional[ModelResponse]: Model or None if not found
        """
        return self.models.get(model_id)
    
    def list_models(self) -> List[ModelResponse]:
        """
        List all models.
        
        Returns:
            List[ModelResponse]: List of models
        """
        return list(self.models.values())
    
    def get_latest_model(self) -> Optional[ModelResponse]:
        """
        Get latest model.
        
        Returns:
            Optional[ModelResponse]: Latest model or None
        """
        if not self.models:
            return None
        return max(self.models.values(), key=lambda m: m.created_at)
    
    def delete_model(self, model_id: str) -> bool:
        """
        Delete a model.
        
        Args:
            model_id: Model ID
            
        Returns:
            bool: True if deleted successfully
        """
        if model_id not in self.models:
            return False
        
        model = self.models[model_id]
        checkpoint_file = Path(model.checkpoint_path)
        if checkpoint_file.exists():
            # Optionally delete checkpoint file
            # checkpoint_file.unlink()
            pass
        
        del self.models[model_id]
        logger.info(f"Deleted model: {model_id}")
        return True


