"""
Model models for model service.
"""
from pydantic import BaseModel, Field, validator
from typing import Dict, Optional, List
from datetime import datetime
from pathlib import Path


class ModelResponse(BaseModel):
    """Response model for model."""
    model_id: str
    name: str
    version: str
    task_id: Optional[str] = None
    checkpoint_path: str
    file_size: int
    created_at: datetime
    metadata: Dict = Field(default_factory=dict)
    metrics: Dict = Field(default_factory=dict)


class ModelListResponse(BaseModel):
    """Response model for model list."""
    models: List[ModelResponse]
    count: int


class InferenceRequest(BaseModel):
    """Request model for inference."""
    protein_path: str
    ligand_path: str
    model_id: Optional[str] = None  # If None, use latest model
    
    @validator('protein_path')
    def validate_protein_path(cls, v):
        """Validate protein file path."""
        if not isinstance(v, str) or len(v) == 0:
            raise ValueError("protein_path is required and must be a non-empty string")
        # Check file extension
        path = Path(v)
        if path.suffix.lower() not in ['.pdb', '.pdbqt']:
            raise ValueError("protein_path must be a PDB or PDBQT file")
        return v
    
    @validator('ligand_path')
    def validate_ligand_path(cls, v):
        """Validate ligand file path."""
        if not isinstance(v, str) or len(v) == 0:
            raise ValueError("ligand_path is required and must be a non-empty string")
        # Check file extension
        path = Path(v)
        if path.suffix.lower() not in ['.sdf', '.mol2', '.mol', '.pdb']:
            raise ValueError("ligand_path must be a SDF, MOL2, MOL, or PDB file")
        return v
    
    @validator('model_id')
    def validate_model_id(cls, v):
        """Validate model ID format."""
        if v is not None and not isinstance(v, str):
            raise ValueError("model_id must be a string")
        if v is not None and len(v) == 0:
            raise ValueError("model_id cannot be empty")
        return v


class InferenceResponse(BaseModel):
    """Response model for inference."""
    binding_affinity: float
    confidence: Optional[float] = None
    model_id: str
    inference_time: float


class BatchInferenceRequest(BaseModel):
    """Request model for batch inference."""
    protein_ligand_pairs: List[Dict[str, str]]  # List of {"protein": path, "ligand": path}
    model_id: Optional[str] = None
    
    @validator('protein_ligand_pairs')
    def validate_pairs(cls, v):
        """Validate protein-ligand pairs."""
        if not isinstance(v, list):
            raise ValueError("protein_ligand_pairs must be a list")
        if len(v) == 0:
            raise ValueError("protein_ligand_pairs cannot be empty")
        if len(v) > 1000:
            raise ValueError("protein_ligand_pairs cannot exceed 1000 pairs")
        
        for idx, pair in enumerate(v):
            if not isinstance(pair, dict):
                raise ValueError(f"Pair {idx} must be a dictionary")
            if 'protein' not in pair or 'ligand' not in pair:
                raise ValueError(f"Pair {idx} must contain 'protein' and 'ligand' keys")
            if not isinstance(pair['protein'], str) or not isinstance(pair['ligand'], str):
                raise ValueError(f"Pair {idx} protein and ligand paths must be strings")
        
        return v
    
    @validator('model_id')
    def validate_model_id(cls, v):
        """Validate model ID format."""
        if v is not None and not isinstance(v, str):
            raise ValueError("model_id must be a string")
        if v is not None and len(v) == 0:
            raise ValueError("model_id cannot be empty")
        return v


class BatchInferenceResponse(BaseModel):
    """Response model for batch inference."""
    results: List[Dict]
    total_count: int
    success_count: int
    failed_count: int

