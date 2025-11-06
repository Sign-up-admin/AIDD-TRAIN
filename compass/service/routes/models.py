"""
Model management routes.
"""
from fastapi import APIRouter, HTTPException, status
from fastapi.responses import FileResponse
from typing import Optional
import os

from compass.service.models.model import ModelListResponse, ModelResponse
from compass.service.services.model_service import ModelService

router = APIRouter()
model_service = ModelService()


@router.get("/api/v1/models", response_model=ModelListResponse)
async def list_models():
    """List all models."""
    models = model_service.list_models()
    return ModelListResponse(models=models, count=len(models))


@router.get("/api/v1/models/{model_id}", response_model=ModelResponse)
async def get_model(model_id: str):
    """Get model by ID."""
    model = model_service.get_model(model_id)
    if not model:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Model {model_id} not found"
        )
    return model


@router.get("/api/v1/models/{model_id}/download")
async def download_model(model_id: str):
    """Download model checkpoint."""
    model = model_service.get_model(model_id)
    if not model:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Model {model_id} not found"
        )
    
    checkpoint_path = model.checkpoint_path
    if not os.path.exists(checkpoint_path):
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Checkpoint file not found: {checkpoint_path}"
        )
    
    return FileResponse(
        checkpoint_path,
        filename=os.path.basename(checkpoint_path),
        media_type='application/octet-stream'
    )


@router.delete("/api/v1/models/{model_id}")
async def delete_model(model_id: str):
    """Delete a model."""
    success = model_service.delete_model(model_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Model {model_id} not found"
        )
    return {"message": f"Model {model_id} deleted", "model_id": model_id}


