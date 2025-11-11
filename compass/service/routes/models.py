"""
Model management routes.
"""

from fastapi import APIRouter, HTTPException, status
from fastapi.responses import FileResponse
from typing import Optional
import os
import logging

from compass.service.models.model import ModelListResponse, ModelResponse
from compass.service.services.model_service import ModelService
from compass.service.exceptions import ServiceException, NotFoundError

router = APIRouter()
model_service = ModelService()
logger = logging.getLogger(__name__)


@router.get("/api/v1/models", response_model=ModelListResponse)
async def list_models():
    """List all models."""
    try:
        models = model_service.list_models()
        return ModelListResponse(models=models, count=len(models))
    except Exception as e:
        logger.error(f"Failed to list models: {e}", exc_info=True)
        raise ServiceException(
            "Failed to list models", status_code=status.HTTP_500_INTERNAL_SERVER_ERROR
        )


@router.get("/api/v1/models/{model_id}", response_model=ModelResponse)
async def get_model(model_id: str):
    """Get model by ID."""
    try:
        model = model_service.get_model(model_id)
        if not model:
            raise NotFoundError("Model", model_id)
        return model
    except NotFoundError:
        raise
    except Exception as e:
        logger.error(f"Failed to get model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to get model {model_id}", status_code=status.HTTP_500_INTERNAL_SERVER_ERROR
        )


@router.get("/api/v1/models/{model_id}/download")
async def download_model(model_id: str):
    """Download model checkpoint."""
    try:
        model = model_service.get_model(model_id)
        if not model:
            raise NotFoundError("Model", model_id)

        checkpoint_path = model.checkpoint_path
        if not os.path.exists(checkpoint_path):
            raise NotFoundError("Checkpoint", checkpoint_path)

        return FileResponse(
            checkpoint_path,
            filename=os.path.basename(checkpoint_path),
            media_type="application/octet-stream",
        )
    except NotFoundError:
        raise
    except Exception as e:
        logger.error(f"Failed to download model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to download model {model_id}", status_code=status.HTTP_500_INTERNAL_SERVER_ERROR
        )


@router.delete("/api/v1/models/{model_id}")
async def delete_model(model_id: str):
    """Delete a model."""
    try:
        success = model_service.delete_model(model_id)
        if not success:
            raise NotFoundError("Model", model_id)
        return {"message": f"Model {model_id} deleted", "model_id": model_id}
    except NotFoundError:
        raise
    except Exception as e:
        logger.error(f"Failed to delete model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to delete model {model_id}", status_code=status.HTTP_500_INTERNAL_SERVER_ERROR
        )
