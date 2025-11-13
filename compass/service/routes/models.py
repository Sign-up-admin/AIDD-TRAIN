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
from compass.service.utils.input_sanitizer import sanitize_string

router = APIRouter()
model_service = ModelService()
logger = logging.getLogger(__name__)


@router.get("/api/v1/models", response_model=ModelListResponse)
async def list_models():
    """List all models."""
    try:
        models = model_service.list_models()
        return ModelListResponse(models=models, count=len(models))
    except (OSError, IOError, PermissionError) as e:
        # File system errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"File system error listing models: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to list models: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    except Exception as e:
        # Unexpected errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"Unexpected error listing models: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to list models: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@router.get("/api/v1/models/{model_id}", response_model=ModelResponse)
async def get_model(model_id: str):
    """Get model by ID."""
    try:
        # Sanitize model_id to prevent injection attacks
        model_id = sanitize_string(model_id, max_length=255)
        model = model_service.get_model(model_id)
        if not model:
            raise NotFoundError("Model", model_id)
        return model
    except NotFoundError:
        raise
    except ValueError as e:
        from compass.service.exceptions import ValidationError

        raise ValidationError(f"Invalid model ID format: {str(e)}")
    except (OSError, IOError, PermissionError) as e:
        # File system errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"File system error getting model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to get model: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    except Exception as e:
        # Unexpected errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"Unexpected error getting model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to get model: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@router.get("/api/v1/models/{model_id}/download")
async def download_model(model_id: str):
    """Download model checkpoint."""
    try:
        # Sanitize model_id to prevent injection attacks
        model_id = sanitize_string(model_id, max_length=255)
        model = model_service.get_model(model_id)
        if not model:
            raise NotFoundError("Model", model_id)

        checkpoint_path = model.checkpoint_path
        # Additional security: validate path doesn't contain path traversal
        if not os.path.exists(checkpoint_path) or ".." in checkpoint_path:
            raise NotFoundError("Checkpoint", checkpoint_path)

        return FileResponse(
            checkpoint_path,
            filename=os.path.basename(checkpoint_path),
            media_type="application/octet-stream",
        )
    except NotFoundError:
        raise
    except ValueError as e:
        from compass.service.exceptions import ValidationError

        raise ValidationError(f"Invalid model ID format: {str(e)}")
    except (OSError, IOError, PermissionError) as e:
        # File system errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"File system error downloading model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to download model: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    except Exception as e:
        # Unexpected errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"Unexpected error downloading model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to download model: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@router.delete("/api/v1/models/{model_id}")
async def delete_model(model_id: str):
    """Delete a model."""
    try:
        # Sanitize model_id to prevent injection attacks
        model_id = sanitize_string(model_id, max_length=255)
        success = model_service.delete_model(model_id)
        if not success:
            raise NotFoundError("Model", model_id)
        return {"message": "Model deleted", "model_id": model_id}
    except NotFoundError:
        raise
    except ValueError as e:
        from compass.service.exceptions import ValidationError

        raise ValidationError(f"Invalid model ID format: {str(e)}")
    except (OSError, IOError, PermissionError) as e:
        # File system errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"File system error deleting model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to delete model: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    except Exception as e:
        # Unexpected errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"Unexpected error deleting model {model_id}: {e}", exc_info=True)
        raise ServiceException(
            f"Failed to delete model: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
