"""
Inference routes.
"""

from fastapi import APIRouter, HTTPException, status
import logging

from compass.service.models.model import (
    InferenceRequest,
    InferenceResponse,
    BatchInferenceRequest,
    BatchInferenceResponse,
)
from compass.service.services.inference_service import InferenceService
from compass.service.services.model_service import ModelService
from compass.service.exceptions import ServiceException, ValidationError, sanitize_error_message

router = APIRouter()
logger = logging.getLogger(__name__)

# Initialize services
model_service = ModelService()
inference_service = InferenceService(model_service)


@router.post("/api/v1/inference/predict", response_model=InferenceResponse)
async def predict(request: InferenceRequest):
    """Perform single prediction."""
    try:
        result = inference_service.predict(request)
        return result
    except ValueError as e:
        raise ValidationError(sanitize_error_message(e))
    except (FileNotFoundError, OSError, IOError) as e:
        # File-related errors
        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"File error during prediction: {e}", exc_info=True)
        raise ServiceException(
            f"Prediction failed: {error_message}",
            status_code=status.HTTP_400_BAD_REQUEST,
        )
    except (RuntimeError, MemoryError) as e:
        # Runtime or memory errors
        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"Runtime error during prediction: {e}", exc_info=True)
        raise ServiceException(
            f"Prediction failed: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    except Exception as e:
        # Unexpected errors
        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"Unexpected error during prediction: {e}", exc_info=True)
        raise ServiceException(
            f"Prediction failed: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@router.post("/api/v1/inference/batch", response_model=BatchInferenceResponse)
async def batch_predict(request: BatchInferenceRequest):
    """Perform batch predictions."""
    try:
        if not request.protein_ligand_pairs:
            raise ValidationError("Empty request list")

        # Convert pairs to requests
        requests = [
            InferenceRequest(
                protein_path=pair["protein"], ligand_path=pair["ligand"], model_id=request.model_id
            )
            for pair in request.protein_ligand_pairs
        ]

        result = inference_service.batch_predict(requests)
        return result
    except ValueError as e:
        raise ValidationError(sanitize_error_message(e))
    except (FileNotFoundError, OSError, IOError) as e:
        # File-related errors
        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"File error during batch prediction: {e}", exc_info=True)
        raise ServiceException(
            f"Batch prediction failed: {error_message}",
            status_code=status.HTTP_400_BAD_REQUEST,
        )
    except (RuntimeError, MemoryError) as e:
        # Runtime or memory errors
        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"Runtime error during batch prediction: {e}", exc_info=True)
        raise ServiceException(
            f"Batch prediction failed: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
    except Exception as e:
        # Unexpected errors
        error_message = sanitize_error_message(e, include_details=False)
        logger.error(f"Unexpected error during batch prediction: {e}", exc_info=True)
        raise ServiceException(
            f"Batch prediction failed: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@router.get("/api/v1/inference/status")
async def inference_status():
    """Get inference service status."""
    latest_model = model_service.get_latest_model()
    cache_info = inference_service.get_cache_info()
    return {
        "status": "ready" if latest_model else "no_model",
        "loaded_models": cache_info["size"],
        "max_cache_size": cache_info["max_size"],
        "cached_model_ids": cache_info["cached_models"],
        "latest_model": latest_model.model_id if latest_model else None,
        "device": str(inference_service.device),
    }
