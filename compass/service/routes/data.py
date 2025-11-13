"""
Data management routes.
"""

from fastapi import APIRouter, HTTPException, status, UploadFile, File
from typing import Optional
import os
import logging
import asyncio

from compass.service.models.dataset import DatasetCreate, DatasetResponse, DatasetListResponse
from compass.service.services.data_service import DataService
from compass.service.services.upload_queue import UploadQueueManager
from compass.service.config import SERVICE_CONFIG
from compass.service.error_codes import ErrorCode
from compass.service.exceptions import ServiceException
from compass.service.utils.file_upload_helpers import (
    validate_file_extension,
    read_and_validate_file_size,
    sanitize_dataset_metadata,
    create_temp_file,
    cleanup_temp_file,
    check_zip_bomb,
    MAX_FILE_SIZE,
)

router = APIRouter()
data_service = DataService()
logger = logging.getLogger(__name__)

# Upload queue manager for concurrency control
MAX_CONCURRENT_UPLOADS = int(os.getenv("MAX_CONCURRENT_UPLOADS", "2"))
upload_queue = UploadQueueManager(max_concurrent=MAX_CONCURRENT_UPLOADS)


@router.post(
    "/api/v1/data/upload", response_model=DatasetResponse, status_code=status.HTTP_201_CREATED
)
async def upload_dataset(
    file: UploadFile = File(...), name: Optional[str] = None, description: Optional[str] = None
):
    """
    Upload a dataset.

    Validates file size, type, and checks for zip bombs before processing.
    """
    import uuid

    # Validate file extension
    file_ext = validate_file_extension(file.filename)

    # Read and validate file size
    content_chunks, file_size = await read_and_validate_file_size(file)

    # Sanitize dataset metadata
    dataset_name, sanitized_description = sanitize_dataset_metadata(
        name, description, file.filename
    )

    # Create dataset
    dataset_id = data_service.create_dataset(name=dataset_name, description=sanitized_description)

    upload_task_id = str(uuid.uuid4())
    temp_file_path = None

    try:
        # Create temporary file
        temp_file_path = create_temp_file(content_chunks, file_ext)

        # Check for zip bomb
        if check_zip_bomb(temp_file_path):
            cleanup_temp_file(temp_file_path)
            raise ServiceException(
                "Archive file appears to be malicious (zip bomb detected)",
                error_code=ErrorCode.ZIP_BOMB_DETECTED,
            )

        # Check queue capacity
        active_count = upload_queue.get_active_count()
        if active_count >= MAX_CONCURRENT_UPLOADS:
            cleanup_temp_file(temp_file_path)
            raise ServiceException(
                f"Upload queue is full. Currently processing {active_count} uploads. Please try again later.",
                error_code=ErrorCode.UPLOAD_QUEUE_FULL,
            )

        # Track temp file for cleanup in upload function
        temp_file_for_upload = temp_file_path

        def upload_func(ds_id: str):
            """Upload function to execute in queue with guaranteed cleanup."""
            try:
                return data_service.upload_dataset(ds_id, temp_file_for_upload)
            finally:
                cleanup_temp_file(temp_file_for_upload)

        # Submit to queue (starts processing in background thread)
        upload_queue.submit_upload(upload_task_id, dataset_id, upload_func)

        # Clear reference so finally block doesn't delete it (upload_func will handle it)
        temp_file_path = None

        # Wait for upload to complete (with timeout)
        max_wait = 300  # 5 minutes
        wait_interval = 0.5  # Check every 500ms
        waited = 0

        while waited < max_wait:
            await asyncio.sleep(wait_interval)
            waited += wait_interval

            # Get updated task status
            task = upload_queue.get_task_status(upload_task_id)
            if task:
                with task.lock:
                    task_status = task.status.value
                    if task_status == "failed":
                        raise HTTPException(
                            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail=f"Upload failed: {task.error}",
                        )
                    elif task_status == "completed":
                        # Upload completed successfully
                        return data_service.get_dataset(dataset_id)
            else:
                # Task not found, might have been cleaned up
                break

        # Timeout or task not found
        # Return dataset info if available, otherwise return processing status
        dataset = data_service.get_dataset(dataset_id)
        if dataset and dataset.status == "ready":
            return dataset

        # Still processing or timeout
        return {
            "dataset_id": dataset_id,
            "status": "processing",
            "message": "Upload is being processed in background",
        }
    except HTTPException:
        raise
    except ServiceException:
        raise
    except Exception as e:
        logger.error(f"Error uploading dataset: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Failed to upload dataset"
        )
    finally:
        # Cleanup temp file if it wasn't passed to upload function
        cleanup_temp_file(temp_file_path)


@router.get("/api/v1/data/datasets", response_model=DatasetListResponse)
async def list_datasets():
    """List all datasets."""
    datasets = data_service.list_datasets()
    return DatasetListResponse(datasets=datasets, count=len(datasets))


@router.get("/api/v1/data/datasets/{dataset_id}", response_model=DatasetResponse)
async def get_dataset(dataset_id: str):
    """Get dataset by ID."""
    # Sanitize dataset_id to prevent injection attacks
    from compass.service.utils.input_sanitizer import sanitize_string

    try:
        dataset_id = sanitize_string(dataset_id, max_length=255)
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid dataset ID format: {str(e)}"
        )

    dataset = data_service.get_dataset(dataset_id)
    if not dataset:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Dataset not found")
    return dataset


@router.delete("/api/v1/data/datasets/{dataset_id}")
async def delete_dataset(dataset_id: str):
    """Delete a dataset."""
    # Sanitize dataset_id to prevent injection attacks
    from compass.service.utils.input_sanitizer import sanitize_string

    try:
        dataset_id = sanitize_string(dataset_id, max_length=255)
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid dataset ID format: {str(e)}"
        )

    success = data_service.delete_dataset(dataset_id)
    if not success:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Dataset not found")
    return {"message": "Dataset deleted", "dataset_id": dataset_id}
