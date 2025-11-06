"""
Data management routes.
"""

from fastapi import APIRouter, HTTPException, status, UploadFile, File
from typing import Optional
import tempfile
import os
import zipfile
import tarfile
from pathlib import Path
import logging
import asyncio

from compass.service.models.dataset import DatasetCreate, DatasetResponse, DatasetListResponse
from compass.service.services.data_service import DataService
from compass.service.services.upload_queue import UploadQueueManager
from compass.service.config import SERVICE_CONFIG
from compass.service.error_codes import ErrorCode

router = APIRouter()
data_service = DataService()
logger = logging.getLogger(__name__)

# Upload queue manager for concurrency control
MAX_CONCURRENT_UPLOADS = int(os.getenv("MAX_CONCURRENT_UPLOADS", "2"))
upload_queue = UploadQueueManager(max_concurrent=MAX_CONCURRENT_UPLOADS)

# File upload configuration
MAX_FILE_SIZE = SERVICE_CONFIG.get("upload_max_size", 10 * 1024 * 1024 * 1024)  # Default 10GB
ALLOWED_EXTENSIONS = {".zip", ".tar", ".tar.gz"}
MAX_COMPRESSION_RATIO = 100  # Maximum compression ratio to detect zip bombs
MAX_FILES_IN_ARCHIVE = 10000  # Maximum files in archive


def check_zip_bomb(
    file_path: str, max_ratio: int = MAX_COMPRESSION_RATIO, max_files: int = MAX_FILES_IN_ARCHIVE
) -> bool:
    """
    Check if zip file is a zip bomb.

    Args:
        file_path: Path to the zip file
        max_ratio: Maximum compression ratio allowed
        max_files: Maximum number of files allowed in archive

    Returns:
        bool: True if zip bomb detected, False otherwise
    """
    try:
        file_ext = Path(file_path).suffix.lower()

        if file_ext == ".zip":
            with zipfile.ZipFile(file_path, "r") as zip_ref:
                file_count = len(zip_ref.namelist())
                if file_count > max_files:
                    logger.warning(f"Zip bomb detected: too many files ({file_count})")
                    return True

                total_size = 0
                compressed_size = 0
                for info in zip_ref.infolist():
                    total_size += info.file_size
                    compressed_size += info.compress_size

                if compressed_size > 0:
                    ratio = total_size / compressed_size
                    if ratio > max_ratio:
                        logger.warning(f"Zip bomb detected: high compression ratio ({ratio:.1f})")
                        return True
        elif file_ext in [".tar", ".gz"] or file_path.endswith(".tar.gz"):
            with tarfile.open(file_path, "r:*") as tar_ref:
                file_count = len(tar_ref.getnames())
                if file_count > max_files:
                    logger.warning(f"Tar bomb detected: too many files ({file_count})")
                    return True

        return False
    except Exception as e:
        logger.error(f"Error checking for zip bomb: {e}")
        # If we can't check, assume it's safe but log the error
        return False


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
    # Validate file extension
    if not file.filename:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Filename is required")

    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in ALLOWED_EXTENSIONS and not file.filename.endswith(".tar.gz"):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid file type. Allowed types: {', '.join(ALLOWED_EXTENSIONS)}",
        )

    # Validate and read file with size check
    file_size = 0
    content_chunks = []
    chunk_size = 8192  # 8KB chunks

    try:
        while True:
            chunk = await file.read(chunk_size)
            if not chunk:
                break
            file_size += len(chunk)

            if file_size > MAX_FILE_SIZE:
                raise HTTPException(
                    status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
                    detail=f"File too large. Maximum size: {MAX_FILE_SIZE / (1024**3):.1f} GB",
                )

            content_chunks.append(chunk)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error reading uploaded file: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Failed to read uploaded file"
        )

    # Reset file pointer for later use if needed
    await file.seek(0)

    # Create dataset
    dataset_name = name or file.filename or "uploaded_dataset"
    dataset_id = data_service.create_dataset(name=dataset_name, description=description)

    # Save uploaded file temporarily
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=file_ext) as tmp_file:
            for chunk in content_chunks:
                tmp_file.write(chunk)
            tmp_path = tmp_file.name

        # Check for zip bomb
        if check_zip_bomb(tmp_path):
            from compass.service.exceptions import ServiceException

            raise ServiceException(
                "Archive file appears to be malicious (zip bomb detected)",
                error_code=ErrorCode.ZIP_BOMB_DETECTED,
            )

        # Check queue capacity
        active_count = upload_queue.get_active_count()
        queue_size = upload_queue.get_queue_size()

        if active_count >= MAX_CONCURRENT_UPLOADS:
            # Queue is full, reject request
            from compass.service.exceptions import ServiceException

            raise ServiceException(
                f"Upload queue is full. Currently processing {active_count} uploads. Please try again later.",
                error_code=ErrorCode.UPLOAD_QUEUE_FULL,
            )

        # Submit to upload queue
        import uuid

        upload_task_id = str(uuid.uuid4())

        # Track temp file for cleanup
        temp_file_path = tmp_path

        def upload_func(ds_id: str):
            """Upload function to execute in queue."""
            try:
                return data_service.upload_dataset(ds_id, temp_file_path)
            finally:
                # Cleanup temp file after upload completes
                if temp_file_path and os.path.exists(temp_file_path):
                    try:
                        os.remove(temp_file_path)
                        logger.debug(f"Cleaned up temp file after upload: {temp_file_path}")
                    except Exception as e:
                        logger.warning(f"Failed to cleanup temp file {temp_file_path}: {e}")

        # Submit to queue (starts processing in background thread)
        upload_task = upload_queue.submit_upload(upload_task_id, dataset_id, upload_func)

        # Don't cleanup here - let the upload function handle it
        tmp_path = None  # Prevent finally block from deleting it

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
                    status = task.status.value
                    if status == "failed":
                        raise HTTPException(
                            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail=f"Upload failed: {task.error}",
                        )
                    elif status == "completed":
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
    except Exception as e:
        logger.error(f"Error uploading dataset: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Failed to upload dataset"
        )
    finally:
        # Cleanup temp file
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except Exception as e:
                logger.warning(f"Failed to cleanup temp file {tmp_path}: {e}")


@router.get("/api/v1/data/datasets", response_model=DatasetListResponse)
async def list_datasets():
    """List all datasets."""
    datasets = data_service.list_datasets()
    return DatasetListResponse(datasets=datasets, count=len(datasets))


@router.get("/api/v1/data/datasets/{dataset_id}", response_model=DatasetResponse)
async def get_dataset(dataset_id: str):
    """Get dataset by ID."""
    dataset = data_service.get_dataset(dataset_id)
    if not dataset:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Dataset {dataset_id} not found"
        )
    return dataset


@router.delete("/api/v1/data/datasets/{dataset_id}")
async def delete_dataset(dataset_id: str):
    """Delete a dataset."""
    success = data_service.delete_dataset(dataset_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Dataset {dataset_id} not found"
        )
    return {"message": f"Dataset {dataset_id} deleted", "dataset_id": dataset_id}
