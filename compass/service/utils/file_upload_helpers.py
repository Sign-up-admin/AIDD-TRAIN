"""
Helper functions for file upload processing.

This module provides utility functions for validating and processing uploaded files.
"""

import os
import tempfile
import logging
from pathlib import Path
from typing import Tuple, List, Optional
from fastapi import UploadFile, HTTPException, status

from compass.service.error_codes import ErrorCode
from compass.service.exceptions import ServiceException
from compass.service.utils.input_sanitizer import (
    sanitize_string,
    sanitize_description,
    sanitize_filename,
)

logger = logging.getLogger(__name__)

# File upload configuration
MAX_FILE_SIZE = 10 * 1024 * 1024 * 1024  # Default 10GB
ALLOWED_EXTENSIONS = {".zip", ".tar", ".tar.gz"}
MAX_COMPRESSION_RATIO = 100  # Maximum compression ratio to detect zip bombs
MAX_FILES_IN_ARCHIVE = 10000  # Maximum files in archive


def validate_file_extension(filename: Optional[str]) -> str:
    """
    Validate file extension.

    Args:
        filename: Name of the uploaded file

    Returns:
        File extension (lowercase)

    Raises:
        HTTPException: If filename is missing or extension is invalid
    """
    if not filename:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Filename is required")

    file_ext = Path(filename).suffix.lower()
    if file_ext not in ALLOWED_EXTENSIONS and not filename.endswith(".tar.gz"):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid file type. Allowed types: {', '.join(ALLOWED_EXTENSIONS)}",
        )

    return file_ext


async def read_and_validate_file_size(
    file: UploadFile, max_size: int = MAX_FILE_SIZE, chunk_size: int = 8192
) -> Tuple[List[bytes], int]:
    """
    Read file content and validate size.

    Args:
        file: Uploaded file
        max_size: Maximum file size in bytes
        chunk_size: Chunk size for reading

    Returns:
        Tuple of (content_chunks, file_size)

    Raises:
        HTTPException: If file is too large or read fails
    """
    file_size = 0
    content_chunks = []

    try:
        while True:
            chunk = await file.read(chunk_size)
            if not chunk:
                break
            file_size += len(chunk)

            if file_size > max_size:
                raise HTTPException(
                    status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
                    detail=f"File too large. Maximum size: {max_size / (1024**3):.1f} GB",
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

    return content_chunks, file_size


def sanitize_dataset_metadata(
    name: Optional[str], description: Optional[str], filename: Optional[str]
) -> Tuple[str, Optional[str]]:
    """
    Sanitize dataset name and description.

    Args:
        name: Dataset name (optional)
        description: Dataset description (optional)
        filename: Original filename (optional)

    Returns:
        Tuple of (sanitized_name, sanitized_description)
    """
    # Sanitize dataset name
    if name:
        dataset_name = sanitize_string(name, max_length=255)
    elif filename:
        dataset_name = sanitize_filename(filename)
        if not dataset_name or len(dataset_name) == 0:
            dataset_name = "uploaded_dataset"
    else:
        dataset_name = "uploaded_dataset"

    # Sanitize description
    sanitized_description = sanitize_description(description) if description else None

    return dataset_name, sanitized_description


def create_temp_file(content_chunks: List[bytes], file_ext: str) -> str:
    """
    Create temporary file from content chunks.

    Args:
        content_chunks: List of file content chunks
        file_ext: File extension

    Returns:
        Path to temporary file

    Raises:
        OSError: If file creation fails
    """
    with tempfile.NamedTemporaryFile(delete=False, suffix=file_ext) as tmp_file:
        for chunk in content_chunks:
            tmp_file.write(chunk)
        return tmp_file.name


def cleanup_temp_file(file_path: Optional[str]) -> bool:
    """
    Clean up temporary file.

    Args:
        file_path: Path to temporary file

    Returns:
        True if cleanup succeeded, False otherwise
    """
    if file_path and os.path.exists(file_path):
        try:
            os.remove(file_path)
            logger.debug(f"Cleaned up temp file: {file_path}")
            return True
        except (OSError, PermissionError) as e:
            logger.warning(f"Failed to cleanup temp file {file_path}: {e}")
            return False
    return True


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
    import zipfile
    import tarfile

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
