"""
Dataset models for data service.
"""

from pydantic import BaseModel, Field
from typing import Dict, Optional, List
from datetime import datetime


class DatasetCreate(BaseModel):
    """Request model for creating a dataset."""

    name: str
    description: Optional[str] = None
    metadata: Dict = Field(default_factory=dict)


class DatasetResponse(BaseModel):
    """Response model for dataset."""

    dataset_id: str
    name: str
    description: Optional[str] = None
    size: int = 0
    file_count: int = 0
    status: str = "pending"  # pending, processing, ready, error
    created_at: datetime
    updated_at: datetime
    metadata: Dict = Field(default_factory=dict)
    error: Optional[str] = None


class DatasetListResponse(BaseModel):
    """Response model for dataset list."""

    datasets: List[DatasetResponse]
    count: int
