"""
Data models for service registry.
"""

from pydantic import BaseModel, Field
from typing import Dict, Optional, List
from datetime import datetime


class ServiceRegistrationRequest(BaseModel):
    """Request model for service registration."""

    service_name: str
    host: str
    port: int
    metadata: Dict = Field(default_factory=dict)
    version: str = "1.0.0"


class ServiceRegistrationResponse(BaseModel):
    """Response model for service registration."""

    service_id: str
    message: str
    registered_at: datetime


class ServiceQueryResponse(BaseModel):
    """Response model for service query."""

    services: List[Dict]
    count: int


class HealthCheckResponse(BaseModel):
    """Response model for health check."""

    status: str
    message: str
    timestamp: datetime = Field(default_factory=datetime.now)
