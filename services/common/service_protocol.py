"""
Service protocol definitions for service registry.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, List
from datetime import datetime
from enum import Enum


class ServiceStatus(Enum):
    """Service health status."""

    HEALTHY = "healthy"
    UNHEALTHY = "unhealthy"
    UNKNOWN = "unknown"


@dataclass
class ServiceInfo:
    """
    Service information model for registry.
    """

    service_name: str
    service_id: str
    host: str
    port: int
    base_url: str
    status: ServiceStatus = ServiceStatus.UNKNOWN
    metadata: Dict = field(default_factory=dict)
    registered_at: Optional[datetime] = None
    last_heartbeat: Optional[datetime] = None
    version: str = "1.0.0"

    def __post_init__(self):
        if self.registered_at is None:
            self.registered_at = datetime.now()
        if self.last_heartbeat is None:
            self.last_heartbeat = datetime.now()

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "service_name": self.service_name,
            "service_id": self.service_id,
            "host": self.host,
            "port": self.port,
            "base_url": self.base_url,
            "status": self.status.value,
            "metadata": self.metadata,
            "registered_at": self.registered_at.isoformat() if self.registered_at else None,
            "last_heartbeat": self.last_heartbeat.isoformat() if self.last_heartbeat else None,
            "version": self.version,
        }

    @classmethod
    def from_dict(cls, data: Dict) -> "ServiceInfo":
        """Create from dictionary."""
        info = cls(
            service_name=data["service_name"],
            service_id=data["service_id"],
            host=data["host"],
            port=data["port"],
            base_url=data["base_url"],
            status=ServiceStatus(data.get("status", "unknown")),
            metadata=data.get("metadata", {}),
            version=data.get("version", "1.0.0"),
        )
        # Parse datetime strings - handle both string and datetime objects
        if "registered_at" in data and data["registered_at"]:
            if isinstance(data["registered_at"], str):
                info.registered_at = datetime.fromisoformat(data["registered_at"])
            elif isinstance(data["registered_at"], datetime):
                info.registered_at = data["registered_at"]
        if "last_heartbeat" in data and data["last_heartbeat"]:
            if isinstance(data["last_heartbeat"], str):
                info.last_heartbeat = datetime.fromisoformat(data["last_heartbeat"])
            elif isinstance(data["last_heartbeat"], datetime):
                info.last_heartbeat = data["last_heartbeat"]
        return info
