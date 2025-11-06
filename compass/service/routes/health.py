"""
Health check routes.
"""
from fastapi import APIRouter, status, Query
from fastapi.responses import JSONResponse
from datetime import datetime
from typing import Dict, Optional
import logging

router = APIRouter(tags=["health"])
logger = logging.getLogger(__name__)


@router.get(
    "/health",
    summary="Health Check",
    description="Basic health check endpoint to verify service is running",
    response_description="Service health status",
    status_code=status.HTTP_200_OK
)
async def health() -> Dict[str, str]:
    """
    Health check endpoint.
    
    Returns:
        Dict containing:
        - status: "healthy" if service is running
        - timestamp: Current ISO timestamp
    
    Example response:
        {
            "status": "healthy",
            "timestamp": "2025-01-15T10:30:00.123456"
        }
    """
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat()
    }


@router.get(
    "/health/ready",
    summary="Readiness Check",
    description="Readiness check endpoint to verify service is ready to accept requests",
    response_description="Service readiness status",
    status_code=status.HTTP_200_OK
)
async def ready() -> Dict[str, str]:
    """
    Readiness check endpoint.
    
    This endpoint checks if the service is ready to accept requests.
    It verifies that all required components are initialized.
    
    Returns:
        Dict containing:
        - status: "ready" if service is ready
        - timestamp: Current ISO timestamp
    
    Example response:
        {
            "status": "ready",
            "timestamp": "2025-01-15T10:30:00.123456"
        }
    """
    return {
        "status": "ready",
        "timestamp": datetime.now().isoformat()
    }


@router.get("/metrics")
async def metrics():
    """
    Get performance metrics.
    
    Returns:
        Dict containing performance metrics including:
        - total_requests: Total number of requests
        - total_errors: Total number of errors
        - error_rate: Error rate (0-1)
        - response_time: Statistics (avg, min, max, p50, p95, p99)
        - status_codes: Count by status code
        - endpoints: Per-endpoint statistics
    """
    from compass.service.middleware.metrics import get_metrics_collector
    collector = get_metrics_collector()
    return collector.get_metrics()


@router.get("/backups", summary="List Backups", description="List all available backups")
async def list_backups():
    """
    List all available backups.
    
    Returns:
        List of backup information
    """
    from compass.service.utils.backup import get_backup_manager
    backup_manager = get_backup_manager()
    if not backup_manager:
        return {"backups": [], "message": "Backup disabled"}
    return {"backups": backup_manager.list_backups()}


@router.post("/backups/create", summary="Create Backup", description="Create a manual backup")
async def create_backup():
    """
    Create a manual backup.
    
    Returns:
        Backup creation result
    """
    from compass.service.utils.backup import get_backup_manager
    from compass.service.config import SERVICE_CONFIG
    import os
    
    backup_manager = get_backup_manager()
    if not backup_manager:
        from fastapi import HTTPException
        raise HTTPException(status_code=503, detail="Backup disabled")
    
    backup_path = backup_manager.create_backup(
        datasets_dir=SERVICE_CONFIG.get('data_dir'),
        models_dir=SERVICE_CONFIG.get('checkpoint_dir'),
        database_file=os.getenv('REGISTRY_DB_PATH', 'registry.db')
    )
    
    if backup_path:
        return {"status": "success", "backup_path": backup_path}
    else:
        from fastapi import HTTPException
        raise HTTPException(status_code=500, detail="Failed to create backup")


@router.post("/benchmark", summary="Run Benchmark", description="Run performance benchmark tests")
async def run_benchmark(
    benchmark_type: str = "inference",
    iterations: int = 10
):
    """
    Run performance benchmark.
    
    Args:
        benchmark_type: Type of benchmark ('inference', 'all')
        iterations: Number of iterations
    
    Returns:
        Benchmark results
    """
    from compass.service.utils.benchmark import BenchmarkRunner, benchmark_inference_service
    from compass.service.services.inference_service import InferenceService
    from compass.service.services.model_service import ModelService
    
    runner = BenchmarkRunner()
    
    if benchmark_type == "inference" or benchmark_type == "all":
        try:
            model_service = ModelService()
            inference_service = InferenceService(model_service)
            result = benchmark_inference_service(
                inference_service,
                num_iterations=iterations
            )
            runner.results.append(result)
        except Exception as e:
            logger.error(f"Benchmark failed: {e}", exc_info=True)
            from fastapi import HTTPException
            raise HTTPException(status_code=500, detail=f"Benchmark failed: {str(e)}")
    
    return {
        "benchmark_type": benchmark_type,
        "iterations": iterations,
        "results": runner.get_results()
    }


@router.get("/alerts", summary="Get Alerts", description="Get monitoring alerts")
async def get_alerts(
    level: Optional[str] = None,
    limit: int = 100
):
    """
    Get monitoring alerts.
    
    Args:
        level: Filter by alert level (info, warning, error, critical)
        limit: Maximum number of alerts to return
    
    Returns:
        List of alerts
    """
    from compass.service.monitoring.alert_manager import get_alert_manager, AlertLevel
    
    alert_manager = get_alert_manager()
    
    alert_level = None
    if level:
        try:
            alert_level = AlertLevel(level.lower())
        except ValueError:
            from fastapi import HTTPException
            raise HTTPException(status_code=400, detail=f"Invalid alert level: {level}")
    
    return {
        "alerts": alert_manager.get_alerts(level=alert_level, limit=limit)
    }


@router.get("/alerts/summary", summary="Get Alert Summary", description="Get alert summary statistics")
async def get_alert_summary():
    """
    Get alert summary statistics.
    
    Returns:
        Alert summary including counts by level and recent alerts
    """
    from compass.service.monitoring.alert_manager import get_alert_manager
    
    alert_manager = get_alert_manager()
    return alert_manager.get_alert_summary()

