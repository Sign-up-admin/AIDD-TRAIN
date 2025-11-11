"""
Debug configuration for Chrome debugging and vulnerability scanning.
"""

# Frontend URLs
FRONTEND_URL = "http://localhost:8501"
FRONTEND_PAGES = [
    "主页",
    "准备配体",
    "口袋预测",
    "分子对接",
    "批量口袋预测与对接",
    "预测亲和力",
    "数据管理",
    "服务监控",
    "训练管理",
]

# Backend URLs
BACKEND_URL = "http://localhost:8080"
REGISTRY_URL = "http://localhost:8500"

# Backend API Endpoints
API_ENDPOINTS = {
    # Health endpoints
    "health": [
        {"method": "GET", "path": "/health", "description": "Health check"},
        {"method": "GET", "path": "/health/ready", "description": "Readiness check"},
        {"method": "GET", "path": "/metrics", "description": "Performance metrics"},
        {"method": "GET", "path": "/backups", "description": "List backups"},
        {"method": "POST", "path": "/backups/create", "description": "Create backup"},
        {"method": "POST", "path": "/benchmark", "description": "Run benchmark"},
        {"method": "GET", "path": "/alerts", "description": "Get alerts"},
        {"method": "GET", "path": "/alerts/summary", "description": "Get alert summary"},
    ],
    # Training endpoints
    "training": [
        {"method": "POST", "path": "/api/v1/training/tasks", "description": "Create training task"},
        {"method": "GET", "path": "/api/v1/training/tasks", "description": "List training tasks"},
        {"method": "GET", "path": "/api/v1/training/tasks/{task_id}", "description": "Get training task"},
        {"method": "GET", "path": "/api/v1/training/tasks/{task_id}/progress", "description": "Get task progress"},
        {"method": "POST", "path": "/api/v1/training/tasks/{task_id}/start", "description": "Start task"},
        {"method": "POST", "path": "/api/v1/training/tasks/{task_id}/pause", "description": "Pause task"},
        {"method": "POST", "path": "/api/v1/training/tasks/{task_id}/stop", "description": "Stop task"},
        {"method": "GET", "path": "/api/v1/training/tasks/{task_id}/logs", "description": "Get task logs"},
        {"method": "GET", "path": "/api/v1/training/tasks/{task_id}/metrics", "description": "Get task metrics"},
    ],
    # Data endpoints
    "data": [
        {"method": "POST", "path": "/api/v1/data/upload", "description": "Upload dataset"},
        {"method": "GET", "path": "/api/v1/data/datasets", "description": "List datasets"},
        {"method": "GET", "path": "/api/v1/data/datasets/{dataset_id}", "description": "Get dataset"},
        {"method": "DELETE", "path": "/api/v1/data/datasets/{dataset_id}", "description": "Delete dataset"},
    ],
    # Model endpoints
    "models": [
        {"method": "GET", "path": "/api/v1/models", "description": "List models"},
        {"method": "GET", "path": "/api/v1/models/{model_id}", "description": "Get model"},
    ],
    # Inference endpoints
    "inference": [
        {"method": "POST", "path": "/api/v1/inference/predict", "description": "Single prediction"},
        {"method": "POST", "path": "/api/v1/inference/batch", "description": "Batch prediction"},
    ],
}

# WebSocket endpoints
WEBSOCKET_ENDPOINTS = [
    {"path": "/api/v1/training/tasks/{task_id}/stream", "description": "Training task stream"},
]

# Security test configurations
SECURITY_TESTS = {
    "cors": {
        "enabled": True,
        "test_origins": ["http://evil.com", "https://attacker.com", "http://localhost:3000"],
    },
    "authentication": {
        "enabled": True,
        "test_endpoints": ["/api/v1/training/tasks", "/api/v1/data/upload"],
    },
    "file_upload": {
        "enabled": True,
        "test_files": [
            {"name": "large_file.zip", "size": 11 * 1024 * 1024 * 1024},  # 11GB
            {"name": "malicious.exe", "extension": ".exe"},
            {"name": "script.js", "extension": ".js"},
            {"name": "zip_bomb.zip", "compression_ratio": 200},
        ],
    },
    "input_validation": {
        "enabled": True,
        "sql_injection": ["' OR '1'='1", "'; DROP TABLE users; --"],
        "xss": ["<script>alert('XSS')</script>", "<img src=x onerror=alert('XSS')>"],
        "path_traversal": ["../../../etc/passwd", "..\\..\\..\\windows\\system32\\config\\sam"],
    },
    "rate_limiting": {
        "enabled": True,
        "requests_per_minute": 150,  # Test rate limit (default is 100)
    },
}

# Performance test configurations
PERFORMANCE_TESTS = {
    "slow_request_threshold": 5.0,  # seconds
    "concurrent_requests": 10,
    "memory_leak_duration": 300,  # seconds
}

# Browser configuration
BROWSER_CONFIG = {
    "headless": False,  # Set to True for headless mode
    "timeout": 30,  # seconds
    "wait_timeout": 10,  # seconds
    "screenshot": True,
    "viewport": {"width": 1920, "height": 1080},
}

# Report configuration
REPORT_CONFIG = {
    "output_dir": "reports",
    "format": "markdown",
    "include_screenshots": True,
    "include_network_logs": True,
    "include_console_logs": True,
}

# Timeout and retry configuration
TIMEOUT_CONFIG = {
    "page_load": 30,
    "api_request": 10,
    "websocket_connection": 5,
    "retry_count": 3,
    "retry_delay": 1,
}





