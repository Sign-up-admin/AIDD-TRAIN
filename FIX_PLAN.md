# 问题修复计划

本文档提供对`AUDIT_REPORT.md`中发现的所有问题的详细修复方案。

---

## P0级别问题修复方案

### P0-1: 硬编码路径修复

**文件**: `compass/data/processing.py`

**问题代码**:
```python
log_directory = "E:/GitHub/AIDD-TRAIN/logs"
```

**修复方案**:
```python
import os
from pathlib import Path

# 使用相对路径或环境变量
log_directory = os.getenv('COMPASS_LOG_DIR', 'logs')
log_directory = Path(log_directory).resolve()
log_directory.mkdir(parents=True, exist_ok=True)
```

**测试方案**:
1. 在不同平台（Windows/Linux/Mac）测试
2. 验证日志文件正常创建
3. 验证环境变量配置生效

**风险评估**: 低风险，向后兼容

---

### P0-2: 文件上传验证修复

**文件**: `compass/service/routes/data.py`

**修复方案**:

1. **添加文件大小验证**:
```python
from fastapi import UploadFile, File, HTTPException
from compass.service.config import SERVICE_CONFIG

MAX_FILE_SIZE = SERVICE_CONFIG['upload_max_size']  # 10GB default
ALLOWED_EXTENSIONS = {'.zip', '.tar', '.tar.gz'}

@router.post("/api/v1/data/upload")
async def upload_dataset(
    file: UploadFile = File(...),
    name: Optional[str] = None,
    description: Optional[str] = None
):
    # 验证文件扩展名
    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in ALLOWED_EXTENSIONS and not file.filename.endswith('.tar.gz'):
        raise HTTPException(
            status_code=400,
            detail=f"Invalid file type. Allowed: {', '.join(ALLOWED_EXTENSIONS)}"
        )
    
    # 验证文件大小（分块读取以支持大文件）
    file_size = 0
    content_chunks = []
    while True:
        chunk = await file.read(8192)  # 8KB chunks
        if not chunk:
            break
        file_size += len(chunk)
        content_chunks.append(chunk)
        
        if file_size > MAX_FILE_SIZE:
            raise HTTPException(
                status_code=413,
                detail=f"File too large. Maximum size: {MAX_FILE_SIZE / (1024**3):.1f} GB"
            )
    
    # 重置文件指针
    await file.seek(0)
    
    # 继续原有逻辑...
```

2. **添加ZIP炸弹检测**:
```python
import zipfile

def check_zip_bomb(zip_path: Path, max_ratio: int = 100, max_files: int = 10000) -> bool:
    """Check if zip file is a zip bomb."""
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            total_size = sum(info.file_size for info in zip_ref.infolist())
            uncompressed_size = sum(info.compress_size for info in zip_ref.infolist())
            
            # Check compression ratio
            if uncompressed_size > 0:
                ratio = total_size / uncompressed_size
                if ratio > max_ratio:
                    return True
            
            # Check file count
            if len(zip_ref.namelist()) > max_files:
                return True
            
            return False
    except Exception:
        return False
```

**测试方案**:
1. 测试正常文件上传
2. 测试超大文件（超过限制）
3. 测试无效文件类型
4. 测试ZIP炸弹文件
5. 测试并发上传

**风险评估**: 中等风险，需要充分测试

---

### P0-3: 服务状态持久化修复

**文件**: `services/registry/server.py`

**修复方案**:

1. **使用SQLite作为持久化存储**:
```python
import sqlite3
from contextlib import contextmanager
from pathlib import Path

class ServiceRegistry:
    def __init__(self, db_path: str = "registry.db"):
        self.db_path = Path(db_path)
        self._init_db()
        self.services: Dict[str, ServiceInfo] = {}
        self._load_services()
    
    def _init_db(self):
        """Initialize database schema."""
        with self._get_connection() as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS services (
                    service_id TEXT PRIMARY KEY,
                    service_name TEXT NOT NULL,
                    host TEXT NOT NULL,
                    port INTEGER NOT NULL,
                    base_url TEXT NOT NULL,
                    status TEXT NOT NULL,
                    metadata TEXT,
                    version TEXT,
                    registered_at TEXT,
                    last_heartbeat TEXT
                )
            """)
    
    @contextmanager
    def _get_connection(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()
    
    def _load_services(self):
        """Load services from database."""
        with self._get_connection() as conn:
            rows = conn.execute("SELECT * FROM services").fetchall()
            for row in rows:
                service = ServiceInfo.from_dict(dict(row))
                self.services[service.service_id] = service
    
    def register_service(self, service_info: ServiceInfo):
        """Register and persist service."""
        self.services[service_info.service_id] = service_info
        with self._get_connection() as conn:
            conn.execute("""
                INSERT OR REPLACE INTO services 
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                service_info.service_id,
                service_info.service_name,
                service_info.host,
                service_info.port,
                service_info.base_url,
                service_info.status.value,
                json.dumps(service_info.metadata),
                service_info.version,
                service_info.registered_at.isoformat(),
                service_info.last_heartbeat.isoformat()
            ))
    
    def update_service_status(self, service_id: str, status: ServiceStatus):
        """Update service status in database."""
        if service_id in self.services:
            self.services[service_id].status = status
            with self._get_connection() as conn:
                conn.execute(
                    "UPDATE services SET status = ? WHERE service_id = ?",
                    (status.value, service_id)
                )
```

2. **在server.py中使用**:
```python
registry = ServiceRegistry(db_path=os.getenv('REGISTRY_DB_PATH', 'registry.db'))

@app.post("/api/v1/services/register")
async def register_service(request: ServiceRegistrationRequest):
    service_info = ServiceInfo(...)
    registry.register_service(service_info)
    return ServiceRegistrationResponse(...)
```

**测试方案**:
1. 注册服务后重启注册中心，验证服务列表恢复
2. 测试并发注册
3. 测试数据库损坏恢复

**风险评估**: 中等风险，需要数据迁移方案

---

### P0-4: 线程资源清理修复

**文件**: `compass/service/services/training_service.py`

**修复方案**:

```python
def start_task(self, task_id: str) -> bool:
    """Start a training task."""
    with self.lock:
        # ... existing validation ...
        
        # Start training in background thread
        def run_training():
            try:
                self._run_training(task_id)
            except Exception as e:
                logger.error(f"Training task {task_id} failed: {e}", exc_info=True)
                with self.lock:
                    if task_id in self.tasks:
                        self.tasks[task_id].status = TaskStatus.FAILED
                        self.tasks[task_id].error = str(e)
                        self.tasks[task_id].updated_at = datetime.now()
                        self._log(task_id, f"ERROR: {e}")
            finally:
                # Clean up thread reference
                with self.lock:
                    if task_id in self.task_threads:
                        del self.task_threads[task_id]
                    if task_id in self.progress_trackers:
                        del self.progress_trackers[task_id]
        
        thread = threading.Thread(target=run_training, daemon=True)
        thread.start()
        self.task_threads[task_id] = thread
        
        logger.info(f"Started training task: {task_id}")
        return True

def stop_task(self, task_id: str) -> bool:
    """Stop a training task."""
    with self.lock:
        if task_id not in self.tasks:
            return False
        
        task = self.tasks[task_id]
        if task.status != TaskStatus.RUNNING:
            return False
        
        task.status = TaskStatus.CANCELLED
        task.updated_at = datetime.now()
        
        # Clean up thread reference
        if task_id in self.task_threads:
            thread = self.task_threads[task_id]
            del self.task_threads[task_id]
            # Note: daemon threads will be cleaned up automatically
            # but we remove the reference to prevent memory leak
    
    self._log(task_id, "Training task cancelled")
    logger.info(f"Stopped training task: {task_id}")
    return True
```

**测试方案**:
1. 创建多个任务，验证线程清理
2. 长时间运行，监控内存使用
3. 测试任务完成后的清理

**风险评估**: 低风险

---

### P0-5: 统一异常处理修复

**文件**: `compass/service/server.py`

**修复方案**:

1. **创建异常处理器**:
```python
from fastapi import FastAPI, Request, status
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
import logging

logger = logging.getLogger(__name__)

class ServiceException(Exception):
    """Base exception for service errors."""
    def __init__(self, message: str, status_code: int = 500, detail: Optional[Dict] = None):
        self.message = message
        self.status_code = status_code
        self.detail = detail or {}

@app.exception_handler(ServiceException)
async def service_exception_handler(request: Request, exc: ServiceException):
    """Handle service exceptions."""
    logger.error(f"Service error: {exc.message}", exc_info=True)
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "error": exc.message,
            "detail": exc.detail,
            "path": str(request.url)
        }
    )

@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    """Handle unexpected exceptions."""
    logger.error(f"Unexpected error: {exc}", exc_info=True)
    # Don't expose internal error details to users
    return JSONResponse(
        status_code=500,
        content={
            "error": "Internal server error",
            "path": str(request.url)
        }
    )

@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Handle validation errors."""
    return JSONResponse(
        status_code=422,
        content={
            "error": "Validation error",
            "detail": exc.errors()
        }
    )
```

2. **在路由中使用**:
```python
from compass.service.exceptions import ServiceException

@router.post("/api/v1/training/tasks")
async def create_task(request: TrainingTaskCreate):
    """Create a new training task."""
    try:
        task_id = training_service.create_task(
            config=request.config,
            dataset_id=request.dataset_id,
            description=request.description
        )
        return training_service.get_task(task_id)
    except ValueError as e:
        raise ServiceException(str(e), status_code=400)
    except Exception as e:
        logger.error(f"Failed to create task: {e}", exc_info=True)
        raise ServiceException("Failed to create task", status_code=500)
```

**测试方案**:
1. 测试各种异常场景
2. 验证错误响应格式
3. 验证日志记录

**风险评估**: 低风险

---

## P1级别问题修复方案（部分）

### P1-1: 统一日志系统

**创建统一日志配置模块**: `compass/service/logging_config.py`

```python
import logging
import logging.handlers
import os
from pathlib import Path
from typing import Optional

def setup_logging(
    log_dir: str = "logs",
    log_level: str = "INFO",
    service_name: str = "compass",
    enable_file: bool = True,
    enable_console: bool = True
):
    """
    Setup unified logging configuration.
    
    Args:
        log_dir: Directory for log files
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        service_name: Name of the service
        enable_file: Enable file logging
        enable_console: Enable console logging
    """
    log_path = Path(log_dir)
    log_path.mkdir(parents=True, exist_ok=True)
    
    # Create formatters
    detailed_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    simple_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )
    
    # Get root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, log_level.upper()))
    
    # Clear existing handlers
    root_logger.handlers.clear()
    
    # File handler with rotation
    if enable_file:
        file_handler = logging.handlers.RotatingFileHandler(
            log_path / f"{service_name}.log",
            maxBytes=10 * 1024 * 1024,  # 10MB
            backupCount=5,
            encoding='utf-8'
        )
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(detailed_formatter)
        root_logger.addHandler(file_handler)
        
        # Error file handler
        error_handler = logging.handlers.RotatingFileHandler(
            log_path / f"{service_name}_errors.log",
            maxBytes=10 * 1024 * 1024,
            backupCount=5,
            encoding='utf-8'
        )
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(detailed_formatter)
        root_logger.addHandler(error_handler)
    
    # Console handler
    if enable_console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(simple_formatter)
        root_logger.addHandler(console_handler)
    
    return root_logger
```

**使用方式**:
```python
# In server.py
from compass.service.logging_config import setup_logging

setup_logging(
    log_dir=SERVICE_CONFIG['log_dir'],
    log_level=os.getenv('LOG_LEVEL', 'INFO'),
    service_name='compass-service'
)
```

---

### P1-2: 错误信息安全修复

**修复方案**:

```python
def sanitize_error_message(error: Exception, include_details: bool = False) -> str:
    """
    Sanitize error message to avoid exposing internal details.
    
    Args:
        error: The exception
        include_details: Whether to include detailed information (for logging)
    
    Returns:
        Sanitized error message
    """
    error_type = type(error).__name__
    error_message = str(error)
    
    # User-friendly messages for common errors
    user_messages = {
        'ValueError': 'Invalid input provided',
        'FileNotFoundError': 'Required file not found',
        'PermissionError': 'Permission denied',
        'ConnectionError': 'Failed to connect to service',
        'TimeoutError': 'Request timed out',
    }
    
    if include_details:
        # For logging, include full details
        return f"{error_type}: {error_message}"
    
    # For user response, use friendly message
    user_message = user_messages.get(error_type, 'An error occurred')
    
    # Log full error details
    logger.error(f"Error: {error_type}: {error_message}", exc_info=True)
    
    return user_message
```

**使用方式**:
```python
@router.post("/api/v1/inference/predict")
async def predict(request: InferenceRequest):
    """Perform single prediction."""
    try:
        result = inference_service.predict(request)
        return result
    except ValueError as e:
        raise HTTPException(
            status_code=400,
            detail=sanitize_error_message(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=sanitize_error_message(e)
        )
```

---

### P1-3: 批量推理错误处理修复

**修复方案**:

```python
def batch_predict(self, requests: List[InferenceRequest]) -> BatchInferenceResponse:
    """Perform batch predictions."""
    if not requests:
        raise ValueError("Empty request list")
    
    # Pre-load model before processing
    model_id = requests[0].model_id if requests else None
    try:
        model = self.load_model(model_id)
    except Exception as e:
        logger.error(f"Failed to load model for batch prediction: {e}")
        raise ServiceException(
            "Failed to load model for batch prediction",
            status_code=500
        )
    
    results = []
    success_count = 0
    failed_count = 0
    
    for idx, request in enumerate(requests):
        try:
            # Use pre-loaded model
            result = self._predict_with_model(model, request)
            results.append({
                'protein_path': request.protein_path,
                'ligand_path': request.ligand_path,
                'binding_affinity': result.binding_affinity,
                'success': True
            })
            success_count += 1
        except Exception as e:
            logger.warning(f"Failed to process request {idx}: {e}")
            results.append({
                'protein_path': request.protein_path,
                'ligand_path': request.ligand_path,
                'error': sanitize_error_message(e),
                'success': False
            })
            failed_count += 1
    
    return BatchInferenceResponse(
        results=results,
        total_count=len(requests),
        success_count=success_count,
        failed_count=failed_count
    )
```

---

## 修复实施检查清单

### 第一阶段检查清单
- [ ] P0-1: 修复硬编码路径
- [ ] P0-2: 添加文件上传验证
- [ ] P0-5: 添加统一异常处理
- [ ] P1-1: 统一日志系统
- [ ] P1-2: 修复错误信息泄露

### 第二阶段检查清单
- [ ] P0-3: 实现服务状态持久化
- [ ] P0-4: 修复线程资源清理
- [ ] P1-3: 修复批量推理错误处理
- [ ] P1-4: 添加服务清理机制
- [ ] P1-5: 添加并发控制
- [ ] P1-6: 实现任务取消机制
- [ ] P1-7: 添加模型缓存管理
- [ ] P1-8: 添加注册重试机制
- [ ] P1-9: 修复负载均衡计数
- [ ] P1-10: 添加参数验证
- [ ] P1-11: 改进临时文件清理
- [ ] P1-12: 修复心跳线程清理

### 第三阶段检查清单
- [ ] P2-1: 添加类型注解
- [ ] P2-2: 消除代码重复
- [ ] P2-3: 完善API文档
- [ ] P2-4: 统一配置管理
- [ ] P2-5: 添加单元测试
- [ ] P2-6: 规范日志级别
- [ ] P2-7: 添加性能监控
- [ ] P2-8: 添加请求限流
- [ ] P2-9: 添加验证中间件
- [ ] P2-10: 实现版本控制

---

## 测试策略

### 单元测试
- 为每个修复的功能编写单元测试
- 目标覆盖率：60%+

### 集成测试
- 测试服务间交互
- 测试错误场景
- 测试并发场景

### 回归测试
- 验证修复不影响现有功能
- 运行完整测试套件

---

## 风险评估

### 高风险修复
- P0-3: 服务状态持久化（需要数据迁移）
- P0-2: 文件上传验证（可能影响现有功能）

### 中风险修复
- P1-1: 日志系统统一（需要测试所有日志输出）
- P1-6: 任务取消机制（需要修改训练循环）

### 低风险修复
- P0-1: 硬编码路径
- P0-4: 线程清理
- P0-5: 异常处理

---

**修复计划版本**: 1.0  
**最后更新**: 2025-01-XX


