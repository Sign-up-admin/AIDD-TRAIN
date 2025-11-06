"""
Unified configuration management for COMPASS service.
"""
import os
import json
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass, field, asdict
from enum import Enum

# Optional YAML support
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


class ConfigSource(str, Enum):
    """Configuration source priority."""
    ENV = "environment"
    FILE = "file"
    DEFAULT = "default"


@dataclass
class ServiceConfig:
    """Service configuration dataclass."""
    # Service settings
    service_name: str = "compass-service"
    service_id: Optional[str] = None
    host: str = "0.0.0.0"
    port: int = 8080
    
    # Registry settings
    registry_url: str = "http://localhost:8500"
    registry_retry_max: int = 5
    check_interval: int = 30
    heartbeat_interval: int = 30
    
    # Resource settings
    max_workers: int = 4
    data_dir: str = "./data"
    checkpoint_dir: str = "./checkpoints"
    log_dir: str = "./logs"
    upload_max_size: int = 10 * 1024 * 1024 * 1024  # 10GB
    
    # Logging settings
    log_level: str = "INFO"
    log_max_bytes: int = 10 * 1024 * 1024  # 10MB
    log_backup_count: int = 5
    log_rotation_strategy: str = "size"
    
    # Rate limiting
    rate_limit_default: int = 100
    rate_limit_window: int = 60
    rate_limit_training: int = 5
    rate_limit_upload: int = 10
    rate_limit_inference: int = 50
    
    # Model cache
    model_cache_size: int = 3
    
    # Database
    registry_db_path: str = "registry.db"
    
    # Upload queue
    max_concurrent_uploads: int = 2


class ConfigManager:
    """Unified configuration manager."""
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration manager.
        
        Args:
            config_file: Optional path to configuration file (YAML or JSON)
        """
        self.config_file = config_file
        self.config: ServiceConfig = ServiceConfig()
        self.config_sources: Dict[str, ConfigSource] = {}
        self.load_config()
    
    def load_config(self):
        """Load configuration from all sources (file, env, defaults)."""
        # Load from file first
        if self.config_file and Path(self.config_file).exists():
            self._load_from_file(self.config_file)
        
        # Override with environment variables
        self._load_from_env()
    
    def _load_from_file(self, file_path: str):
        """Load configuration from file."""
        path = Path(file_path)
        file_config = {}
        
        try:
            if path.suffix.lower() in ['.yaml', '.yml']:
                if not YAML_AVAILABLE:
                    raise ImportError("YAML support not available. Install PyYAML: pip install pyyaml")
                with open(path, 'r', encoding='utf-8') as f:
                    file_config = yaml.safe_load(f) or {}
            elif path.suffix.lower() == '.json':
                with open(path, 'r', encoding='utf-8') as f:
                    file_config = json.load(f)
            else:
                raise ValueError(f"Unsupported config file format: {path.suffix}")
            
            # Update config from file
            for key, value in file_config.items():
                if hasattr(self.config, key):
                    setattr(self.config, key, value)
                    self.config_sources[key] = ConfigSource.FILE
            
        except Exception as e:
            print(f"Warning: Failed to load config file {file_path}: {e}")
    
    def _load_from_env(self):
        """Load configuration from environment variables."""
        env_mapping = {
            # Service settings
            'COMPASS_SERVICE_NAME': 'service_name',
            'COMPASS_HOST': 'host',
            'COMPASS_PORT': 'port',
            'COMPASS_SERVICE_ID': 'service_id',
            
            # Registry settings
            'REGISTRY_URL': 'registry_url',
            'REGISTRY_RETRY_MAX': 'registry_retry_max',
            'COMPASS_CHECK_INTERVAL': 'check_interval',
            'COMPASS_HEARTBEAT_INTERVAL': 'heartbeat_interval',
            
            # Resource settings
            'COMPASS_MAX_WORKERS': 'max_workers',
            'COMPASS_DATA_DIR': 'data_dir',
            'COMPASS_CHECKPOINT_DIR': 'checkpoint_dir',
            'COMPASS_LOG_DIR': 'log_dir',
            'COMPASS_UPLOAD_MAX_SIZE': 'upload_max_size',
            
            # Logging settings
            'LOG_LEVEL': 'log_level',
            'LOG_MAX_BYTES': 'log_max_bytes',
            'LOG_BACKUP_COUNT': 'log_backup_count',
            'LOG_ROTATION_STRATEGY': 'log_rotation_strategy',
            
            # Rate limiting
            'RATE_LIMIT_DEFAULT': 'rate_limit_default',
            'RATE_LIMIT_WINDOW': 'rate_limit_window',
            'RATE_LIMIT_TRAINING': 'rate_limit_training',
            'RATE_LIMIT_UPLOAD': 'rate_limit_upload',
            'RATE_LIMIT_INFERENCE': 'rate_limit_inference',
            
            # Model cache
            'MODEL_CACHE_SIZE': 'model_cache_size',
            
            # Database
            'REGISTRY_DB_PATH': 'registry_db_path',
            
            # Upload queue
            'MAX_CONCURRENT_UPLOADS': 'max_concurrent_uploads',
        }
        
        for env_var, attr_name in env_mapping.items():
            value = os.getenv(env_var)
            if value is not None:
                # Type conversion
                attr = getattr(self.config, attr_name)
                if isinstance(attr, int):
                    value = int(value)
                elif isinstance(attr, bool):
                    value = value.lower() in ('true', '1', 'yes', 'on')
                
                setattr(self.config, attr_name, value)
                self.config_sources[attr_name] = ConfigSource.ENV
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value.
        
        Args:
            key: Configuration key
            default: Default value if key not found
        
        Returns:
            Configuration value
        """
        return getattr(self.config, key, default)
    
    def get_dict(self) -> Dict[str, Any]:
        """Get configuration as dictionary."""
        return asdict(self.config)
    
    def get_source(self, key: str) -> ConfigSource:
        """Get configuration source for a key."""
        return self.config_sources.get(key, ConfigSource.DEFAULT)
    
    def save_to_file(self, file_path: str, format: str = "yaml"):
        """
        Save current configuration to file.
        
        Args:
            file_path: Path to save configuration
            format: File format ('yaml' or 'json')
        """
        path = Path(file_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        config_dict = self.get_dict()
        
        try:
            if format.lower() == 'yaml':
                if not YAML_AVAILABLE:
                    raise ImportError("YAML support not available. Install PyYAML: pip install pyyaml")
                with open(path, 'w', encoding='utf-8') as f:
                    yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
            elif format.lower() == 'json':
                with open(path, 'w', encoding='utf-8') as f:
                    json.dump(config_dict, f, indent=2, sort_keys=False)
            else:
                raise ValueError(f"Unsupported format: {format}")
        except Exception as e:
            raise Exception(f"Failed to save config to {file_path}: {e}")
    
    def validate(self) -> list:
        """
        Validate configuration values.
        
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Validate port
        if not (1 <= self.config.port <= 65535):
            errors.append(f"Invalid port: {self.config.port} (must be 1-65535)")
        
        # Validate max_workers
        if self.config.max_workers < 1:
            errors.append(f"Invalid max_workers: {self.config.max_workers} (must be >= 1)")
        
        # Validate log_level
        valid_log_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
        if self.config.log_level.upper() not in valid_log_levels:
            errors.append(f"Invalid log_level: {self.config.log_level}")
        
        # Validate upload_max_size
        if self.config.upload_max_size < 0:
            errors.append(f"Invalid upload_max_size: {self.config.upload_max_size}")
        
        # Validate rate limits
        if self.config.rate_limit_default < 1:
            errors.append("rate_limit_default must be >= 1")
        
        if self.config.rate_limit_window < 1:
            errors.append("rate_limit_window must be >= 1")
        
        return errors


# Global configuration manager instance
_config_manager: Optional[ConfigManager] = None


def get_config_manager(config_file: Optional[str] = None) -> ConfigManager:
    """
    Get global configuration manager instance.
    
    Args:
        config_file: Optional path to configuration file
    
    Returns:
        ConfigManager instance
    """
    global _config_manager
    if _config_manager is None:
        config_file = config_file or os.getenv('COMPASS_CONFIG_FILE')
        _config_manager = ConfigManager(config_file)
    return _config_manager


def get_service_config() -> Dict[str, Any]:
    """
    Get service configuration as dictionary (backward compatibility).
    
    Returns:
        Configuration dictionary
    """
    manager = get_config_manager()
    return manager.get_dict()

