"""
端到端测试配置
"""
import pytest
import sys
import time
from pathlib import Path

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

# 标记所有E2E测试
pytestmark = pytest.mark.e2e


@pytest.fixture(scope="session")
def flashdock_url():
    """FlashDock应用URL"""
    return "http://localhost:8501"


@pytest.fixture(scope="session")
def compass_url():
    """COMPASS服务URL"""
    return "http://localhost:8080"


@pytest.fixture(scope="session")
def registry_url():
    """服务注册中心URL"""
    return "http://localhost:8500"


@pytest.fixture
def wait_time():
    """等待时间（秒）"""
    return 2.0

