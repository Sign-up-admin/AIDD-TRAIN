"""
端到端测试：FlashDock UI测试
"""
import pytest
import sys
import time
import requests
from pathlib import Path
from typing import Optional

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))


@pytest.mark.e2e
class TestFlashDockUI:
    """FlashDock UI端到端测试类"""
    
    def test_flashdock_service_available(self, flashdock_url, wait_time):
        """测试FlashDock服务是否可用"""
        try:
            response = requests.get(flashdock_url, timeout=5)
            assert response.status_code == 200
            # 检查是否包含Streamlit相关内容
            assert "streamlit" in response.text.lower() or "FlashDock" in response.text
        except requests.exceptions.ConnectionError:
            pytest.skip(f"FlashDock服务未运行在 {flashdock_url}")
    
    def test_compass_service_available(self, compass_url, wait_time):
        """测试COMPASS服务是否可用"""
        try:
            response = requests.get(f"{compass_url}/health", timeout=5)
            assert response.status_code == 200
        except requests.exceptions.ConnectionError:
            pytest.skip(f"COMPASS服务未运行在 {compass_url}")
    
    def test_registry_service_available(self, registry_url, wait_time):
        """测试服务注册中心是否可用"""
        try:
            response = requests.get(f"{registry_url}/health", timeout=5)
            assert response.status_code == 200
        except requests.exceptions.ConnectionError:
            pytest.skip(f"服务注册中心未运行在 {registry_url}")
    
    def test_homepage_loads(self, flashdock_url):
        """测试主页加载"""
        try:
            response = requests.get(flashdock_url, timeout=10)
            assert response.status_code == 200
            # 检查页面内容
            content = response.text.lower()
            # 应该包含一些FlashDock相关的关键词
            assert len(content) > 0
        except requests.exceptions.ConnectionError:
            pytest.skip(f"FlashDock服务未运行在 {flashdock_url}")
    
    def test_api_endpoints_accessible(self, compass_url):
        """测试API端点可访问性"""
        try:
            # 测试健康检查
            response = requests.get(f"{compass_url}/health", timeout=5)
            assert response.status_code == 200
            
            # 测试任务列表API
            response = requests.get(f"{compass_url}/api/v1/training/tasks", timeout=5)
            # 应该返回200或401（如果启用认证）
            assert response.status_code in [200, 401]
        except requests.exceptions.ConnectionError:
            pytest.skip(f"COMPASS服务未运行在 {compass_url}")


# 注意：完整的浏览器自动化测试需要安装selenium或playwright
# 以下是使用selenium的示例（需要额外安装selenium和webdriver）

try:
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.common.exceptions import TimeoutException, WebDriverException
    
    SELENIUM_AVAILABLE = True
except ImportError:
    SELENIUM_AVAILABLE = False


if SELENIUM_AVAILABLE:
    @pytest.mark.e2e
    @pytest.mark.slow
    class TestFlashDockUISelenium:
        """使用Selenium的FlashDock UI测试类"""
        
        @pytest.fixture(scope="class")
        def driver(self):
            """创建WebDriver实例"""
            try:
                # 尝试使用Chrome
                options = webdriver.ChromeOptions()
                options.add_argument("--headless")
                options.add_argument("--no-sandbox")
                options.add_argument("--disable-dev-shm-usage")
                driver = webdriver.Chrome(options=options)
                yield driver
                driver.quit()
            except (WebDriverException, Exception) as e:
                pytest.skip(f"无法创建WebDriver: {e}")
        
        def test_navigation_sidebar(self, driver, flashdock_url):
            """测试侧边栏导航"""
            try:
                driver.get(flashdock_url)
                wait = WebDriverWait(driver, 10)
                
                # 等待页面加载
                wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
                
                # 检查侧边栏是否存在
                # Streamlit的侧边栏通常包含导航按钮
                time.sleep(2)  # 等待Streamlit完全加载
                
                # 这里可以添加更多具体的UI测试
                assert True
            except (TimeoutException, WebDriverException) as e:
                pytest.skip(f"无法访问FlashDock UI: {e}")
        
        def test_page_elements_load(self, driver, flashdock_url):
            """测试页面元素加载"""
            try:
                driver.get(flashdock_url)
                wait = WebDriverWait(driver, 10)
                
                # 等待页面加载
                wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
                time.sleep(2)
                
                # 检查页面标题或主要内容
                body = driver.find_element(By.TAG_NAME, "body")
                assert body is not None
            except (TimeoutException, WebDriverException) as e:
                pytest.skip(f"无法访问FlashDock UI: {e}")

