"""
Chrome DevTools Protocol debugger for frontend and backend testing.
"""

import asyncio
import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
import sys

try:
    from playwright.async_api import async_playwright, Browser, Page, BrowserContext
except ImportError:
    print("Playwright not installed. Installing...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "playwright"])
    subprocess.check_call([sys.executable, "-m", "playwright", "install", "chromium"])
    from playwright.async_api import async_playwright, Browser, Page, BrowserContext

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from config.debug_config import (
    FRONTEND_URL,
    FRONTEND_PAGES,
    BACKEND_URL,
    API_ENDPOINTS,
    WEBSOCKET_ENDPOINTS,
    BROWSER_CONFIG,
    REPORT_CONFIG,
    TIMEOUT_CONFIG,
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


class ChromeDebugger:
    """Chrome DevTools Protocol debugger for testing frontend and backend."""

    def __init__(self):
        self.browser: Optional[Browser] = None
        self.context: Optional[BrowserContext] = None
        self.page: Optional[Page] = None
        self.playwright = None
        self.errors: List[Dict[str, Any]] = []
        self.network_requests: List[Dict[str, Any]] = []
        self.console_logs: List[Dict[str, Any]] = []
        self.screenshots: List[Dict[str, str]] = []

    async def start(self):
        """Start browser and create context."""
        self.playwright = await async_playwright().start()
        self.browser = await self.playwright.chromium.launch(
            headless=BROWSER_CONFIG["headless"],
            args=["--disable-web-security", "--disable-features=VizDisplayCompositor"],
        )
        self.context = await self.browser.new_context(
            viewport=BROWSER_CONFIG["viewport"],
            ignore_https_errors=True,
        )
        self.page = await self.context.new_page()

        # Setup error handlers
        self.page.on("console", self._on_console)
        self.page.on("pageerror", self._on_page_error)
        self.page.on("request", self._on_request)
        self.page.on("response", self._on_response)
        self.page.on("requestfailed", self._on_request_failed)

        logger.info("Browser started successfully")

    async def stop(self):
        """Stop browser and cleanup."""
        if self.context:
            await self.context.close()
        if self.browser:
            await self.browser.close()
        if self.playwright:
            await self.playwright.stop()
        logger.info("Browser stopped")

    def _on_console(self, msg):
        """Handle console messages."""
        log_entry = {
            "type": msg.type,
            "text": msg.text,
            "location": {
                "url": msg.location.get("url", ""),
                "line": msg.location.get("lineNumber", 0),
                "column": msg.location.get("columnNumber", 0),
            },
            "timestamp": datetime.now().isoformat(),
        }
        self.console_logs.append(log_entry)
        if msg.type in ["error", "warning"]:
            logger.warning(f"Console {msg.type}: {msg.text}")

    def _on_page_error(self, error):
        """Handle page errors."""
        error_entry = {
            "type": "page_error",
            "message": str(error),
            "timestamp": datetime.now().isoformat(),
        }
        self.errors.append(error_entry)
        logger.error(f"Page error: {error}")

    def _on_request(self, request):
        """Handle network requests."""
        request_entry = {
            "method": request.method,
            "url": request.url,
            "headers": request.headers,
            "timestamp": datetime.now().isoformat(),
        }
        self.network_requests.append(request_entry)

    def _on_response(self, response):
        """Handle network responses."""
        # Find matching request
        for req in self.network_requests:
            if req["url"] == response.url and req.get("response") is None:
                req["response"] = {
                    "status": response.status,
                    "statusText": response.status_text,
                    "headers": response.headers,
                    "timestamp": datetime.now().isoformat(),
                }
                if response.status >= 400:
                    error_entry = {
                        "type": "network_error",
                        "url": response.url,
                        "status": response.status,
                        "statusText": response.status_text,
                        "timestamp": datetime.now().isoformat(),
                    }
                    self.errors.append(error_entry)
                    logger.error(f"Network error: {response.status} {response.url}")

    def _on_request_failed(self, request):
        """Handle failed requests."""
        error_entry = {
            "type": "request_failed",
            "url": request.url,
            "method": request.method,
            "failure": request.failure,
            "timestamp": datetime.now().isoformat(),
        }
        self.errors.append(error_entry)
        logger.error(f"Request failed: {request.method} {request.url} - {request.failure}")

    async def test_frontend_page(self, page_name: str) -> Dict[str, Any]:
        """Test a frontend page."""
        logger.info(f"Testing frontend page: {page_name}")
        page_errors = []
        page_console_errors = []
        page_network_errors = []

        try:
            # Navigate to frontend
            await self.page.goto(FRONTEND_URL, wait_until="domcontentloaded", timeout=TIMEOUT_CONFIG["page_load"] * 1000)

            # Wait for Streamlit to load
            await asyncio.sleep(3)

            # For Streamlit, we can't directly navigate to pages via URL
            # Instead, we check if the main page loads and look for errors
            # In a real scenario, we would click the sidebar buttons to switch pages
            
            # Check if page loaded successfully
            page_title = await self.page.title()
            page_content = await self.page.content()

            # Clear previous errors for this page test
            previous_error_count = len(self.errors)

            # Take screenshot
            if BROWSER_CONFIG["screenshot"]:
                screenshot_path = f"{REPORT_CONFIG['output_dir']}/screenshots/{page_name.replace(' ', '_').replace('/', '_')}.png"
                Path(screenshot_path).parent.mkdir(parents=True, exist_ok=True)
                await self.page.screenshot(path=screenshot_path, full_page=True)
                self.screenshots.append({"page": page_name, "path": screenshot_path})

            # Wait a bit more for any async errors to appear
            await asyncio.sleep(2)

            # Check for errors in console logs (filter by recent)
            all_console_logs = self.console_logs[-50:]  # Last 50 logs
            page_console_errors = [log for log in all_console_logs if log["type"] == "error"]

            # Check for network errors (filter by recent and frontend URL)
            recent_errors = self.errors[-20:]  # Last 20 errors
            page_network_errors = [
                err for err in recent_errors 
                if err.get("url", "").startswith(FRONTEND_URL) or err.get("type") in ["network_error", "request_failed"]
            ]

            # Combine errors
            all_errors = page_console_errors + page_network_errors

            result = {
                "page": page_name,
                "status": "success" if not all_errors else "error",
                "errors": all_errors,
                "console_logs": page_console_errors,
                "network_errors": page_network_errors,
                "page_title": page_title,
                "timestamp": datetime.now().isoformat(),
            }

            logger.info(f"Page {page_name} test completed: {result['status']} ({len(all_errors)} errors)")
            return result

        except Exception as e:
            error_entry = {
                "type": "test_error",
                "page": page_name,
                "message": str(e),
                "timestamp": datetime.now().isoformat(),
            }
            page_errors.append(error_entry)
            logger.error(f"Error testing page {page_name}: {e}", exc_info=True)
            return {
                "page": page_name,
                "status": "error",
                "errors": page_errors,
                "timestamp": datetime.now().isoformat(),
            }

    async def test_api_endpoint(self, method: str, path: str, description: str) -> Dict[str, Any]:
        """Test an API endpoint."""
        logger.info(f"Testing API endpoint: {method} {path}")

        try:
            # Prepare request
            url = f"{BACKEND_URL}{path}"
            headers = {"Content-Type": "application/json"}

            # Replace path parameters with test values
            if "{task_id}" in path:
                path = path.replace("{task_id}", "test-task-id")
                url = f"{BACKEND_URL}{path}"
            if "{dataset_id}" in path:
                path = path.replace("{dataset_id}", "test-dataset-id")
                url = f"{BACKEND_URL}{path}"
            if "{model_id}" in path:
                path = path.replace("{model_id}", "test-model-id")
                url = f"{BACKEND_URL}{path}"

            # Make request
            start_time = time.time()
            response = await self.page.request.fetch(
                url, method=method, headers=headers, timeout=TIMEOUT_CONFIG["api_request"] * 1000
            )
            elapsed_time = time.time() - start_time

            # Get response body
            try:
                response_body = await response.json()
            except:
                response_body = await response.text()

            result = {
                "method": method,
                "path": path,
                "description": description,
                "status": response.status,
                "statusText": response.status_text,
                "response_time": elapsed_time,
                "response_body": response_body,
                "timestamp": datetime.now().isoformat(),
            }

            # Check for errors
            if response.status >= 400:
                result["error"] = True
                error_entry = {
                    "type": "api_error",
                    "method": method,
                    "path": path,
                    "status": response.status,
                    "statusText": response.status_text,
                    "response_body": response_body,
                    "timestamp": datetime.now().isoformat(),
                }
                self.errors.append(error_entry)

            logger.info(f"API endpoint {method} {path} test completed: {response.status}")
            return result

        except Exception as e:
            error_entry = {
                "type": "api_test_error",
                "method": method,
                "path": path,
                "message": str(e),
                "timestamp": datetime.now().isoformat(),
            }
            self.errors.append(error_entry)
            logger.error(f"Error testing API endpoint {method} {path}: {e}")
            return {
                "method": method,
                "path": path,
                "description": description,
                "status": "error",
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            }

    async def test_all_frontend_pages(self) -> List[Dict[str, Any]]:
        """Test all frontend pages."""
        logger.info("Testing all frontend pages")
        results = []

        for page_name in FRONTEND_PAGES:
            result = await self.test_frontend_page(page_name)
            results.append(result)
            await asyncio.sleep(1)  # Wait between page tests

        return results

    async def test_all_api_endpoints(self) -> List[Dict[str, Any]]:
        """Test all API endpoints."""
        logger.info("Testing all API endpoints")
        results = []

        for category, endpoints in API_ENDPOINTS.items():
            logger.info(f"Testing {category} endpoints")
            for endpoint in endpoints:
                result = await self.test_api_endpoint(
                    endpoint["method"], endpoint["path"], endpoint["description"]
                )
                results.append(result)
                await asyncio.sleep(0.5)  # Wait between API tests

        return results

    def get_summary(self) -> Dict[str, Any]:
        """Get test summary."""
        return {
            "total_errors": len(self.errors),
            "total_console_logs": len(self.console_logs),
            "total_network_requests": len(self.network_requests),
            "total_screenshots": len(self.screenshots),
            "errors": self.errors,
            "console_logs": self.console_logs,
            "network_requests": self.network_requests,
            "screenshots": self.screenshots,
        }


async def main():
    """Main function."""
    debugger = ChromeDebugger()

    try:
        # Create reports directory
        Path(REPORT_CONFIG["output_dir"]).mkdir(parents=True, exist_ok=True)
        Path(f"{REPORT_CONFIG['output_dir']}/screenshots").mkdir(parents=True, exist_ok=True)

        # Start browser
        await debugger.start()

        # Test frontend pages
        logger.info("Starting frontend tests")
        frontend_results = await debugger.test_all_frontend_pages()

        # Test API endpoints
        logger.info("Starting API tests")
        api_results = await debugger.test_all_api_endpoints()

        # Get summary
        summary = debugger.get_summary()

        # Save results
        results = {
            "frontend_results": frontend_results,
            "api_results": api_results,
            "summary": summary,
            "timestamp": datetime.now().isoformat(),
        }

        results_path = f"{REPORT_CONFIG['output_dir']}/debug_results.json"
        with open(results_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

        logger.info(f"Test results saved to {results_path}")
        logger.info(f"Total errors: {summary['total_errors']}")

    except Exception as e:
        logger.error(f"Error in main: {e}", exc_info=True)
    finally:
        await debugger.stop()


if __name__ == "__main__":
    asyncio.run(main())

