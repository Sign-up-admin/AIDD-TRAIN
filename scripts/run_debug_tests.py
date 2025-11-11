"""
Main script to run all debug tests and generate reports.
"""

import asyncio
import logging
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from scripts.chrome_debugger import ChromeDebugger, main as run_chrome_debugger
from scripts.vulnerability_scanner import VulnerabilityScanner
from scripts.report_generator import ReportGenerator

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def check_services_running():
    """Check if required services are running."""
    import requests
    from config.debug_config import BACKEND_URL, FRONTEND_URL, REGISTRY_URL

    services_ok = True

    # Check backend
    try:
        response = requests.get(f"{BACKEND_URL}/health", timeout=5)
        if response.status_code == 200:
            logger.info("Backend service is running")
        else:
            logger.warning(f"Backend service returned status {response.status_code}")
            services_ok = False
    except Exception as e:
        logger.error(f"Backend service is not running: {e}")
        services_ok = False

    # Check frontend
    try:
        response = requests.get(FRONTEND_URL, timeout=5)
        if response.status_code == 200:
            logger.info("Frontend service is running")
        else:
            logger.warning(f"Frontend service returned status {response.status_code}")
            services_ok = False
    except Exception as e:
        logger.error(f"Frontend service is not running: {e}")
        services_ok = False

    # Check registry (optional)
    try:
        response = requests.get(f"{REGISTRY_URL}/health", timeout=5)
        if response.status_code == 200:
            logger.info("Registry service is running")
        else:
            logger.warning(f"Registry service returned status {response.status_code}")
    except Exception as e:
        logger.warning(f"Registry service is not running: {e}")

    return services_ok


def main():
    """Main function to run all tests and generate reports."""
    logger.info("Starting debug tests and vulnerability scanning")

    # Check if services are running
    if not check_services_running():
        logger.error("Required services are not running. Please start them first.")
        logger.info("Start services with: start_all_services.bat")
        return 1

    try:
        # Run Chrome debugger
        logger.info("Running Chrome debugger tests...")
        asyncio.run(run_chrome_debugger())

        # Run vulnerability scanner
        logger.info("Running vulnerability scanner...")
        scanner = VulnerabilityScanner()
        scanner.scan_all()

        # Generate reports
        logger.info("Generating reports...")
        generator = ReportGenerator()
        generator.generate_all_reports()

        logger.info("All tests completed successfully!")
        logger.info("Check reports/ directory for detailed results")

        return 0

    except Exception as e:
        logger.error(f"Error running tests: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())





