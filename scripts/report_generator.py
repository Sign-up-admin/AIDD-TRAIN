"""
Report generator for debugging and vulnerability scan results.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any
import sys

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from config.debug_config import REPORT_CONFIG

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


class ReportGenerator:
    """Generator for debug and vulnerability reports."""

    def __init__(self):
        self.output_dir = Path(REPORT_CONFIG["output_dir"])
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_debug_results(self) -> Dict[str, Any]:
        """Load debug results from JSON file."""
        results_path = self.output_dir / "debug_results.json"
        if results_path.exists():
            with open(results_path, "r", encoding="utf-8") as f:
                return json.load(f)
        return {}

    def load_vulnerability_results(self) -> Dict[str, Any]:
        """Load vulnerability scan results from JSON file."""
        results_path = self.output_dir / "vulnerability_scan_results.json"
        if results_path.exists():
            with open(results_path, "r", encoding="utf-8") as f:
                return json.load(f)
        return {}

    def generate_frontend_error_report(self, debug_results: Dict[str, Any]) -> str:
        """Generate frontend error report."""
        report = ["# Frontend Error Report\n"]
        report.append(f"**Generated:** {datetime.now().isoformat()}\n\n")

        frontend_results = debug_results.get("frontend_results", [])
        if not frontend_results:
            report.append("No frontend test results found.\n")
            return "".join(report)

        # Count errors by page
        error_count = 0
        for result in frontend_results:
            if result.get("status") == "error" or result.get("errors"):
                error_count += 1

        report.append(f"## Summary\n\n")
        report.append(f"- **Total pages tested:** {len(frontend_results)}\n")
        report.append(f"- **Pages with errors:** {error_count}\n")
        report.append(f"- **Total errors:** {debug_results.get('summary', {}).get('total_errors', 0)}\n\n")

        # Detailed errors by page
        report.append("## Errors by Page\n\n")
        for result in frontend_results:
            if result.get("status") == "error" or result.get("errors"):
                report.append(f"### {result['page']}\n\n")
                report.append(f"**Status:** {result.get('status', 'unknown')}\n\n")

                errors = result.get("errors", [])
                if errors:
                    report.append("**Errors:**\n\n")
                    for error in errors:
                        report.append(f"- **Type:** {error.get('type', 'unknown')}\n")
                        report.append(f"  - **Message:** {error.get('message', 'N/A')}\n")
                        report.append(f"  - **Timestamp:** {error.get('timestamp', 'N/A')}\n\n")

                console_logs = result.get("console_logs", [])
                if console_logs:
                    report.append("**Console Logs:**\n\n")
                    for log in console_logs:
                        report.append(f"- **{log.get('type', 'unknown')}:** {log.get('text', 'N/A')}\n")
                        report.append(f"  - **Location:** {log.get('location', {}).get('url', 'N/A')}\n\n")

                network_errors = result.get("network_errors", [])
                if network_errors:
                    report.append("**Network Errors:**\n\n")
                    for error in network_errors:
                        report.append(f"- **URL:** {error.get('url', 'N/A')}\n")
                        report.append(f"  - **Status:** {error.get('status', 'N/A')}\n")
                        report.append(f"  - **Status Text:** {error.get('statusText', 'N/A')}\n\n")

        return "".join(report)

    def generate_backend_error_report(self, debug_results: Dict[str, Any]) -> str:
        """Generate backend API error report."""
        report = ["# Backend API Error Report\n"]
        report.append(f"**Generated:** {datetime.now().isoformat()}\n\n")

        api_results = debug_results.get("api_results", [])
        if not api_results:
            report.append("No API test results found.\n")
            return "".join(report)

        # Count errors
        error_count = sum(1 for r in api_results if r.get("status") >= 400 or r.get("error"))
        slow_requests = [r for r in api_results if r.get("response_time", 0) > 5.0]

        report.append(f"## Summary\n\n")
        report.append(f"- **Total endpoints tested:** {len(api_results)}\n")
        report.append(f"- **Endpoints with errors:** {error_count}\n")
        report.append(f"- **Slow requests (>5s):** {len(slow_requests)}\n\n")

        # Errors by endpoint
        report.append("## Errors by Endpoint\n\n")
        for result in api_results:
            if result.get("status") >= 400 or result.get("error"):
                report.append(f"### {result.get('method', 'UNKNOWN')} {result.get('path', 'N/A')}\n\n")
                report.append(f"**Description:** {result.get('description', 'N/A')}\n\n")
                report.append(f"**Status:** {result.get('status', 'N/A')}\n")
                report.append(f"**Status Text:** {result.get('statusText', 'N/A')}\n")
                if result.get("error"):
                    report.append(f"**Error:** {result.get('error', 'N/A')}\n")
                if result.get("response_body"):
                    report.append(f"**Response Body:**\n```json\n{json.dumps(result.get('response_body'), indent=2)}\n```\n\n")

        # Slow requests
        if slow_requests:
            report.append("## Slow Requests\n\n")
            for result in slow_requests:
                report.append(f"- **{result.get('method', 'UNKNOWN')} {result.get('path', 'N/A')}:** {result.get('response_time', 0):.2f}s\n")

        return "".join(report)

    def generate_security_vulnerability_report(self, vuln_results: Dict[str, Any]) -> str:
        """Generate security vulnerability report."""
        report = ["# Security Vulnerability Report\n"]
        report.append(f"**Generated:** {datetime.now().isoformat()}\n\n")

        vulnerabilities = vuln_results.get("vulnerabilities", [])
        if not vulnerabilities:
            report.append("No security vulnerabilities found.\n")
            return "".join(report)

        # Group by severity
        high_severity = [v for v in vulnerabilities if v.get("severity") == "high"]
        medium_severity = [v for v in vulnerabilities if v.get("severity") == "medium"]
        low_severity = [v for v in vulnerabilities if v.get("severity") == "low"]

        report.append(f"## Summary\n\n")
        report.append(f"- **Total vulnerabilities:** {len(vulnerabilities)}\n")
        report.append(f"- **High severity:** {len(high_severity)}\n")
        report.append(f"- **Medium severity:** {len(medium_severity)}\n")
        report.append(f"- **Low severity:** {len(low_severity)}\n\n")

        # Vulnerabilities by category
        report.append("## Vulnerabilities by Category\n\n")

        categories = {}
        for vuln in vulnerabilities:
            category = vuln.get("category", "unknown")
            if category not in categories:
                categories[category] = []
            categories[category].append(vuln)

        for category, vulns in categories.items():
            report.append(f"### {category.upper()}\n\n")
            for vuln in vulns:
                report.append(f"#### {vuln.get('title', 'Unknown')}\n\n")
                report.append(f"**Severity:** {vuln.get('severity', 'unknown')}\n\n")
                report.append(f"**Description:** {vuln.get('description', 'N/A')}\n\n")
                report.append(f"**Location:** {vuln.get('location', 'N/A')}\n\n")
                report.append(f"**Recommendation:** {vuln.get('recommendation', 'N/A')}\n\n")
                report.append("---\n\n")

        return "".join(report)

    def generate_performance_issues_report(self, vuln_results: Dict[str, Any]) -> str:
        """Generate performance issues report."""
        report = ["# Performance Issues Report\n"]
        report.append(f"**Generated:** {datetime.now().isoformat()}\n\n")

        performance_issues = vuln_results.get("performance_issues", [])
        if not performance_issues:
            report.append("No performance issues found.\n")
            return "".join(report)

        report.append(f"## Summary\n\n")
        report.append(f"- **Total performance issues:** {len(performance_issues)}\n\n")

        # Issues by category
        report.append("## Performance Issues\n\n")
        for issue in performance_issues:
            report.append(f"### {issue.get('title', 'Unknown')}\n\n")
            report.append(f"**Severity:** {issue.get('severity', 'unknown')}\n\n")
            report.append(f"**Description:** {issue.get('description', 'N/A')}\n\n")
            report.append(f"**Location:** {issue.get('location', 'N/A')}\n\n")
            if issue.get("response_time"):
                report.append(f"**Response Time:** {issue.get('response_time', 0):.2f}s\n\n")
            report.append(f"**Recommendation:** {issue.get('recommendation', 'N/A')}\n\n")
            report.append("---\n\n")

        return "".join(report)

    def generate_comprehensive_report(
        self, debug_results: Dict[str, Any], vuln_results: Dict[str, Any]
    ) -> str:
        """Generate comprehensive report."""
        report = ["# Comprehensive Debug and Vulnerability Report\n"]
        report.append(f"**Generated:** {datetime.now().isoformat()}\n\n")

        report.append("## Executive Summary\n\n")
        report.append("This report contains the results of automated debugging and vulnerability scanning.\n\n")

        # Summary statistics
        debug_summary = debug_results.get("summary", {})
        report.append("### Debug Results\n\n")
        report.append(f"- **Total errors:** {debug_summary.get('total_errors', 0)}\n")
        report.append(f"- **Total console logs:** {debug_summary.get('total_console_logs', 0)}\n")
        report.append(f"- **Total network requests:** {debug_summary.get('total_network_requests', 0)}\n\n")

        vulnerabilities = vuln_results.get("vulnerabilities", [])
        performance_issues = vuln_results.get("performance_issues", [])
        report.append("### Vulnerability Scan Results\n\n")
        report.append(f"- **Total vulnerabilities:** {len(vulnerabilities)}\n")
        report.append(f"- **Total performance issues:** {len(performance_issues)}\n\n")

        # Priority issues
        high_severity_vulns = [v for v in vulnerabilities if v.get("severity") == "high"]
        if high_severity_vulns:
            report.append("## Priority Issues (P0)\n\n")
            report.append("The following high-severity issues require immediate attention:\n\n")
            for vuln in high_severity_vulns:
                report.append(f"1. **{vuln.get('title', 'Unknown')}** - {vuln.get('category', 'unknown')}\n")
                report.append(f"   - Location: {vuln.get('location', 'N/A')}\n")
                report.append(f"   - Recommendation: {vuln.get('recommendation', 'N/A')}\n\n")

        # Detailed sections
        report.append("## Detailed Reports\n\n")
        report.append("For detailed information, see the following reports:\n\n")
        report.append("- [Frontend Error Report](frontend_errors.md)\n")
        report.append("- [Backend Error Report](backend_errors.md)\n")
        report.append("- [Security Vulnerability Report](security_vulnerabilities.md)\n")
        report.append("- [Performance Issues Report](performance_issues.md)\n\n")

        return "".join(report)

    def generate_all_reports(self):
        """Generate all reports."""
        logger.info("Generating reports")

        # Load results
        debug_results = self.load_debug_results()
        vuln_results = self.load_vulnerability_results()

        # Generate reports
        frontend_report = self.generate_frontend_error_report(debug_results)
        backend_report = self.generate_backend_error_report(debug_results)
        security_report = self.generate_security_vulnerability_report(vuln_results)
        performance_report = self.generate_performance_issues_report(vuln_results)
        comprehensive_report = self.generate_comprehensive_report(debug_results, vuln_results)

        # Save reports
        reports = {
            "frontend_errors.md": frontend_report,
            "backend_errors.md": backend_report,
            "security_vulnerabilities.md": security_report,
            "performance_issues.md": performance_report,
            "comprehensive_report.md": comprehensive_report,
        }

        for filename, content in reports.items():
            report_path = self.output_dir / filename
            with open(report_path, "w", encoding="utf-8") as f:
                f.write(content)
            logger.info(f"Report saved: {report_path}")

        logger.info("All reports generated successfully")


def main():
    """Main function."""
    generator = ReportGenerator()
    generator.generate_all_reports()


if __name__ == "__main__":
    main()





