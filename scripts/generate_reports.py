#!/usr/bin/env python3
"""
Generate unified reports from lint check results.
Combines individual tool reports into a comprehensive HTML and Markdown report.
"""

import os
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List

PROJECT_ROOT = Path(__file__).parent.parent
LINT_REPORTS_DIR = PROJECT_ROOT / "lint_reports"


def parse_report_file(report_file: Path) -> Dict:
    """Parse a report file and extract key information."""
    if not report_file.exists():
        return None
    
    content = report_file.read_text(encoding="utf-8", errors="replace")
    
    # Extract tool name from filename
    tool_name = report_file.stem.split("_")[-1].replace("_report", "")
    
    # Try to extract exit code
    exit_code = 0
    if "Exit Code:" in content:
        try:
            exit_code_line = [line for line in content.split("\n") if "Exit Code:" in line][0]
            exit_code = int(exit_code_line.split(":")[-1].strip())
        except (ValueError, IndexError):
            pass
    
    # Count issues
    issue_count = 0
    if "Summary of Issues:" in content:
        issue_section = content.split("Summary of Issues:")[1]
        issue_count = len([line for line in issue_section.split("\n") if line.strip().startswith("-")])
    
    return {
        "tool": tool_name,
        "file": report_file.name,
        "exit_code": exit_code,
        "issue_count": issue_count,
        "content": content,
        "timestamp": report_file.stat().st_mtime,
    }


def generate_html_report(reports: List[Dict], output_file: Path):
    """Generate an HTML report from all reports."""
    html = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Code Quality Report</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            color: #333;
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 10px;
        }
        h2 {
            color: #555;
            margin-top: 30px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #4CAF50;
            color: white;
        }
        tr:hover {
            background-color: #f5f5f5;
        }
        .status-pass {
            color: #4CAF50;
            font-weight: bold;
        }
        .status-fail {
            color: #f44336;
            font-weight: bold;
        }
        .pass-icon::before {
            content: "✓ ";
        }
        .fail-icon::before {
            content: "✗ ";
        }
        .report-section {
            margin: 30px 0;
            padding: 20px;
            background-color: #f9f9f9;
            border-radius: 4px;
            border-left: 4px solid #2196F3;
        }
        .report-content {
            background: white;
            padding: 15px;
            border-radius: 4px;
            font-family: 'Courier New', monospace;
            font-size: 12px;
            overflow-x: auto;
            white-space: pre-wrap;
            max-height: 500px;
            overflow-y: auto;
        }
        .summary {
            display: flex;
            gap: 20px;
            margin: 20px 0;
        }
        .summary-card {
            flex: 1;
            padding: 20px;
            border-radius: 4px;
            text-align: center;
        }
        .summary-card.total { background-color: #e3f2fd; }
        .summary-card.passed { background-color: #e8f5e9; }
        .summary-card.failed { background-color: #ffebee; }
        .summary-card h3 {
            margin: 0;
            font-size: 36px;
            color: #333;
        }
        .summary-card p {
            margin: 5px 0 0 0;
            color: #666;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Code Quality Report</h1>
        <p><strong>Generated:</strong> {timestamp}</p>
        
        <div class="summary">
            <div class="summary-card total">
                <h3>{total_tools}</h3>
                <p>Total Tools</p>
            </div>
            <div class="summary-card passed">
                <h3>{passed_tools}</h3>
                <p>Passed</p>
            </div>
            <div class="summary-card failed">
                <h3>{failed_tools}</h3>
                <p>Failed</p>
            </div>
        </div>
        
        <h2>Tool Status</h2>
        <table>
            <thead>
                <tr>
                    <th>Tool</th>
                    <th>Status</th>
                    <th>Issues</th>
                    <th>Report File</th>
                </tr>
            </thead>
            <tbody>
                {tool_rows}
            </tbody>
        </table>
        
        <h2>Detailed Reports</h2>
        {report_sections}
    </div>
</body>
</html>"""
    
    total_tools = len(reports)
    passed_tools = sum(1 for r in reports if r['exit_code'] == 0)
    failed_tools = total_tools - passed_tools
    
    tool_rows = ""
    for report in reports:
        status_class = "status-pass" if report['exit_code'] == 0 else "status-fail"
        status_text = "PASS" if report['exit_code'] == 0 else "FAIL"
        tool_rows += f"""
                <tr>
                    <td><strong>{report['tool'].upper()}</strong></td>
                    <td class="{status_class}"><span class="{'pass-icon' if report['exit_code'] == 0 else 'fail-icon'}"></span>{status_text}</td>
                    <td>{report['issue_count']}</td>
                    <td>{report['file']}</td>
                </tr>"""
    
    report_sections = ""
    for report in reports:
        # Truncate content for display
        content = report['content'][:5000] + "..." if len(report['content']) > 5000 else report['content']
        report_sections += f"""
        <div class="report-section">
            <h3>{report['tool'].upper()}</h3>
            <div class="report-content">{content}</div>
        </div>"""
    
    html_content = html.format(
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        total_tools=total_tools,
        passed_tools=passed_tools,
        failed_tools=failed_tools,
        tool_rows=tool_rows,
        report_sections=report_sections,
    )
    
    output_file.write_text(html_content, encoding="utf-8")


def main():
    """Generate unified reports."""
    # Find all report files
    report_files = list(LINT_REPORTS_DIR.glob("*_report.txt"))
    
    if not report_files:
        print("No report files found in lint_reports/")
        return
    
    # Parse all reports
    reports = []
    for report_file in sorted(report_files, key=lambda x: x.stat().st_mtime, reverse=True):
        # Get the most recent set of reports (same timestamp)
        report = parse_report_file(report_file)
        if report:
            reports.append(report)
    
    if not reports:
        print("No valid reports found")
        return
    
    # Generate HTML report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    html_file = LINT_REPORTS_DIR / f"{timestamp}_unified_report.html"
    generate_html_report(reports, html_file)
    
    print(f"[OK] Generated unified HTML report: {html_file}")
    print(f"   Open in browser to view: file:///{html_file.absolute()}")


if __name__ == "__main__":
    main()

