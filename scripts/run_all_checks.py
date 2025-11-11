#!/usr/bin/env python3
"""
Code Quality Check Script
Runs all code quality checks (Black, Flake8, Pylint, MyPy, Bandit, Pytest)
and generates unified reports.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional

# Fix Windows console encoding
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent
LINT_REPORTS_DIR = PROJECT_ROOT / "lint_reports"
LINT_REPORTS_DIR.mkdir(exist_ok=True)

# Directories and files to check
CHECK_DIRS = {
    "compass": ["compass"],
    "services": ["services"],
    "flashdock_services": ["FLASH_DOCK-main/services"],
    "flashdock_pages": ["FLASH_DOCK-main/pages"],
    "flashdock_main": ["FLASH_DOCK-main/FlashDock.py"],
    "scripts": ["scripts"],
    "root_scripts": [
        "check_and_start_services.py",
        "check_ports.py",
        "check_service_status.py",
    ],
}

# Tool configurations
TOOLS = {
    "black": {
        "cmd": ["black", "--check"],
        "report_file": "black_report.txt",
        "fix_cmd": ["black"],
        "description": "Code formatting check",
    },
    "flake8": {
        "cmd": ["flake8"],
        "report_file": "flake8_report.txt",
        "description": "Code style and complexity check",
    },
    "pylint": {
        "cmd": ["pylint"],
        "report_file": "pylint_report.txt",
        "description": "Code quality analysis",
    },
    "mypy": {
        "cmd": ["mypy"],
        "report_file": "mypy_report.txt",
        "description": "Type checking",
    },
    "bandit": {
        "cmd": ["bandit", "-r", "-f", "txt"],
        "report_file": "bandit_report.txt",
        "description": "Security vulnerability scan",
    },
    "pytest": {
        "cmd": ["pytest", "--cov=compass", "--cov-report=html", "--cov-report=xml", "--cov-report=term"],
        "report_file": "pytest_report.txt",
        "description": "Unit tests and coverage",
    },
}


def run_command(cmd: List[str], cwd: Path = None, capture_output: bool = True) -> Tuple[int, str, str]:
    """Run a shell command and return exit code, stdout, and stderr."""
    try:
        # For Python-based tools, always use sys.executable as module (most reliable in conda environments)
        if cmd[0] in ["black", "flake8", "pylint", "mypy", "bandit", "pytest"]:
            module_name = cmd[0]
            # Use python -m module_name instead of direct command
            module_cmd = [sys.executable, "-m", module_name] + cmd[1:]
            result = subprocess.run(
                module_cmd,
                cwd=cwd or PROJECT_ROOT,
                capture_output=capture_output,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=300,  # 5 minute timeout
            )
            return result.returncode, result.stdout, result.stderr
        
        # For other commands, run directly
        result = subprocess.run(
            cmd,
            cwd=cwd or PROJECT_ROOT,
            capture_output=capture_output,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=300,  # 5 minute timeout
        )
        return result.returncode, result.stdout, result.stderr
    except FileNotFoundError:
        return 1, "", f"Command not found: {cmd[0]}"
    except subprocess.TimeoutExpired:
        return 1, "", f"Command timed out: {cmd[0]}"
    except Exception as e:
        return 1, "", str(e)


def detect_conda_env() -> Optional[str]:
    """
    Detect if we're in a conda environment and return the environment name.
    Returns None if not in conda environment.
    """
    # Check if CONDA_DEFAULT_ENV is set
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    if conda_env:
        return conda_env
    
    # Check if Python is in a conda environment path
    python_path = Path(sys.executable)
    python_path_str = str(python_path).lower()
    if "conda" in python_path_str or "anaconda" in python_path_str:
        # Try to extract environment name from path
        for part in python_path.parts:
            if part in ["envs", "env"]:
                idx = python_path.parts.index(part)
                if idx + 1 < len(python_path.parts):
                    return python_path.parts[idx + 1]
    
    return None


def find_conda_python(env_name: str) -> Optional[str]:
    """
    Find Python executable in a conda environment.
    Returns the path to Python executable or None if not found.
    """
    # Common conda installation paths on Windows
    possible_conda_paths = [
        Path("C:/ProgramData/Anaconda3"),
        Path("C:/ProgramData/Miniconda3"),
        Path(os.path.expanduser("~/Anaconda3")),
        Path(os.path.expanduser("~/Miniconda3")),
    ]
    
    # Also check if conda is in PATH
    try:
        result = subprocess.run(
            ["conda", "info", "--base"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            conda_base = Path(result.stdout.strip())
            possible_conda_paths.insert(0, conda_base)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    # Check each possible path
    for conda_base in possible_conda_paths:
        env_path = conda_base / "envs" / env_name
        python_exe = env_path / "python.exe"
        if python_exe.exists():
            return str(python_exe)
    
    return None


def ensure_conda_environment() -> Tuple[bool, Optional[str]]:
    """
    Ensure we're running in a conda environment.
    Returns (True, None) if in conda environment, 
    (False, python_path) if found conda environment but not activated,
    (False, None) if no conda environment found.
    """
    conda_env = detect_conda_env()
    
    if conda_env:
        print(f"[INFO] Running in conda environment: {conda_env}")
        return True, None
    
    # Not in conda environment, try to find and use one
    print("[WARNING] Not running in a conda environment.")
    print("[INFO] Attempting to find conda environment...")
    
    # Try common environment names
    env_names = ["AIDDTRAIN", "flash_dock_new", "flash_dock", "AIDD"]
    
    for env_name in env_names:
        python_path = find_conda_python(env_name)
        if python_path:
            print(f"[INFO] Found conda environment: {env_name}")
            print(f"[INFO] Python path: {python_path}")
            print(f"[INFO] Will use conda environment Python for running checks.")
            # Update sys.executable to use conda Python
            return False, python_path
    
    print("[ERROR] Could not find any conda environment.")
    print("[INFO] Please activate a conda environment before running this script.")
    return False, None


def check_tool_available(tool_name: str) -> bool:
    """Check if a tool is available in PATH."""
    try:
        # Always use python -m module_name for Python tools (works best in conda environments)
        if tool_name in ["black", "flake8", "pylint", "mypy", "bandit", "pytest"]:
            subprocess.run(
                [sys.executable, "-m", tool_name, "--version"],
                capture_output=True,
                check=True,
                timeout=10
            )
        else:
            subprocess.run([tool_name, "--version"], capture_output=True, check=True, timeout=10)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False


def run_black_check(directories: Dict[str, List[str]], fix: bool = False) -> Dict:
    """Run Black formatting check."""
    tool_config = TOOLS["black"]
    cmd = tool_config["fix_cmd"] if fix else tool_config["cmd"]
    
    results = {}
    all_output = []
    all_errors = []
    total_exit_code = 0
    
    for dir_name, dirs in directories.items():
        for directory in dirs:
            dir_path = PROJECT_ROOT / directory
            if not dir_path.exists():
                continue
                
            full_cmd = cmd + [str(dir_path)]
            exit_code, stdout, stderr = run_command(full_cmd)
            
            output = f"\n{'='*80}\n"
            output += f"Checking: {directory}\n"
            output += f"{'='*80}\n"
            output += stdout + stderr
            
            all_output.append(output)
            if exit_code != 0:
                all_errors.append(f"{directory}: Formatting issues found")
                total_exit_code = max(total_exit_code, exit_code)
            
            results[directory] = {
                "exit_code": exit_code,
                "output": stdout + stderr,
            }
    
    return {
        "tool": "black",
        "exit_code": total_exit_code,
        "output": "\n".join(all_output),
        "errors": all_errors,
        "results": results,
    }


def run_flake8_check(directories: Dict[str, List[str]]) -> Dict:
    """Run Flake8 style check."""
    tool_config = TOOLS["flake8"]
    all_output = []
    all_errors = []
    total_exit_code = 0
    
    for dir_name, dirs in directories.items():
        for directory in dirs:
            dir_path = PROJECT_ROOT / directory
            if not dir_path.exists():
                continue
                
            cmd = tool_config["cmd"] + [str(dir_path)]
            exit_code, stdout, stderr = run_command(cmd)
            
            output = f"\n{'='*80}\n"
            output += f"Checking: {directory}\n"
            output += f"{'='*80}\n"
            output += stdout + stderr
            
            all_output.append(output)
            if exit_code != 0:
                error_count = len([line for line in stdout.split('\n') if line.strip()])
                all_errors.append(f"{directory}: {error_count} style issues")
                total_exit_code = max(total_exit_code, exit_code)
    
    return {
        "tool": "flake8",
        "exit_code": total_exit_code,
        "output": "\n".join(all_output),
        "errors": all_errors,
    }


def run_pylint_check(directories: Dict[str, List[str]]) -> Dict:
    """Run Pylint quality check."""
    tool_config = TOOLS["pylint"]
    all_output = []
    all_errors = []
    total_exit_code = 0
    
    for dir_name, dirs in directories.items():
        for directory in dirs:
            dir_path = PROJECT_ROOT / directory
            if not dir_path.exists():
                continue
                
            cmd = tool_config["cmd"] + [str(dir_path)]
            exit_code, stdout, stderr = run_command(cmd)
            
            output = f"\n{'='*80}\n"
            output += f"Checking: {directory}\n"
            output += f"{'='*80}\n"
            output += stdout + stderr
            
            all_output.append(output)
            if exit_code != 0:
                all_errors.append(f"{directory}: Code quality issues found")
                total_exit_code = max(total_exit_code, exit_code)
    
    return {
        "tool": "pylint",
        "exit_code": total_exit_code,
        "output": "\n".join(all_output),
        "errors": all_errors,
    }


def run_mypy_check(directories: Dict[str, List[str]]) -> Dict:
    """Run MyPy type check."""
    tool_config = TOOLS["mypy"]
    all_output = []
    all_errors = []
    total_exit_code = 0
    
    # MyPy should check Python packages and key scripts
    check_dirs = {
        "compass": ["compass"],
        "services": ["services"],
        "flashdock_services": ["FLASH_DOCK-main/services"],
        "scripts": ["scripts"],
        "root_scripts": [
            "check_and_start_services.py",
            "check_ports.py",
            "check_service_status.py",
        ],
    }
    
    for dir_name, dirs in check_dirs.items():
        for directory in dirs:
            dir_path = PROJECT_ROOT / directory
            if not dir_path.exists():
                continue
                
            cmd = tool_config["cmd"] + [str(dir_path)]
            exit_code, stdout, stderr = run_command(cmd)
            
            output = f"\n{'='*80}\n"
            output += f"Checking: {directory}\n"
            output += f"{'='*80}\n"
            output += stdout + stderr
            
            all_output.append(output)
            if exit_code != 0:
                error_count = len([line for line in stdout.split('\n') if 'error:' in line])
                all_errors.append(f"{directory}: {error_count} type errors")
                total_exit_code = max(total_exit_code, exit_code)
    
    return {
        "tool": "mypy",
        "exit_code": total_exit_code,
        "output": "\n".join(all_output),
        "errors": all_errors,
    }


def run_bandit_check(directories: Dict[str, List[str]]) -> Dict:
    """Run Bandit security check."""
    tool_config = TOOLS["bandit"]
    all_output = []
    all_errors = []
    total_exit_code = 0
    
    check_dirs = {
        "compass": ["compass"],
        "services": ["services"],
        "flashdock_services": ["FLASH_DOCK-main/services"],
        "flashdock_main": ["FLASH_DOCK-main/FlashDock.py"],
        "scripts": ["scripts"],
        "root_scripts": [
            "check_and_start_services.py",
            "check_ports.py",
            "check_service_status.py",
        ],
    }
    
    for dir_name, dirs in check_dirs.items():
        for directory in dirs:
            dir_path = PROJECT_ROOT / directory
            if not dir_path.exists():
                continue
                
            # Exclude others directory for FlashDock if checking directory
            if directory == "FLASH_DOCK-main":
                cmd = tool_config["cmd"] + ["-x", "others", str(dir_path)]
            else:
                cmd = tool_config["cmd"] + [str(dir_path)]
            
            exit_code, stdout, stderr = run_command(cmd)
            
            output = f"\n{'='*80}\n"
            output += f"Checking: {directory}\n"
            output += f"{'='*80}\n"
            output += stdout + stderr
            
            all_output.append(output)
            if exit_code != 0:
                # Check for high severity issues
                high_severity = stdout.count("Severity: High")
                medium_severity = stdout.count("Severity: Medium")
                if high_severity > 0:
                    all_errors.append(f"{directory}: {high_severity} high severity issues")
                    total_exit_code = max(total_exit_code, 2)
                elif medium_severity > 0:
                    total_exit_code = max(total_exit_code, 1)
    
    return {
        "tool": "bandit",
        "exit_code": total_exit_code,
        "output": "\n".join(all_output),
        "errors": all_errors,
    }


def run_pytest() -> Dict:
    """Run Pytest tests."""
    tool_config = TOOLS["pytest"]
    cmd = tool_config["cmd"]
    
    exit_code, stdout, stderr = run_command(cmd, capture_output=False)
    
    # Capture output for report
    exit_code_capture, stdout_capture, stderr_capture = run_command(cmd)
    
    output = f"\n{'='*80}\n"
    output += "Running Pytest Tests\n"
    output += f"{'='*80}\n"
    output += stdout_capture + stderr_capture
    
    errors = []
    if exit_code != 0:
        errors.append("Some tests failed")
    
    return {
        "tool": "pytest",
        "exit_code": exit_code,
        "output": output,
        "errors": errors,
    }


def save_report(tool_name: str, report_data: Dict, timestamp: str):
    """Save individual tool report to file."""
    tool_config = TOOLS[tool_name]
    report_file = LINT_REPORTS_DIR / f"{timestamp}_{tool_config['report_file']}"
    
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(f"Report generated at: {datetime.now().isoformat()}\n")
        f.write(f"Tool: {tool_name}\n")
        f.write(f"Description: {tool_config['description']}\n")
        f.write(f"Exit Code: {report_data['exit_code']}\n")
        f.write("\n" + "="*80 + "\n\n")
        f.write(report_data["output"])
        
        if report_data.get("errors"):
            f.write("\n" + "="*80 + "\n")
            f.write("Summary of Issues:\n")
            for error in report_data["errors"]:
                f.write(f"  - {error}\n")
    
    return report_file


def generate_summary_report(all_reports: List[Dict], timestamp: str):
    """Generate a summary report of all checks."""
    summary_file = LINT_REPORTS_DIR / f"{timestamp}_summary.md"
    
    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("# Code Quality Check Summary\n\n")
        f.write(f"**Generated at:** {datetime.now().isoformat()}\n\n")
        f.write("## Overview\n\n")
        
        total_issues = 0
        failed_tools = []
        
        f.write("| Tool | Status | Issues | Description |\n")
        f.write("|------|--------|--------|-------------|\n")
        
        for report in all_reports:
            tool_name = report["tool"]
            tool_config = TOOLS[tool_name]
            exit_code = report["exit_code"]
            errors = report.get("errors", [])
            error_count = len(errors)
            
            status = "[PASS]" if exit_code == 0 else "[FAIL]"
            if exit_code != 0:
                failed_tools.append(tool_name)
                total_issues += error_count
            
            f.write(f"| {tool_name} | {status} | {error_count} | {tool_config['description']} |\n")
        
        f.write("\n## Details\n\n")
        
        for report in all_reports:
            tool_name = report["tool"]
            f.write(f"### {tool_name.upper()}\n\n")
            
            if report.get("errors"):
                f.write("**Issues found:**\n")
                for error in report["errors"]:
                    f.write(f"- {error}\n")
                f.write("\n")
            else:
                f.write("[OK] No issues found\n\n")
            
            report_file = LINT_REPORTS_DIR / f"{timestamp}_{TOOLS[tool_name]['report_file']}"
            f.write(f"**Full report:** `{report_file.name}`\n\n")
        
        f.write("\n## Recommendations\n\n")
        
        if failed_tools:
            f.write("### Immediate Actions:\n\n")
            if "bandit" in failed_tools:
                f.write("- **[P0] Security**: Review and fix high-severity Bandit findings\n")
            if "black" in failed_tools:
                f.write("- **[P2] Formatting**: Run `python scripts/run_all_checks.py --format` to auto-fix\n")
            if "flake8" in failed_tools:
                f.write("- **[P2] Style**: Fix Flake8 style issues\n")
            if "pylint" in failed_tools:
                f.write("- **[P1] Quality**: Review Pylint findings and improve code quality\n")
            if "mypy" in failed_tools:
                f.write("- **[P1] Types**: Add type annotations to fix MyPy errors\n")
            if "pytest" in failed_tools:
                f.write("- **[P0] Tests**: Fix failing tests\n")
        else:
            f.write("[OK] All checks passed! Code quality is good.\n")
        
        f.write(f"\n**Total Issues:** {total_issues}\n")
        f.write(f"**Failed Tools:** {len(failed_tools)}\n")
    
    return summary_file


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Run all code quality checks")
    parser.add_argument(
        "--format",
        action="store_true",
        help="Auto-fix formatting issues with Black",
    )
    parser.add_argument(
        "--tool",
        choices=list(TOOLS.keys()),
        help="Run only a specific tool",
    )
    parser.add_argument(
        "--no-tests",
        action="store_true",
        help="Skip running pytest",
    )
    parser.add_argument(
        "--ci",
        action="store_true",
        help="CI mode: non-interactive, exit with error code if any check fails",
    )
    
    args = parser.parse_args()
    
    # Handle conflicting arguments: --tool pytest and --no-tests
    if args.tool == "pytest" and args.no_tests:
        print("[ERROR] Cannot use --tool pytest with --no-tests. These options are conflicting.")
        sys.exit(1)
    
    # Check conda environment
    in_conda, conda_python = ensure_conda_environment()
    if not in_conda:
        if conda_python:
            # Use conda environment Python for all operations
            print(f"[INFO] Switching to conda Python: {conda_python}")
            sys.executable = conda_python
            # Update PATH to include conda environment's Scripts directory
            conda_scripts = str(Path(conda_python).parent / "Scripts")
            if conda_scripts not in os.environ.get("PATH", ""):
                os.environ["PATH"] = conda_scripts + os.pathsep + os.environ.get("PATH", "")
        else:
            print("\n[ERROR] This script must be run in a conda environment.")
            print("[INFO] Please activate a conda environment and try again.")
            sys.exit(1)
    
    # Check if tools are available
    missing_tools = []
    tools_to_run = [args.tool] if args.tool else list(TOOLS.keys())
    # Remove pytest from tools_to_run if --no-tests is set
    # This ensures pytest is completely excluded when --no-tests is used
    if args.no_tests:
        if "pytest" in tools_to_run:
            tools_to_run.remove("pytest")
        # Also check if user explicitly requested pytest with --no-tests
        if args.tool == "pytest":
            print("[ERROR] Cannot use --tool pytest with --no-tests. These options are conflicting.")
            sys.exit(1)
    
    for tool_name in tools_to_run:
        if tool_name == "pytest":
            continue  # pytest availability checked separately
        if not check_tool_available(tool_name):
            missing_tools.append(tool_name)
    
    if missing_tools:
        print(f"[WARNING] The following tools are not available: {', '.join(missing_tools)}")
        # Check if running in interactive mode
        is_interactive = sys.stdin.isatty() and not args.ci
        
        if args.ci or not is_interactive:
            print("[INFO] Non-interactive mode: Installing missing tools automatically...")
            install_cmd = [sys.executable, "-m", "pip", "install", "-q"] + missing_tools
            result = subprocess.run(install_cmd, check=False)
            if result.returncode != 0:
                print(f"[ERROR] Failed to install missing tools. Please install manually:")
                print(f"    {sys.executable} -m pip install {' '.join(missing_tools)}")
                sys.exit(1)
            print("[INFO] Missing tools installed successfully.")
        else:
            try:
                response = input("Install missing tools? (y/n): ")
                if response.lower() == 'y':
                    subprocess.run([sys.executable, "-m", "pip", "install", "-q"] + missing_tools, check=False)
            except (EOFError, KeyboardInterrupt):
                print("\n[INFO] Non-interactive mode detected. Installing missing tools automatically...")
                subprocess.run([sys.executable, "-m", "pip", "install", "-q"] + missing_tools, check=False)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    all_reports = []
    
    print("="*80)
    print("Code Quality Check")
    print("="*80)
    print(f"Timestamp: {timestamp}\n")
    
    # Run checks
    if "black" in tools_to_run:
        print("Running Black...")
        report = run_black_check(CHECK_DIRS, fix=args.format)
        save_report("black", report, timestamp)
        all_reports.append(report)
        print(f"  {'[OK]' if report['exit_code'] == 0 else '[FAIL]'} Black check completed")
    
    if "flake8" in tools_to_run:
        print("Running Flake8...")
        report = run_flake8_check(CHECK_DIRS)
        save_report("flake8", report, timestamp)
        all_reports.append(report)
        print(f"  {'[OK]' if report['exit_code'] == 0 else '[FAIL]'} Flake8 check completed")
    
    if "pylint" in tools_to_run:
        print("Running Pylint...")
        report = run_pylint_check(CHECK_DIRS)
        save_report("pylint", report, timestamp)
        all_reports.append(report)
        print(f"  {'[OK]' if report['exit_code'] == 0 else '[FAIL]'} Pylint check completed")
    
    if "mypy" in tools_to_run:
        print("Running MyPy...")
        report = run_mypy_check(CHECK_DIRS)
        save_report("mypy", report, timestamp)
        all_reports.append(report)
        print(f"  {'[OK]' if report['exit_code'] == 0 else '[FAIL]'} MyPy check completed")
    
    if "bandit" in tools_to_run:
        print("Running Bandit...")
        report = run_bandit_check(CHECK_DIRS)
        save_report("bandit", report, timestamp)
        all_reports.append(report)
        print(f"  {'[OK]' if report['exit_code'] == 0 else '[FAIL]'} Bandit check completed")
    
    # Run pytest only if it's in tools_to_run and --no-tests is not set
    # Double-check to ensure pytest is not run when --no-tests is specified
    if "pytest" in tools_to_run and not args.no_tests:
        print("Running Pytest...")
        report = run_pytest()
        save_report("pytest", report, timestamp)
        all_reports.append(report)
        print(f"  {'[OK]' if report['exit_code'] == 0 else '[FAIL]'} Pytest completed")
    
    # Generate summary
    print("\nGenerating summary report...")
    summary_file = generate_summary_report(all_reports, timestamp)
    print(f"  [OK] Summary saved to: {summary_file}")
    
    # Print summary to console
    print("\n" + "="*80)
    print("Summary")
    print("="*80)
    for report in all_reports:
        status = "[PASS]" if report['exit_code'] == 0 else "[FAIL]"
        print(f"{report['tool']:10s} {status}")
    
    # Exit with appropriate code for CI
    if args.ci:
        failed_count = sum(1 for r in all_reports if r['exit_code'] != 0)
        sys.exit(failed_count)
    else:
        failed_count = sum(1 for r in all_reports if r['exit_code'] != 0)
        if failed_count > 0:
            print(f"\n[WARNING] {failed_count} check(s) failed. See reports in {LINT_REPORTS_DIR}")
        sys.exit(failed_count)


if __name__ == "__main__":
    main()

