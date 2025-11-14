#!/usr/bin/env python3
"""
前端代码质量检查脚本
整合 ESLint、Stylelint、Prettier、HTMLHint 等工具
"""

import json
import os
import subprocess
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict


@dataclass
class LintResult:
    """检查结果"""
    tool: str
    file: str
    status: str  # 'success', 'error', 'warning'
    message: str
    line: Optional[int] = None
    column: Optional[int] = None


class FrontendLinter:
    """前端代码检查器"""

    def __init__(self, output_dir: str = "lint_reports/frontend"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results: List[LintResult] = []
        self.temp_code_dir = Path("temp_frontend_code")

    def check_node_installed(self) -> bool:
        """检查 Node.js 是否已安装"""
        try:
            result = subprocess.run(
                ["node", "--version"], capture_output=True, text=True, timeout=5
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def check_npm_packages(self) -> bool:
        """检查 npm 包是否已安装"""
        node_modules = Path("node_modules")
        return node_modules.exists() and any(node_modules.iterdir())

    def install_npm_packages(self) -> bool:
        """安装 npm 包"""
        print("正在安装 npm 依赖包...")
        try:
            result = subprocess.run(
                ["npm", "install"],
                capture_output=True,
                text=True,
                timeout=300,
            )
            if result.returncode == 0:
                print("✓ npm 包安装成功")
                return True
            else:
                print(f"✗ npm 包安装失败: {result.stderr}")
                return False
        except Exception as e:
            print(f"✗ 安装 npm 包时出错: {e}")
            return False

    def run_eslint(self, files: List[str]) -> List[LintResult]:
        """运行 ESLint 检查 JavaScript"""
        results = []
        if not files:
            return results

        print("\n运行 ESLint 检查 JavaScript...")
        try:
            # 创建临时配置文件，只检查提取的文件
            cmd = ["npx", "eslint", "--format", "json"] + files
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=60, cwd=Path.cwd()
            )

            if result.returncode == 0:
                print("✓ ESLint 检查通过")
                results.append(
                    LintResult(
                        tool="ESLint",
                        file="all",
                        status="success",
                        message="JavaScript 代码检查通过",
                    )
                )
            else:
                # 解析 ESLint JSON 输出
                try:
                    eslint_results = json.loads(result.stdout)
                    for file_result in eslint_results:
                        file_path = file_result.get("filePath", "unknown")
                        for message in file_result.get("messages", []):
                            results.append(
                                LintResult(
                                    tool="ESLint",
                                    file=file_path,
                                    status="error" if message["severity"] == 2 else "warning",
                                    message=message["message"],
                                    line=message.get("line"),
                                    column=message.get("column"),
                                )
                            )
                except json.JSONDecodeError:
                    results.append(
                        LintResult(
                            tool="ESLint",
                            file="unknown",
                            status="error",
                            message=result.stderr or result.stdout,
                        )
                    )
        except subprocess.TimeoutExpired:
            results.append(
                LintResult(
                    tool="ESLint", file="all", status="error", message="检查超时"
                )
            )
        except Exception as e:
            results.append(
                LintResult(
                    tool="ESLint", file="all", status="error", message=f"运行错误: {e}"
                )
            )

        return results

    def run_stylelint(self, files: List[str]) -> List[LintResult]:
        """运行 Stylelint 检查 CSS"""
        results = []
        if not files:
            return results

        print("\n运行 Stylelint 检查 CSS...")
        try:
            cmd = ["npx", "stylelint", "--formatter", "json"] + files
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=60, cwd=Path.cwd()
            )

            if result.returncode == 0:
                print("✓ Stylelint 检查通过")
                results.append(
                    LintResult(
                        tool="Stylelint",
                        file="all",
                        status="success",
                        message="CSS 代码检查通过",
                    )
                )
            else:
                try:
                    stylelint_results = json.loads(result.stdout)
                    for file_result in stylelint_results:
                        file_path = file_result.get("source", "unknown")
                        for warning in file_result.get("warnings", []):
                            results.append(
                                LintResult(
                                    tool="Stylelint",
                                    file=file_path,
                                    status="warning",
                                    message=warning["text"],
                                    line=warning.get("line"),
                                    column=warning.get("column"),
                                )
                            )
                except json.JSONDecodeError:
                    results.append(
                        LintResult(
                            tool="Stylelint",
                            file="unknown",
                            status="error",
                            message=result.stderr or result.stdout,
                        )
                    )
        except Exception as e:
            results.append(
                LintResult(
                    tool="Stylelint", file="all", status="error", message=f"运行错误: {e}"
                )
            )

        return results

    def run_htmlhint(self, files: List[str]) -> List[LintResult]:
        """运行 HTMLHint 检查 HTML"""
        results = []
        if not files:
            return results

        print("\n运行 HTMLHint 检查 HTML...")
        try:
            # HTMLHint 不支持 JSON 输出，使用自定义格式
            for file_path in files:
                cmd = ["npx", "htmlhint", file_path]
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=30, cwd=Path.cwd()
                )

                if result.returncode == 0:
                    continue
                else:
                    # 解析 HTMLHint 输出
                    for line in result.stdout.split("\n"):
                        if ":" in line and "line" in line.lower():
                            # 尝试解析行号
                            parts = line.split(":")
                            if len(parts) >= 3:
                                line_num = None
                                try:
                                    line_num = int(parts[1].strip())
                                except ValueError:
                                    pass

                                results.append(
                                    LintResult(
                                        tool="HTMLHint",
                                        file=file_path,
                                        status="warning",
                                        message=line.strip(),
                                        line=line_num,
                                    )
                                )
                            else:
                                results.append(
                                    LintResult(
                                        tool="HTMLHint",
                                        file=file_path,
                                        status="warning",
                                        message=line.strip(),
                                    )
                                )

            if not any(r.tool == "HTMLHint" for r in results):
                results.append(
                    LintResult(
                        tool="HTMLHint",
                        file="all",
                        status="success",
                        message="HTML 代码检查通过",
                    )
                )
        except Exception as e:
            results.append(
                LintResult(
                    tool="HTMLHint", file="all", status="error", message=f"运行错误: {e}"
                )
            )

        return results

    def run_prettier_check(self, files: List[str]) -> List[LintResult]:
        """运行 Prettier 格式化检查"""
        results = []
        if not files:
            return results

        print("\n运行 Prettier 格式化检查...")
        try:
            # 检查所有支持的文件类型
            all_files = []
            for ext in [".js", ".css", ".html", ".json"]:
                all_files.extend([f for f in files if f.endswith(ext)])

            if not all_files:
                return results

            cmd = ["npx", "prettier", "--check"] + all_files
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=60, cwd=Path.cwd()
            )

            if result.returncode == 0:
                print("✓ Prettier 格式化检查通过")
                results.append(
                    LintResult(
                        tool="Prettier",
                        file="all",
                        status="success",
                        message="代码格式符合规范",
                    )
                )
            else:
                # Prettier 会列出需要格式化的文件
                for line in result.stdout.split("\n"):
                    if line.strip() and not line.startswith("Checking"):
                        file_path = line.strip()
                        if file_path:
                            results.append(
                                LintResult(
                                    tool="Prettier",
                                    file=file_path,
                                    status="warning",
                                    message="代码格式需要调整，运行 'npm run format' 可自动修复",
                                )
                            )
        except Exception as e:
            results.append(
                LintResult(
                    tool="Prettier", file="all", status="error", message=f"运行错误: {e}"
                )
            )

        return results

    def run_all_checks(self, extract_dir: Optional[Path] = None) -> List[LintResult]:
        """运行所有检查"""
        if extract_dir is None:
            extract_dir = self.temp_code_dir

        if not extract_dir.exists():
            print(f"错误: 临时代码目录不存在: {extract_dir}")
            print("请先运行 extract_frontend_code.py 提取代码")
            return []

        # 检查 Node.js 环境
        if not self.check_node_installed():
            print("错误: 未检测到 Node.js，请先安装 Node.js")
            return []

        if not self.check_npm_packages():
            if not self.install_npm_packages():
                print("错误: 无法安装 npm 依赖包")
                return []

        # 收集文件
        js_files = list(extract_dir.glob("*.js"))
        css_files = list(extract_dir.glob("*.css"))
        html_files = list(extract_dir.glob("*.html"))

        all_results = []

        # 运行各项检查
        if js_files:
            all_results.extend(self.run_eslint([str(f) for f in js_files]))

        if css_files:
            all_results.extend(self.run_stylelint([str(f) for f in css_files]))

        if html_files:
            all_results.extend(self.run_htmlhint([str(f) for f in html_files]))

        # Prettier 检查所有文件
        all_files = [str(f) for f in list(js_files) + list(css_files) + list(html_files)]
        if all_files:
            all_results.extend(self.run_prettier_check(all_files))

        self.results = all_results
        return all_results

    def generate_report(self) -> str:
        """生成检查报告"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # 统计信息
        total = len(self.results)
        errors = len([r for r in self.results if r.status == "error"])
        warnings = len([r for r in self.results if r.status == "warning"])
        success = len([r for r in self.results if r.status == "success"])

        # 按工具分组
        by_tool: Dict[str, List[LintResult]] = {}
        for result in self.results:
            if result.tool not in by_tool:
                by_tool[result.tool] = []
            by_tool[result.tool].append(result)

        # 生成 Markdown 报告
        report_lines = [
            "# 前端代码质量检查报告",
            "",
            f"**生成时间**: {timestamp}",
            "",
            "## 检查摘要",
            "",
            f"- 总计: {total} 项",
            f"- 错误: {errors} 项",
            f"- 警告: {warnings} 项",
            f"- 通过: {success} 项",
            "",
        ]

        # 按工具详细报告
        for tool, results in by_tool.items():
            report_lines.extend([
                f"## {tool}",
                "",
            ])

            tool_errors = [r for r in results if r.status == "error"]
            tool_warnings = [r for r in results if r.status == "warning"]
            tool_success = [r for r in results if r.status == "success"]

            if tool_errors or tool_warnings:
                report_lines.append(f"### 问题列表 ({len(tool_errors)} 错误, {len(tool_warnings)} 警告)")
                report_lines.append("")

                for result in tool_errors + tool_warnings:
                    location = ""
                    if result.line:
                        location = f" (行 {result.line}"
                        if result.column:
                            location += f", 列 {result.column}"
                        location += ")"
                    report_lines.append(
                        f"- **{result.status.upper()}** - `{Path(result.file).name}`{location}: {result.message}"
                    )
                report_lines.append("")
            else:
                report_lines.append("✓ 检查通过，未发现问题")
                report_lines.append("")

        # 保存报告
        report_file = self.output_dir / f"frontend_lint_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
        report_content = "\n".join(report_lines)

        with open(report_file, "w", encoding="utf-8") as f:
            f.write(report_content)

        # 同时保存 JSON 格式
        json_file = self.output_dir / f"frontend_lint_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(json_file, "w", encoding="utf-8") as f:
            json.dump(
                {
                    "timestamp": timestamp,
                    "summary": {
                        "total": total,
                        "errors": errors,
                        "warnings": warnings,
                        "success": success,
                    },
                    "results": [asdict(r) for r in self.results],
                },
                f,
                indent=2,
                ensure_ascii=False,
            )

        print(f"\n报告已保存:")
        print(f"  - Markdown: {report_file}")
        print(f"  - JSON: {json_file}")

        return str(report_file)


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description="前端代码质量检查")
    parser.add_argument(
        "-d", "--dir", type=str, default="temp_frontend_code", help="提取的代码目录"
    )
    parser.add_argument(
        "-o", "--output", type=str, default="lint_reports/frontend", help="报告输出目录"
    )

    args = parser.parse_args()

    linter = FrontendLinter(output_dir=args.output)
    extract_dir = Path(args.dir)

    print("=" * 60)
    print("前端代码质量检查")
    print("=" * 60)

    results = linter.run_all_checks(extract_dir)

    if results:
        report_file = linter.generate_report()

        # 统计摘要
        errors = len([r for r in results if r.status == "error"])
        warnings = len([r for r in results if r.status == "warning"])

        print("\n" + "=" * 60)
        print("检查完成")
        print("=" * 60)
        print(f"总计: {len(results)} 项")
        print(f"错误: {errors} 项")
        print(f"警告: {warnings} 项")
        print(f"通过: {len(results) - errors - warnings} 项")

        if errors > 0:
            sys.exit(1)
    else:
        print("未找到需要检查的文件")


if __name__ == "__main__":
    main()



