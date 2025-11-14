"""
COMPASS自动化测试运行器
提供统一的测试执行入口和报告生成
"""

import os
import sys
import argparse
import subprocess
import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional


class CompassTestRunner:
    """COMPASS测试运行器"""
    
    def __init__(self, project_root: Optional[Path] = None):
        self.project_root = project_root or Path(__file__).parent
        self.test_dir = self.project_root / "tests"
        self.reports_dir = self.project_root / "reports" / "test_reports"
        self.reports_dir.mkdir(parents=True, exist_ok=True)
    
    def run_tests(self, 
                  test_type: str = "all",
                  verbose: bool = True,
                  coverage: bool = True,
                  markers: Optional[List[str]] = None,
                  output_format: str = "both") -> Dict:
        """
        运行测试
        
        Args:
            test_type: 测试类型 (all, unit, integration, e2e)
            verbose: 是否显示详细输出
            coverage: 是否生成覆盖率报告
            markers: pytest标记列表
            output_format: 输出格式 (console, json, both)
        
        Returns:
            测试结果字典
        """
        print("\n" + "=" * 70)
        print("COMPASS自动化测试运行器")
        print("=" * 70)
        print(f"项目根目录: {self.project_root}")
        print(f"测试目录: {self.test_dir}")
        print(f"测试类型: {test_type}")
        print(f"开始时间: {datetime.now()}")
        print("=" * 70 + "\n")
        
        # 构建pytest命令
        cmd = ["pytest"]
        
        # 添加测试路径
        if test_type == "all":
            cmd.append(str(self.test_dir))
        elif test_type == "unit":
            cmd.extend([str(self.test_dir), "-m", "unit"])
        elif test_type == "integration":
            cmd.extend([str(self.test_dir), "-m", "integration"])
        elif test_type == "e2e":
            cmd.extend([str(self.test_dir), "-m", "e2e"])
        else:
            # 特定测试文件
            test_file = self.test_dir / f"test_{test_type}.py"
            if test_file.exists():
                cmd.append(str(test_file))
            else:
                print(f"错误: 找不到测试文件 {test_file}")
                return {"success": False, "error": f"测试文件不存在: {test_file}"}
        
        # 添加标记
        if markers:
            for marker in markers:
                cmd.extend(["-m", marker])
        
        # 详细输出
        if verbose:
            cmd.append("-v")
        
        # 覆盖率
        if coverage:
            cmd.extend([
                "--cov=compass",
                "--cov-report=term-missing",
                "--cov-report=html:" + str(self.reports_dir / "htmlcov"),
                "--cov-report=xml:" + str(self.reports_dir / "coverage.xml"),
                "--cov-report=json:" + str(self.reports_dir / "coverage.json")
            ])
        
        # JSON报告
        json_report = self.reports_dir / f"test_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        cmd.extend(["--json-report", "--json-report-file", str(json_report)])
        
        # 其他选项
        cmd.extend([
            "--tb=short",
            "--strict-markers",
            "-ra"  # 显示所有测试结果摘要
        ])
        
        print(f"执行命令: {' '.join(cmd)}\n")
        
        # 运行测试
        try:
            result = subprocess.run(
                cmd,
                cwd=str(self.project_root),
                capture_output=True,
                text=True
            )
            
            # 解析结果
            test_results = {
                "success": result.returncode == 0,
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "command": " ".join(cmd),
                "timestamp": datetime.now().isoformat(),
                "json_report": str(json_report) if json_report.exists() else None
            }
            
            # 读取JSON报告
            if json_report.exists():
                try:
                    with open(json_report, 'r', encoding='utf-8') as f:
                        test_results["json_data"] = json.load(f)
                except Exception as e:
                    print(f"警告: 无法读取JSON报告: {e}")
            
            # 输出结果
            print("\n" + "=" * 70)
            print("测试执行完成")
            print("=" * 70)
            print(f"返回码: {result.returncode}")
            print(f"成功: {'是' if test_results['success'] else '否'}")
            
            if output_format in ["console", "both"]:
                print("\n标准输出:")
                print(result.stdout)
                if result.stderr:
                    print("\n标准错误:")
                    print(result.stderr)
            
            if output_format in ["json", "both"] and test_results.get("json_data"):
                print(f"\nJSON报告: {json_report}")
                summary = test_results["json_data"].get("summary", {})
                print(f"总测试数: {summary.get('total', 0)}")
                print(f"通过: {summary.get('passed', 0)}")
                print(f"失败: {summary.get('failed', 0)}")
                print(f"跳过: {summary.get('skipped', 0)}")
            
            if coverage:
                print(f"\n覆盖率报告:")
                print(f"  HTML: {self.reports_dir / 'htmlcov' / 'index.html'}")
                print(f"  XML: {self.reports_dir / 'coverage.xml'}")
                print(f"  JSON: {self.reports_dir / 'coverage.json'}")
            
            print("=" * 70 + "\n")
            
            return test_results
            
        except Exception as e:
            error_msg = f"测试执行失败: {e}"
            print(f"错误: {error_msg}")
            return {
                "success": False,
                "error": error_msg,
                "timestamp": datetime.now().isoformat()
            }
    
    def run_specific_test(self, test_path: str, verbose: bool = True) -> Dict:
        """运行特定测试"""
        test_file = Path(test_path)
        if not test_file.is_absolute():
            test_file = self.test_dir / test_path
        
        if not test_file.exists():
            return {
                "success": False,
                "error": f"测试文件不存在: {test_file}"
            }
        
        return self.run_tests(
            test_type=str(test_file),
            verbose=verbose,
            coverage=False
        )
    
    def generate_test_summary(self) -> str:
        """生成测试摘要报告"""
        # 查找最新的JSON报告
        json_reports = sorted(self.reports_dir.glob("test_report_*.json"))
        if not json_reports:
            return "没有找到测试报告"
        
        latest_report = json_reports[-1]
        try:
            with open(latest_report, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            summary = data.get("summary", {})
            report = []
            report.append("=" * 70)
            report.append("COMPASS测试摘要报告")
            report.append("=" * 70)
            report.append(f"报告时间: {data.get('created', 'N/A')}")
            report.append(f"总测试数: {summary.get('total', 0)}")
            report.append(f"通过: {summary.get('passed', 0)}")
            report.append(f"失败: {summary.get('failed', 0)}")
            report.append(f"跳过: {summary.get('skipped', 0)}")
            report.append("=" * 70)
            
            return "\n".join(report)
        except Exception as e:
            return f"无法读取报告: {e}"


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="COMPASS自动化测试运行器")
    parser.add_argument(
        "--type",
        choices=["all", "unit", "integration", "e2e"],
        default="all",
        help="测试类型"
    )
    parser.add_argument(
        "--test",
        type=str,
        help="运行特定测试文件"
    )
    parser.add_argument(
        "--no-coverage",
        action="store_true",
        help="不生成覆盖率报告"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="安静模式（不显示详细输出）"
    )
    parser.add_argument(
        "--markers",
        nargs="+",
        help="pytest标记（例如: unit integration）"
    )
    parser.add_argument(
        "--output",
        choices=["console", "json", "both"],
        default="both",
        help="输出格式"
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="显示测试摘要"
    )
    
    args = parser.parse_args()
    
    runner = CompassTestRunner()
    
    if args.summary:
        print(runner.generate_test_summary())
        return
    
    if args.test:
        results = runner.run_specific_test(args.test, verbose=not args.quiet)
    else:
        results = runner.run_tests(
            test_type=args.type,
            verbose=not args.quiet,
            coverage=not args.no_coverage,
            markers=args.markers,
            output_format=args.output
        )
    
    # 返回适当的退出码
    sys.exit(0 if results.get("success") else 1)


if __name__ == "__main__":
    main()

