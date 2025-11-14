#!/usr/bin/env python3
"""
从 Python 文件中提取内嵌的 HTML/CSS/JavaScript 代码
用于前端代码质量检查
"""

import re
import os
import json
from pathlib import Path
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class ExtractedCode:
    """提取的代码片段"""
    source_file: str
    line_start: int
    line_end: int
    code_type: str  # 'html', 'css', 'javascript'
    content: str
    context: str  # 函数名或上下文


class FrontendCodeExtractor:
    """前端代码提取器"""

    def __init__(self, output_dir: str = "temp_frontend_code"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.extracted_files: List[ExtractedCode] = []

    def extract_from_python_file(self, file_path: Path) -> List[ExtractedCode]:
        """从 Python 文件中提取前端代码"""
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
            lines = content.split("\n")

        extracted = []

        # 匹配 HTML 字符串（包含 <html>, <!DOCTYPE html> 等）
        html_pattern = r'return\s+f?"""(.*?)"""'
        html_matches = list(re.finditer(html_pattern, content, re.DOTALL))

        for match in html_matches:
            html_content = match.group(1)
            if "<html" in html_content or "<!DOCTYPE" in html_content:
                # 计算行号
                line_start = content[: match.start()].count("\n") + 1
                line_end = content[: match.end()].count("\n") + 1

                # 提取函数名作为上下文
                context_match = re.search(
                    r"def\s+(\w+).*?return\s+f?\"\"\"", content[: match.start()], re.DOTALL
                )
                context = context_match.group(1) if context_match else "unknown"

                # 分离 HTML、CSS 和 JavaScript
                html_part, css_part, js_part = self._separate_html_css_js(html_content)

                if html_part:
                    extracted.append(
                        ExtractedCode(
                            source_file=str(file_path),
                            line_start=line_start,
                            line_end=line_end,
                            code_type="html",
                            content=html_part,
                            context=context,
                        )
                    )

                if css_part:
                    extracted.append(
                        ExtractedCode(
                            source_file=str(file_path),
                            line_start=line_start,
                            line_end=line_end,
                            code_type="css",
                            content=css_part,
                            context=context,
                        )
                    )

                if js_part:
                    extracted.append(
                        ExtractedCode(
                            source_file=str(file_path),
                            line_start=line_start,
                            line_end=line_end,
                            code_type="javascript",
                            content=js_part,
                            context=context,
                        )
                    )

        return extracted

    def _separate_html_css_js(self, html_content: str) -> Tuple[str, str, str]:
        """分离 HTML、CSS 和 JavaScript 代码"""
        html_part = ""
        css_part = ""
        js_part = ""

        # 提取 <style> 标签内容
        style_pattern = r"<style[^>]*>(.*?)</style>"
        style_matches = re.finditer(style_pattern, html_content, re.DOTALL | re.IGNORECASE)
        css_parts = [match.group(1) for match in style_matches]
        if css_parts:
            css_part = "\n\n".join(css_parts)

        # 提取 <script> 标签内容
        script_pattern = r"<script[^>]*>(.*?)</script>"
        script_matches = re.finditer(script_pattern, html_content, re.DOTALL | re.IGNORECASE)
        js_parts = [match.group(1) for match in script_matches]
        if js_parts:
            js_part = "\n\n".join(js_parts)

        # HTML 部分（移除 style 和 script 标签）
        html_part = html_content
        html_part = re.sub(style_pattern, "", html_part, flags=re.DOTALL | re.IGNORECASE)
        html_part = re.sub(script_pattern, "", html_part, flags=re.DOTALL | re.IGNORECASE)

        return html_part.strip(), css_part.strip(), js_part.strip()

    def save_extracted_code(self, extracted: List[ExtractedCode]) -> Dict[str, str]:
        """保存提取的代码到临时文件"""
        saved_files = {}

        for idx, code in enumerate(extracted):
            # 生成文件名
            source_name = Path(code.source_file).stem
            filename = f"{source_name}_{code.context}_{idx+1}.{self._get_extension(code.code_type)}"
            file_path = self.output_dir / filename

            # 保存文件
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(code.content)

            saved_files[str(file_path)] = {
                "source": code.source_file,
                "line_start": code.line_start,
                "line_end": code.line_end,
                "type": code.code_type,
                "context": code.context,
            }

            self.extracted_files.append(code)

        return saved_files

    def _get_extension(self, code_type: str) -> str:
        """获取文件扩展名"""
        extensions = {
            "html": "html",
            "css": "css",
            "javascript": "js",
        }
        return extensions.get(code_type, "txt")

    def extract_from_directory(self, directory: Path, pattern: str = "*.py") -> Dict[str, str]:
        """从目录中提取所有 Python 文件的前端代码"""
        all_extracted = []
        python_files = list(directory.rglob(pattern))

        for py_file in python_files:
            # 跳过 __pycache__ 和虚拟环境
            if "__pycache__" in str(py_file) or "venv" in str(py_file):
                continue

            try:
                extracted = self.extract_from_python_file(py_file)
                all_extracted.extend(extracted)
            except Exception as e:
                print(f"警告: 无法处理文件 {py_file}: {e}")

        return self.save_extracted_code(all_extracted)

    def cleanup(self):
        """清理临时文件"""
        if self.output_dir.exists():
            import shutil

            shutil.rmtree(self.output_dir)
            print(f"已清理临时目录: {self.output_dir}")


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description="从 Python 文件中提取前端代码")
    parser.add_argument(
        "directory", type=str, help="要扫描的目录", default="FLASH_DOCK-main", nargs="?"
    )
    parser.add_argument(
        "-o", "--output", type=str, default="temp_frontend_code", help="输出目录"
    )
    parser.add_argument(
        "-c", "--cleanup", action="store_true", help="检查完成后清理临时文件"
    )
    parser.add_argument(
        "-j", "--json", type=str, help="保存文件映射到 JSON 文件"
    )

    args = parser.parse_args()

    extractor = FrontendCodeExtractor(output_dir=args.output)
    directory = Path(args.directory) if args.directory else Path("FLASH_DOCK-main")

    print(f"正在从 {directory} 提取前端代码...")
    file_map = extractor.extract_from_directory(directory)

    print(f"\n提取完成！共提取 {len(extractor.extracted_files)} 个代码片段")
    print(f"临时文件保存在: {extractor.output_dir}")

    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(file_map, f, indent=2, ensure_ascii=False)
        print(f"文件映射已保存到: {args.json}")

    if args.cleanup:
        extractor.cleanup()


if __name__ == "__main__":
    main()



