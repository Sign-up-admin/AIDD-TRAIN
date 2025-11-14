#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用 Docker 运行前端代码质量检查
避免需要本地安装 Node.js
"""

import subprocess
import sys
import io
from pathlib import Path
from typing import Tuple

# Windows 控制台编码修复
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')


def check_docker_installed() -> bool:
    """检查 Docker 是否已安装"""
    try:
        result = subprocess.run(
            ["docker", "--version"], capture_output=True, text=True, timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def check_docker_compose_installed() -> bool:
    """检查 Docker Compose 是否已安装"""
    try:
        # 尝试 docker compose (v2)
        result = subprocess.run(
            ["docker", "compose", "version"], capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            return True
    except:
        pass

    try:
        # 尝试 docker-compose (v1)
        result = subprocess.run(
            ["docker-compose", "--version"], capture_output=True, text=True, timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def build_docker_image():
    """构建 Docker 镜像"""
    print("正在构建 Docker 镜像...")
    try:
        result = subprocess.run(
            [
                "docker",
                "build",
                "-f",
                "Dockerfile.frontend-lint",
                "-t",
                "aidd-frontend-linter:latest",
                ".",
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            print("[OK] Docker 镜像构建成功")
            return True
        else:
            print(f"[ERROR] Docker 镜像构建失败: {result.stderr}")
            return False
    except Exception as e:
        print(f"[ERROR] 构建 Docker 镜像时出错: {e}")
        return False


def run_docker_command(command: list, work_dir: str = ".") -> Tuple[bool, str]:
    """在 Docker 容器中运行命令"""
    docker_cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{Path.cwd()}:/app",
        "-v",
        "aidd-frontend-lint_node_modules:/app/node_modules",
        "-w",
        "/app",
        "aidd-frontend-linter:latest",
    ] + command

    try:
        result = subprocess.run(
            docker_cmd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace',
            cwd=Path.cwd()
        )
        stdout = result.stdout or ""
        stderr = result.stderr or ""
        return result.returncode == 0, stdout + stderr
    except Exception as e:
        return False, str(e)


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description="使用 Docker 运行前端代码检查")
    parser.add_argument(
        "command",
        nargs="?",
        choices=["lint:js", "lint:css", "lint:html", "lint:all", "format", "format:check"],
        default="lint:all",
        help="要运行的命令",
    )
    parser.add_argument(
        "--build", action="store_true", help="重新构建 Docker 镜像"
    )
    parser.add_argument(
        "--extract-only",
        action="store_true",
        help="仅提取代码，不运行检查",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("Docker 前端代码检查")
    print("=" * 60)

    # 检查 Docker
    if not check_docker_installed():
        print("错误: 未检测到 Docker，请先安装 Docker")
        print("下载地址: https://www.docker.com/get-started")
        sys.exit(1)

    if not check_docker_compose_installed():
        print("警告: 未检测到 Docker Compose，将使用 docker run")

    # 构建镜像
    if args.build or not _image_exists():
        print("\n提示: 如果网络连接 Docker Hub 困难，请先配置镜像加速器")
        print("参考文档: docs/DOCKER_SETUP.md")
        print("\n正在构建镜像...")
        if not build_docker_image():
            print("\n构建失败！可能的解决方案：")
            print("1. 配置 Docker 镜像加速器（推荐）")
            print("   - Windows: Docker Desktop > Settings > Docker Engine")
            print("   - 添加镜像加速器地址")
            print("2. 使用国内镜像源构建：")
            print("   docker build -f Dockerfile.frontend-lint.cn -t aidd-frontend-linter:latest .")
            print("3. 检查网络连接或使用代理")
            sys.exit(1)
    else:
        print("[OK] Docker 镜像已存在，跳过构建")

    # 提取代码（如果需要）
    if not args.extract_only:
        print("\n提取前端代码...")
        extract_result = subprocess.run(
            [
                sys.executable,
                "scripts/extract_frontend_code.py",
                "FLASH_DOCK-main",
                "-o",
                "temp_frontend_code",
            ],
            capture_output=True,
            text=True,
        )
        if extract_result.returncode == 0:
            print("[OK] 代码提取完成")
        else:
            print(f"警告: 代码提取可能失败: {extract_result.stderr}")

    if args.extract_only:
        return

    # 运行检查
    print(f"\n运行命令: {args.command}")
    npm_commands = {
        "lint:js": ["npm", "run", "lint:js"],
        "lint:css": ["npm", "run", "lint:css"],
        "lint:html": ["npm", "run", "lint:html"],
        "lint:all": ["npm", "run", "lint:all"],
        "format": ["npm", "run", "format"],
        "format:check": ["npm", "run", "format:check"],
    }

    command = npm_commands.get(args.command, ["npm", "run", "lint:all"])
    success, output = run_docker_command(command)

    print("\n" + "=" * 60)
    print("输出:")
    print("=" * 60)
    print(output)

    if success:
        print("\n[OK] 检查完成")
    else:
        print("\n[ERROR] 检查失败")
        sys.exit(1)


def _image_exists() -> bool:
    """检查镜像是否存在"""
    try:
        result = subprocess.run(
            ["docker", "images", "-q", "aidd-frontend-linter:latest"],
            capture_output=True,
            text=True,
        )
        return bool(result.stdout.strip())
    except:
        return False


if __name__ == "__main__":
    main()

