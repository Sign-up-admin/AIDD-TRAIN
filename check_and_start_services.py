"""
检查并启动COMPASS服务
"""

import sys
import os
import time
import requests
import subprocess
from pathlib import Path


def check_service(url, name):
    """检查服务是否运行"""
    try:
        r = requests.get(url, timeout=2)
        return r.status_code == 200
    except:
        return False


def main():
    print("=" * 60)
    print("COMPASS服务状态检查")
    print("=" * 60)
    print()

    # 检查注册中心
    registry_ok = check_service("http://localhost:8500/health", "服务注册中心")
    print(f"服务注册中心 (8500): {'[OK] 运行中' if registry_ok else '[FAIL] 未运行'}")

    # 检查COMPASS服务
    compass_ok = check_service("http://localhost:8080/health", "COMPASS服务")
    print(f"COMPASS服务 (8080): {'[OK] 运行中' if compass_ok else '[FAIL] 未运行'}")
    print()

    if not registry_ok:
        print("错误: 服务注册中心未运行，请先启动服务注册中心")
        print("运行: python services/registry/server.py --host 0.0.0.0 --port 8500")
        return 1

    if not compass_ok:
        print("=" * 60)
        print("正在启动COMPASS服务...")
        print("=" * 60)

        project_root = Path(__file__).parent.resolve()  # 使用resolve()确保绝对路径
        python_cmd = sys.executable

        # 设置环境变量 - 确保路径是绝对路径且有效
        env = os.environ.copy()
        project_root_str = str(project_root)

        # 清理PYTHONPATH，只保留有效路径
        existing_paths = env.get("PYTHONPATH", "").split(os.pathsep)
        valid_paths = []
        for path in existing_paths:
            if path and os.path.exists(path):
                valid_paths.append(os.path.abspath(path))

        # 添加项目根目录
        if project_root_str not in valid_paths:
            valid_paths.insert(0, project_root_str)

        # 设置PYTHONPATH
        env["PYTHONPATH"] = os.pathsep.join(valid_paths)

        # 启动服务
        cmd = [
            python_cmd,
            str(project_root / "compass" / "service_main.py"),
            "--host",
            "0.0.0.0",
            "--port",
            "8080",
            "--registry-url",
            "http://localhost:8500",
        ]

        print(f"命令: {' '.join(cmd)}")
        print()
        print("提示: 服务将在新窗口中启动")
        print("请检查新窗口中的服务状态")
        print()

        # 在Windows上使用start命令在新窗口启动
        if sys.platform == "win32":
            # 确保工作目录是绝对路径
            os.chdir(str(project_root.resolve()))
            # 确保所有路径都是绝对路径
            cmd = [str(Path(c).resolve()) if Path(c).exists() else c for c in cmd]
            subprocess.Popen(
                cmd,
                creationflags=subprocess.CREATE_NEW_CONSOLE,
                env=env,
                cwd=str(project_root.resolve()),
            )
            print("服务已在新窗口中启动")
            print("等待5秒后检查服务状态...")
            time.sleep(5)

            if check_service("http://localhost:8080/health", "COMPASS服务"):
                print("[OK] COMPASS服务已成功启动!")
                return 0
            else:
                print("[WARNING] 服务可能还在启动中，请检查新窗口中的状态")
                return 0
        else:
            # Linux/Mac: 直接启动
            subprocess.Popen(cmd, env=env)
            print("服务已在后台启动")
            return 0
    else:
        print("=" * 60)
        print("[OK] 所有服务都在运行!")
        print("=" * 60)
        print()
        print("服务地址:")
        print("  - 服务注册中心: http://localhost:8500")
        print("  - COMPASS服务: http://localhost:8080")
        print("  - API文档: http://localhost:8080/docs")
        return 0


if __name__ == "__main__":
    import os

    sys.exit(main())
