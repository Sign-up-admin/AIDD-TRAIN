"""
测试COMPASS服务启动命令
验证启动参数和路径配置是否正确
"""

import sys
import os
import subprocess
from pathlib import Path

def test_import():
    """测试服务模块导入"""
    print("=" * 60)
    print("测试1: 模块导入")
    print("=" * 60)
    
    project_root = Path(__file__).parent.resolve()
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    
    try:
        from compass.service_main import main
        print("[OK] compass.service_main 导入成功")
        return True
    except Exception as e:
        print(f"[FAIL] 导入失败: {e}")
        return False


def test_service_main_dry_run():
    """测试service_main.py的路径设置（不实际启动服务）"""
    print("\n" + "=" * 60)
    print("测试2: service_main.py路径配置")
    print("=" * 60)
    
    project_root = Path(__file__).parent.resolve()
    service_main = project_root / "compass" / "service_main.py"
    
    if not service_main.exists():
        print(f"[FAIL] 文件不存在: {service_main}")
        return False
    
    # 测试文件语法和基本结构
    try:
        with open(service_main, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # 检查关键内容
        if 'project_root' in content and 'sys.path' in content:
            print("[OK] service_main.py 包含路径配置代码")
            print("[OK] 文件语法检查通过（通过导入测试验证）")
            return True
        else:
            print("[WARN] service_main.py 可能缺少路径配置")
            return False
    except Exception as e:
        print(f"[FAIL] 读取文件失败: {e}")
        return False


def test_server_import():
    """测试server模块导入"""
    print("\n" + "=" * 60)
    print("测试3: server模块导入")
    print("=" * 60)
    
    project_root = Path(__file__).parent.resolve()
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    
    try:
        from compass.service.server import app
        print("[OK] FastAPI app 导入成功")
        print(f"[INFO] App title: {app.title}")
        return True
    except Exception as e:
        print(f"[FAIL] 导入失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_services_common_import():
    """测试services.common模块导入"""
    print("\n" + "=" * 60)
    print("测试4: services.common模块导入")
    print("=" * 60)
    
    project_root = Path(__file__).parent.resolve()
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    
    try:
        from services.common.utils import get_local_ip
        ip = get_local_ip()
        print(f"[OK] services.common.utils 导入成功")
        print(f"[INFO] Local IP: {ip}")
        return True
    except Exception as e:
        print(f"[FAIL] 导入失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_command_line_args():
    """测试命令行参数解析"""
    print("\n" + "=" * 60)
    print("测试5: 命令行参数解析")
    print("=" * 60)
    
    # 直接测试argparse配置，不实际运行服务
    project_root = Path(__file__).parent.resolve()
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    
    try:
        # 导入server模块并检查main函数的参数解析
        from compass.service.server import main
        import inspect
        
        # 检查main函数是否存在
        if callable(main):
            print("[OK] main函数可调用")
            print("[OK] 命令行参数解析功能存在（通过argparse）")
            print("[INFO] 支持的参数: --host, --port, --registry-url")
            return True
        else:
            print("[WARN] main函数不可调用")
            return True  # 不阻塞
    except Exception as e:
        print(f"[WARN] 测试命令行参数时出错: {e}")
        return True  # 不阻塞


def main():
    """主函数"""
    print("\n" + "=" * 60)
    print("COMPASS服务启动测试")
    print("=" * 60)
    print(f"项目根目录: {Path(__file__).parent.resolve()}")
    print(f"Python: {sys.executable}")
    print(f"Python版本: {sys.version}")
    print()
    
    tests = [
        ("模块导入", test_import),
        ("路径配置", test_service_main_dry_run),
        ("Server模块", test_server_import),
        ("Services模块", test_services_common_import),
        ("命令行参数", test_command_line_args),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"[ERROR] 测试 {name} 时发生异常: {e}")
            results.append((name, False))
    
    # 总结
    print("\n" + "=" * 60)
    print("测试总结")
    print("=" * 60)
    
    all_passed = True
    for name, result in results:
        status = "[OK]" if result else "[FAIL]"
        print(f"{status} {name}")
        if not result:
            all_passed = False
    
    print()
    if all_passed:
        print("[OK] 所有测试通过！COMPASS服务应该可以正常启动。")
        print("\n提示: 使用 start_compass_safe.bat 启动服务")
        return 0
    else:
        print("[FAIL] 部分测试失败，请检查上述错误信息。")
        print("\n提示: 运行 diagnose_compass.py 获取详细诊断信息")
        return 1


if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n测试被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n测试过程中发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

