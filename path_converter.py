#!/usr/bin/env python3
"""
Windows/WSL2 路径转换工具
用于在 Windows 和 WSL2 之间转换文件路径
"""

import os
import sys
import platform


def is_wsl2():
    """检测是否在 WSL2 环境中运行"""
    try:
        if os.path.exists('/proc/version'):
            with open('/proc/version', 'r') as f:
                version_info = f.read().lower()
                if 'microsoft' in version_info or 'wsl' in version_info:
                    return True
        if os.environ.get('WSL_DISTRO_NAME') or os.environ.get('WSL_INTEROP'):
            return True
    except:
        pass
    return False


def convert_to_wsl_path(windows_path):
    """将 Windows 路径转换为 WSL2 路径"""
    if not windows_path:
        return windows_path
    
    # 如果已经是 Linux 路径格式，直接返回
    if windows_path.startswith('/'):
        return windows_path
    
    # 转换 Windows 路径到 WSL2 路径
    path = os.path.abspath(windows_path)
    
    # 提取驱动器字母
    if len(path) >= 2 and path[1] == ':':
        drive_letter = path[0].lower()
        rest_path = path[2:].replace('\\', '/')
        wsl_path = f'/mnt/{drive_letter}{rest_path}'
        return wsl_path
    
    return path


def convert_to_windows_path(wsl_path):
    """将 WSL2 路径转换为 Windows 路径"""
    if not wsl_path:
        return wsl_path
    
    # 如果是 Windows 路径格式，直接返回
    if ':' in wsl_path and (wsl_path[1] == ':' or wsl_path.startswith('\\\\')):
        return wsl_path
    
    # 转换 WSL2 路径到 Windows 路径
    if wsl_path.startswith('/mnt/'):
        parts = wsl_path[5:].split('/', 1)
        if len(parts) == 2:
            drive_letter = parts[0].upper()
            rest_path = parts[1].replace('/', '\\')
            windows_path = f'{drive_letter}:\\{rest_path}'
            return windows_path
    
    return wsl_path


def main():
    if len(sys.argv) < 2:
        print("用法:")
        print("  python path_converter.py <路径>")
        print("  python path_converter.py --to-wsl <Windows路径>")
        print("  python path_converter.py --to-windows <WSL路径>")
        print("")
        print("示例:")
        print("  python path_converter.py E:\\Qinchaojun\\AIDD-TRAIN")
        print("  python path_converter.py --to-wsl E:\\Qinchaojun\\AIDD-TRAIN")
        print("  python path_converter.py --to-windows /mnt/e/Qinchaojun/AIDD-TRAIN")
        sys.exit(1)
    
    is_wsl = is_wsl2()
    is_windows = platform.system() == 'Windows'
    
    if sys.argv[1] == '--to-wsl':
        if len(sys.argv) < 3:
            print("错误: 需要提供 Windows 路径")
            sys.exit(1)
        windows_path = sys.argv[2]
        wsl_path = convert_to_wsl_path(windows_path)
        print(f"Windows 路径: {windows_path}")
        print(f"WSL2 路径:    {wsl_path}")
    elif sys.argv[1] == '--to-windows':
        if len(sys.argv) < 3:
            print("错误: 需要提供 WSL 路径")
            sys.exit(1)
        wsl_path = sys.argv[2]
        windows_path = convert_to_windows_path(wsl_path)
        print(f"WSL2 路径:    {wsl_path}")
        print(f"Windows 路径: {windows_path}")
    else:
        # 自动检测并转换
        path = sys.argv[1]
        
        if is_wsl:
            # 在 WSL2 中，如果是 Windows 路径格式，转换为 WSL 路径
            if ':' in path and path[1] == ':':
                converted = convert_to_wsl_path(path)
                print(f"原始路径: {path}")
                print(f"WSL2 路径: {converted}")
            else:
                # 已经是 WSL 路径，转换为 Windows 路径
                converted = convert_to_windows_path(path)
                print(f"原始路径: {path}")
                print(f"Windows 路径: {converted}")
        elif is_windows:
            # 在 Windows 中，如果是 WSL 路径格式，转换为 Windows 路径
            if path.startswith('/mnt/'):
                converted = convert_to_windows_path(path)
                print(f"原始路径: {path}")
                print(f"Windows 路径: {converted}")
            else:
                # 已经是 Windows 路径，转换为 WSL 路径
                converted = convert_to_wsl_path(path)
                print(f"原始路径: {path}")
                print(f"WSL2 路径: {converted}")
        else:
            print(f"当前环境: {platform.system()}")
            print(f"路径: {path}")


if __name__ == '__main__':
    main()



