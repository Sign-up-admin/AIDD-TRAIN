"""
Java环境诊断脚本
用于检查口袋预测功能所需的Java环境配置
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Windows编码兼容性处理
if sys.platform == "win32":
    import io
    # 设置标准输出为UTF-8编码
    if sys.stdout.encoding != 'utf-8':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    if sys.stderr.encoding != 'utf-8':
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')


def check_java_environment():
    """检查Java环境"""
    # 获取脚本所在目录作为基准路径
    script_dir = Path(__file__).parent.absolute()
    
    print("=" * 60)
    print("Java环境诊断")
    print("=" * 60)
    print(f"工作目录: {script_dir}")
    print()
    
    all_ok = True
    
    # 检查java命令
    java_path = shutil.which("java")
    if java_path:
        print(f"✓ Java可执行文件找到: {java_path}")
        try:
            result = subprocess.run(
                ["java", "-version"],
                capture_output=True,
                text=True,
                stderr=subprocess.STDOUT,
                timeout=5
            )
            print(f"✓ Java版本信息:")
            # Java版本信息通常在stderr中
            version_output = result.stdout or result.stderr
            print(version_output)
            
            # 检查版本号
            if "version" in version_output.lower():
                # 尝试提取版本号
                import re
                version_match = re.search(r'version\s+"?(\d+)', version_output)
                if version_match:
                    version = int(version_match.group(1))
                    if 17 <= version <= 23:
                        print(f"✓ Java版本 {version} 符合要求（17-23）")
                    else:
                        print(f"✗ Java版本 {version} 不符合要求（需要17-23）")
                        all_ok = False
        except subprocess.TimeoutExpired:
            print("✗ Java命令执行超时")
            all_ok = False
        except Exception as e:
            print(f"✗ 无法执行java -version: {e}")
            all_ok = False
    else:
        print("✗ Java可执行文件未找到")
        all_ok = False
    
    print()
    
    # 检查JAVA_HOME
    java_home = os.environ.get("JAVA_HOME")
    if java_home:
        print(f"✓ JAVA_HOME: {java_home}")
        java_exe = Path(java_home) / "bin" / "java.exe"
        if java_exe.exists():
            print(f"✓ JAVA_HOME/bin/java.exe 存在")
        else:
            print(f"✗ JAVA_HOME/bin/java.exe 不存在")
            all_ok = False
    else:
        print("✗ JAVA_HOME环境变量未设置")
        if not java_path:
            all_ok = False
    
    print()
    
    # 检查P2Rank
    p2rank_home = script_dir / "others" / "p2rank_2.5"
    if p2rank_home.exists():
        print(f"✓ P2Rank目录存在: {p2rank_home.absolute()}")
        prank_bat = p2rank_home / "prank.bat"
        prank_script = p2rank_home / "prank"
        p2rank_jar = p2rank_home / "bin" / "p2rank.jar"
        
        if sys.platform == "win32":
            if prank_bat.exists():
                print(f"✓ prank.bat 存在")
            else:
                print(f"✗ prank.bat 不存在")
                all_ok = False
        else:
            if prank_script.exists():
                print(f"✓ prank 脚本存在")
            else:
                print(f"✗ prank 脚本不存在")
                all_ok = False
        
        if p2rank_jar.exists():
            print(f"✓ p2rank.jar 存在")
        else:
            print(f"✗ p2rank.jar 不存在")
            all_ok = False
    else:
        print(f"✗ P2Rank目录不存在: {p2rank_home.absolute()}")
        all_ok = False
    
    print()
    
    # 检查示例文件
    example_file = script_dir / "examples" / "pocket" / "protein.pdb"
    if example_file.exists():
        print(f"✓ 示例文件存在: {example_file.absolute()}")
        file_size = example_file.stat().st_size
        print(f"  文件大小: {file_size} 字节")
    else:
        print(f"✗ 示例文件不存在: {example_file.absolute()}")
        # 示例文件不存在不影响Java环境检查
    
    print()
    print("=" * 60)
    
    if all_ok:
        print("✓ 所有检查通过！Java环境配置正确。")
        return 0
    else:
        print("✗ 发现问题，请参考《口袋预测问题诊断与解决方案.md》进行修复。")
        return 1


if __name__ == "__main__":
    try:
        exit_code = check_java_environment()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n诊断被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n诊断过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

