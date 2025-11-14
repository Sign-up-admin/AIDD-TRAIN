"""快速检查测试环境是否就绪"""
import requests
import sys

def check_service(url, name):
    try:
        response = requests.get(url, timeout=2)
        if response.status_code == 200:
            print(f"[OK] {name} service is running")
            return True
        else:
            print(f"[X] {name} service returned status {response.status_code}")
            return False
    except Exception as e:
        print(f"[X] {name} service is not running: {e}")
        return False

print("=" * 50)
print("Testing Environment Check")
print("=" * 50)

registry_ok = check_service("http://localhost:8500/health", "Registry")
compass_ok = check_service("http://localhost:8080/health", "COMPASS")

print("\n" + "=" * 50)
if registry_ok and compass_ok:
    print("All services are ready for testing!")
    print("\nYou can now run:")
    print("  1. python test_training_stop_api_only.py")
    print("  2. python debug_training_stop.py")
    sys.exit(0)
else:
    print("Some services are not running.")
    print("\nPlease start services:")
    print("  1. Run start_all_services.bat")
    print("  2. Or start services manually")
    sys.exit(1)










