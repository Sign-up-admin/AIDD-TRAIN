"""Check if services are running."""
import requests
import sys

services = {
    "Registry": "http://localhost:8500/health",
    "COMPASS": "http://localhost:8080/health",
    "Frontend": "http://localhost:8501",
}

all_running = True
for name, url in services.items():
    try:
        r = requests.get(url, timeout=2)
        if r.status_code == 200:
            print(f"[OK] {name} is running (status: {r.status_code})")
        else:
            print(f"[FAIL] {name} returned status: {r.status_code}")
            all_running = False
    except Exception as e:
        print(f"[FAIL] {name} is not running: {str(e)[:50]}")
        all_running = False

sys.exit(0 if all_running else 1)

