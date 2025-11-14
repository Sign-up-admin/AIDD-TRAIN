"""
Entry point for COMPASS service.
"""

import sys
import os
from pathlib import Path

# Get project root (absolute path)
project_root = Path(__file__).parent.parent.resolve()

# Add project root to Python path if not already present
project_root_str = str(project_root)
if project_root_str not in sys.path:
    sys.path.insert(0, project_root_str)

# Also ensure PYTHONPATH environment variable includes project root
pythonpath = os.environ.get("PYTHONPATH", "")
pythonpath_paths = [p for p in pythonpath.split(os.pathsep) if p] if pythonpath else []
if project_root_str not in pythonpath_paths:
    # Update PYTHONPATH environment variable
    if pythonpath:
        os.environ["PYTHONPATH"] = os.pathsep.join([project_root_str] + pythonpath_paths)
    else:
        os.environ["PYTHONPATH"] = project_root_str

# Verify critical imports can be resolved
try:
    from compass.service.server import main
except ImportError as e:
    print(f"Error: Failed to import compass.service.server: {e}")
    print(f"Project root: {project_root_str}")
    print(f"Python path: {sys.path[:3]}...")  # Show first 3 entries
    print(f"PYTHONPATH: {os.environ.get('PYTHONPATH', 'Not set')}")
    sys.exit(1)

if __name__ == "__main__":
    main()
