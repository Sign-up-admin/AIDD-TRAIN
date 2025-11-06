"""
Entry point for COMPASS service.
"""

import sys
import os
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from compass.service.server import main

if __name__ == "__main__":
    main()
