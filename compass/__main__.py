"""
This file makes the 'compass' package executable.

You can run the main training script from the project root directory using:
python -m compass
"""

import platform
from multiprocessing import set_start_method

from .main import main # noqa: E402

if __name__ == '__main__':
    # Set the multiprocessing start method for non-Linux systems to 'spawn'.
    # This is crucial for preventing deadlocks with CUDA and multiprocessing
    # and must be called within the `if __name__ == '__main__':` block.
    try:
        if platform.system() != 'Linux':
            set_start_method('spawn')
    except RuntimeError:
        # This can happen if the method has already been set.
        pass

    # This entry point is executed when the package is run as a script.
    main()