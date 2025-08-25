import os
import random
import torch
import numpy as np
import hashlib

def get_file_hash(filepath):
    """Computes the SHA256 hash of a file to be used as a version identifier."""
    if not os.path.exists(filepath):
        return None
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def set_seed(seed, logger=None):
    """Sets the random seed for reproducibility across all relevant libraries."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    if logger:
        logger.log(f"Global random seed set to {seed} for reproducibility.")
