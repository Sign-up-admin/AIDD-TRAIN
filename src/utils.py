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

def save_checkpoint(state, directory, filename="checkpoint.pth.tar", logger=None):
    """Saves checkpoint to file."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory, filename)
    torch.save(state, filepath)
    if logger:
        logger.log(f"Checkpoint saved to {filepath}")

def load_checkpoint(directory, device, logger=None):
    """Loads the most recent checkpoint from a directory, ensuring robust recovery."""
    if not os.path.isdir(directory):
        if logger: logger.log(f"=> No checkpoint directory found at '{directory}'")
        return None

    checkpoint_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.pth.tar')]
    
    if not checkpoint_files:
        if logger: logger.log(f"=> No checkpoint files (.pth.tar) found in '{directory}'")
        return None

    latest_checkpoint_path = max(checkpoint_files, key=os.path.getmtime)

    if latest_checkpoint_path is None:
        if logger: logger.log(f"=> No valid checkpoints found after checking modification times.")
        return None

    if logger: logger.log(f"=> Attempting to load latest checkpoint: '{os.path.basename(latest_checkpoint_path)}'")
    try:
        checkpoint = torch.load(latest_checkpoint_path, map_location=device)
        if logger: logger.log(f"=> Successfully loaded checkpoint from epoch {checkpoint.get('epoch', 'N/A')}.")
        return checkpoint
    except Exception as e:
        if logger:
            logger.log_error(f"Error loading checkpoint '{os.path.basename(latest_checkpoint_path)}': {e}")
            logger.log_error("The latest checkpoint may be corrupt. Please check the file or remove it to load a previous one.")
        return None

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
