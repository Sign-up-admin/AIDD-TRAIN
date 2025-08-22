import os
import logging
import random
import torch
import numpy as np

# ==============================================================================
# PART 5B: CHECKPOINTING HELPERS
# ==============================================================================

def save_checkpoint(state, directory, filename="checkpoint.pth.tar"):
    """Saves checkpoint to file."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory, filename)
    torch.save(state, filepath)
    logging.info(f"Checkpoint saved to {filepath}")

def load_checkpoint(directory, device):
    """Loads the most recent checkpoint from a directory, ensuring robust recovery."""
    if not os.path.isdir(directory):
        logging.info(f"=> No checkpoint directory found at '{directory}'")
        return None

    checkpoint_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.pth.tar')]
    
    if not checkpoint_files:
        logging.info(f"=> No checkpoint files (.pth.tar) found in '{directory}'")
        return None

    latest_checkpoint_path = max(checkpoint_files, key=os.path.getmtime)

    if latest_checkpoint_path is None:
        logging.info(f"=> No valid checkpoints found after checking modification times.")
        return None

    logging.info(f"=> Attempting to load latest checkpoint: '{os.path.basename(latest_checkpoint_path)}'")
    try:
        checkpoint = torch.load(latest_checkpoint_path, map_location=device)
        logging.info(f"=> Successfully loaded checkpoint from epoch {checkpoint.get('epoch', 'N/A')}.")
        return checkpoint
    except Exception as e:
        logging.error(f"Error loading checkpoint '{os.path.basename(latest_checkpoint_path)}': {e}")
        logging.error("The latest checkpoint may be corrupt. Please check the file or remove it to load a previous one.")
        return None

# ==============================================================================
# PART 6: REPRODUCIBILITY
# ==============================================================================

def set_seed(seed):
    """Sets the random seed for reproducibility across all relevant libraries."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    logging.info(f"Global random seed set to {seed} for reproducibility.")
