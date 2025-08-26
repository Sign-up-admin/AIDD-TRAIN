import os
import torch

def _save_checkpoint(state, directory, filename, logger=None):
    """Saves checkpoint to file."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory, filename)
    torch.save(state, filepath)
    if logger:
        logger.log(f"Checkpoint saved to {filepath}")


def _load_checkpoint(directory, device, logger=None):
    """Loads the most recent checkpoint from a directory."""
    if not os.path.isdir(directory):
        if logger: logger.log(f"=> No checkpoint directory found at '{directory}'")
        return None

    checkpoint_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.pth.tar')]
    if not checkpoint_files:
        if logger: logger.log(f"=> No checkpoint files (.pth.tar) found in '{directory}'")
        return None

    latest_checkpoint_path = max(checkpoint_files, key=os.path.getmtime)
    if not latest_checkpoint_path:
        if logger: logger.log(f"=> No valid checkpoints found.")
        return None

    if logger: logger.log(f"=> Attempting to load latest checkpoint: '{os.path.basename(latest_checkpoint_path)}'")
    try:
        checkpoint = torch.load(latest_checkpoint_path, map_location=device)
        if logger: logger.log(f"=> Successfully loaded checkpoint from epoch {checkpoint.get('epoch', 'N/A')}.")
        return checkpoint
    except Exception as e:
        if logger:
            logger.log_error(f"Error loading checkpoint '{os.path.basename(latest_checkpoint_path)}': {e}")
        return None
