import os
import torch


def create_checkpoint_data(
    epoch, batch_idx, model, optimizer, scaler, scheduler, val_loss, interrupted=False
):
    """
    Creates a dictionary containing the checkpoint state.

    Args:
        epoch (int): The current epoch number.
        batch_idx (int): The current batch index.
        model (torch.nn.Module): The model.
        optimizer (torch.optim.Optimizer): The optimizer.
        scaler (torch.cuda.amp.GradScaler): The gradient scaler.
        scheduler (torch.optim.lr_scheduler._LRScheduler): The learning rate scheduler.
        val_loss (float): The validation loss.
        interrupted (bool): Whether the training was interrupted.

    Returns:
        dict: A dictionary containing the checkpoint state.
    """
    return {
        "epoch": epoch,
        "batch_idx": batch_idx,
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "scaler_state_dict": scaler.state_dict(),
        "scheduler_state_dict": scheduler.state_dict(),
        "val_loss": val_loss,
        "interrupted": interrupted,
    }


def save_checkpoint(state, directory, filename, logger=None):
    """
    Saves checkpoint to file.

    Args:
        state (dict): The state to save.
        directory (str): The directory to save the checkpoint in.
        filename (str): The name of the checkpoint file.
        logger (Logger, optional): The logger for logging messages. Defaults to None.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory, filename)
    torch.save(state, filepath)
    if logger:
        logger.log(f"Checkpoint saved to {filepath}")


def load_latest_checkpoint_file(directory, device, logger=None):
    """
    Loads the most recent checkpoint from a directory.

    Args:
        directory (str): The directory to load the checkpoint from.
        device (torch.device): The device to map the loaded checkpoint to.
        logger (Logger, optional): The logger for logging messages. Defaults to None.

    Returns:
        dict or None: The loaded checkpoint dictionary, or None if no checkpoint is found.
    """
    if not os.path.isdir(directory):
        if logger:
            logger.log(f"=> No checkpoint directory found at '{directory}'")
        return None

    checkpoint_files = [
        os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".pth.tar")
    ]
    if not checkpoint_files:
        if logger:
            logger.log(f"=> No checkpoint files (.pth.tar) found in '{directory}'")
        return None

    try:
        latest_checkpoint_path = max(checkpoint_files, key=os.path.getmtime)
    except FileNotFoundError:
        if logger:
            logger.log(f"=> Checkpoint file not found, it may have been deleted.")
        return None

    if not latest_checkpoint_path:
        if logger:
            logger.log(f"=> No valid checkpoints found.")
        return None

    if logger:
        logger.log(
            f"=> Attempting to load latest checkpoint: '{os.path.basename(latest_checkpoint_path)}'"
        )
    try:
        checkpoint = torch.load(latest_checkpoint_path, map_location=device)
        if logger:
            logger.log(
                f"=> Successfully loaded checkpoint from epoch {checkpoint.get('epoch', 'N/A')}."
            )
        return checkpoint
    except Exception as e:
        if logger:
            logger.log_error(
                f"Error loading checkpoint '{os.path.basename(latest_checkpoint_path)}': {e}"
            )
        return None


def resume_from_checkpoint(trainer):
    """
    Loads model and training state from the latest checkpoint.

    Args:
        trainer (Trainer): The trainer instance.
    """
    logger = trainer.logger
    checkpoint = load_latest_checkpoint_file(
        trainer.config["checkpoint_dir"], trainer.device, logger
    )
    if not checkpoint:
        logger.log("--> No checkpoint found, starting from scratch.")
        return

    try:
        trainer.model.load_state_dict(checkpoint["model_state_dict"])
        trainer.optimizer.load_state_dict(checkpoint["optimizer_state_dict"])

        if checkpoint.get("interrupted", False):
            trainer.start_epoch = checkpoint["epoch"]
            trainer.start_batch = checkpoint.get("batch_idx", 0)
            logger.log_warning(
                f"--> Resuming from an interrupted run. Restarting epoch {trainer.start_epoch} from batch {trainer.start_batch}."
            )
        else:
            trainer.start_epoch = checkpoint["epoch"] + 1
            trainer.start_batch = 0

        trainer.best_val_loss = checkpoint.get("val_loss", float("inf"))
        if "scaler_state_dict" in checkpoint:
            trainer.scaler.load_state_dict(checkpoint["scaler_state_dict"])
            logger.log("--> Resumed GradScaler state.")
        if "scheduler_state_dict" in checkpoint:
            trainer.scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
            logger.log("--> Resumed LR Scheduler state.")

        logger.log(
            f"--> Checkpoint loaded. Starting at epoch {trainer.start_epoch}. Best val loss: {trainer.best_val_loss:.4f}"
        )

    except (KeyError, RuntimeError) as e:
        logger.log_error(f"Could not load checkpoint due to an error: {e}. Starting from scratch.")
        trainer.start_epoch = 1
        trainer.start_batch = 0
        trainer.best_val_loss = float("inf")
