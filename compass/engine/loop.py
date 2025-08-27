from tqdm import tqdm
from torch.amp import autocast
from torch.profiler import record_function
import torch
import os

from .trainer import _save_checkpoint

def train_epoch(model, loader, optimizer, device, scaler, config, trainer, logger=None):
    """
    Runs a single training epoch.

    Args:
        model (torch.nn.Module): The model to be trained.
        loader (torch.utils.data.DataLoader): The data loader for training data.
        optimizer (torch.optim.Optimizer): The optimizer.
        device (torch.device): The device to run the training on.
        scaler (torch.cuda.amp.GradScaler): The gradient scaler for mixed-precision training.
        config (dict): The configuration dictionary.
        trainer (Trainer): The main trainer object, used to access training state.
        logger (logging.Logger, optional): The logger for logging information. Defaults to None.

    Returns:
        float: The average training loss for the epoch.
    """
    model.train()
    total_loss, processed_graphs = 0, 0
    optimizer.zero_grad()


    trigger_file = config.get('manual_save_trigger_file')
    checkpoint_dir = config.get('checkpoint_dir')
    grad_accum_steps = config.get('gradient_accumulation_steps')
    gradient_clip_val = config.get('gradient_clip_val', 0)

    epoch = trainer.training_state['epoch']
    start_batch = trainer.start_batch
    
    pbar = tqdm(enumerate(loader), desc=f"Training Epoch {epoch}", leave=False, total=len(loader), initial=start_batch)
    if start_batch > 0:
        # This description is shown only for the first few batches of a resumed epoch
        pbar.set_description(f"Resuming Epoch {epoch} (from batch {start_batch})")

    for i, batch in pbar:

        if i < start_batch:
            continue
        
        if i == start_batch and start_batch > 0:
            # After the first resumed batch, revert to the normal description
            pbar.set_description(f"Training Epoch {epoch}")

        trainer.training_state['batch_idx'] = i

        data = batch.to(device)
        if data.num_graphs == 0:
            continue

        try:
            with record_function("model_forward_pass"):
                with autocast(device_type=device.type, dtype=torch.float16):
                    output = model(data)
                    loss = torch.nn.functional.huber_loss(output, data.y.view(-1, 1)) / grad_accum_steps

            with record_function("model_backward_pass"):
                scaler.scale(loss).backward()

            with record_function("optimizer_step"):
                if (i + 1) % grad_accum_steps == 0 or (i + 1) == len(loader):
                    scaler.unscale_(optimizer)
                    if gradient_clip_val > 0:
                        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=gradient_clip_val)
                    scaler.step(optimizer)
                    scaler.update()
                    optimizer.zero_grad()

                    save_every_n_batches = config.get('save_every_n_batches', 0)
                    if save_every_n_batches > 0:
                        optimizer_step_idx = (i + 1) // grad_accum_steps
                        is_last_batch_of_epoch = (i + 1) >= len(loader)

                        if optimizer_step_idx > 0 and optimizer_step_idx % save_every_n_batches == 0 and not is_last_batch_of_epoch:
                            checkpoint_data = {
                                'epoch': epoch, 'batch_idx': i + 1, # Save as next batch to start from
                                'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(),
                                'scaler_state_dict': scaler.state_dict(), 'scheduler_state_dict': trainer.scheduler.state_dict(),
                                'val_loss': trainer.best_val_loss, 'interrupted': True,
                            }
                            checkpoint_filename = f"checkpoint_epoch_{epoch}_step_{optimizer_step_idx}.pth.tar"
                            _save_checkpoint(checkpoint_data, checkpoint_dir, checkpoint_filename, logger)
            
            total_loss += loss.item() * grad_accum_steps * data.num_graphs
            processed_graphs += data.num_graphs

        except torch.cuda.OutOfMemoryError as e:
            if checkpoint_dir:
                checkpoint_data = {
                    'epoch': epoch, 'batch_idx': i, # Save current (failed) batch index
                    'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(),
                    'scaler_state_dict': scaler.state_dict(), 'scheduler_state_dict': trainer.scheduler.state_dict(),
                    'val_loss': trainer.best_val_loss, 'interrupted': True,
                }
                oom_save_filename = f"oom_checkpoint_epoch_{epoch}_batch_{i}.pth.tar"
                _save_checkpoint(checkpoint_data, checkpoint_dir, oom_save_filename, logger)
            raise e

        if trigger_file and checkpoint_dir and os.path.exists(trigger_file):
            checkpoint_data = {
                'epoch': epoch, 'batch_idx': i + 1, # Save as next batch to start from
                'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(),
                'scaler_state_dict': scaler.state_dict(), 'scheduler_state_dict': trainer.scheduler.state_dict(),
                'val_loss': trainer.best_val_loss, 'interrupted': True,
            }
            manual_save_filename = f"manual_checkpoint_epoch_{epoch}_batch_{i}.pth.tar"
            _save_checkpoint(checkpoint_data, checkpoint_dir, manual_save_filename, logger)
            
            try:
                os.remove(trigger_file)
            except OSError:
                pass


    return total_loss / processed_graphs if processed_graphs > 0 else 0

def validate_epoch(model, loader, device, logger=None):
    """
    Runs a single validation epoch.

    Args:
        model (torch.nn.Module): The model to be validated.
        loader (torch.utils.data.DataLoader): The data loader for validation data.
        device (torch.device): The device to run the validation on.
        logger (logging.Logger, optional): The logger for logging information. Defaults to None.

    Returns:
        float: The average validation loss for the epoch.
    """
    model.eval()
    total_loss, processed_graphs = 0, 0
    with torch.no_grad():
        for batch in tqdm(loader, desc="Validation", leave=False):
            if batch is None: continue
            data = batch.to(device)
            if data.num_graphs == 0: continue
            with autocast(device_type=device.type, dtype=torch.float16):
                output = model(data)
                loss = torch.nn.functional.huber_loss(output, data.y.view(-1, 1))
            
            if not torch.isnan(loss):
                total_loss += loss.item() * data.num_graphs
                processed_graphs += data.num_graphs
    return total_loss / processed_graphs if processed_graphs > 0 else 0
