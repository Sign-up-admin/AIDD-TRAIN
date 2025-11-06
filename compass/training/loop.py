from tqdm import tqdm
from torch.amp import autocast
import torch
import os

from .checkpoint import save_checkpoint, create_checkpoint_data
from compass.utils import report_gpu_memory


def _apply_diffusion_noise(data, config, current_stage):
    """Applies diffusion noise to the data based on the current training stage."""
    if not config.get('diffusion', {}).get('use_two_stage_diffusion', False) or current_stage == 0:
        return

    # Stage 1: Perturb atom positions slightly to diffuse bond lengths
    if current_stage == 1 and config['diffusion']['stage1']['enabled']:
        noise_level = config['diffusion']['stage1']['noise_level']
        noise = torch.randn_like(data.pos) * noise_level
        data.pos += noise
    
    # Stage 2: Perturb atom types and coordinates
    elif current_stage == 2 and config['diffusion']['stage2']['enabled']:
        # Perturb coordinates more significantly
        noise_level = config['diffusion']['stage2']['noise_level']
        pos_noise = torch.randn_like(data.pos) * noise_level
        data.pos += pos_noise
        
        # Perturb atom features (types).
        x_noise = torch.randn_like(data.x) * noise_level
        data.x += x_noise

def train_epoch(model, loader, optimizer, device, scaler, config, trainer, logger=None, current_stage=0):
    """
    Runs a single training epoch.
    """
    model.train()
    total_loss, processed_graphs = 0, 0
    optimizer.zero_grad()

    trigger_file = config.get('manual_save_trigger_file')
    checkpoint_dir = config.get('checkpoint_dir')
    grad_accum_steps = config.get('gradient_accumulation_steps')
    gradient_clip_val = config.get('gradient_clip_val', 0)
    debug_mode = config.get('debug_mode', False)

    epoch = trainer.training_state['epoch']
    start_batch = trainer.start_batch

    if debug_mode:
        report_gpu_memory(f"Start of Epoch {epoch}", logger)
    
    pbar = tqdm(enumerate(loader), desc=f"Training Epoch {epoch}", leave=False, total=len(loader), initial=start_batch)
    if start_batch > 0:
        pbar.set_description(f"Resuming Epoch {epoch} (from batch {start_batch})")

    for i, batch in pbar:
        # Check for cancellation
        if hasattr(logger, 'progress_tracker') and logger.progress_tracker.is_cancelled():
            from .exceptions import TrainingCancelled
            raise TrainingCancelled("Training cancelled by user")
        
        if i < start_batch:
            continue
        
        if i == start_batch and start_batch > 0:
            pbar.set_description(f"Training Epoch {epoch}")
        
        trainer.training_state['batch_idx'] = i
        
        # Update progress if logger has progress tracker
        if hasattr(logger, 'progress_tracker'):
            logger.progress_tracker.update_training(
                epoch=epoch,
                total_epochs=config.get('epochs', 1),
                batch=i + 1,
                total_batches=len(loader),
                train_loss=0.0,  # Will be updated after loss calculation
                message=f"Training Epoch {epoch}/{config.get('epochs', 1)}, Batch {i+1}/{len(loader)}"
            )

        data = batch.to(device)
        if data.num_graphs == 0:
            continue

        _apply_diffusion_noise(data, config, current_stage)

        try:
            with autocast(device_type=device.type, dtype=torch.float16):
                output = model(data)
                loss = torch.nn.functional.huber_loss(output, data.y.view(-1, 1)) / grad_accum_steps

            scaler.scale(loss).backward()

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
                        checkpoint_data = create_checkpoint_data(
                            epoch, i + 1, model, optimizer, scaler, 
                            trainer.scheduler, trainer.best_val_loss, interrupted=True
                        )
                        checkpoint_filename = f"checkpoint_epoch_{epoch}_step_{optimizer_step_idx}.pth.tar"
                        save_checkpoint(checkpoint_data, checkpoint_dir, checkpoint_filename, logger)
            
            total_loss += loss.item() * grad_accum_steps * data.num_graphs
            processed_graphs += data.num_graphs
            
            # Update progress with current loss (every N batches to avoid too frequent updates)
            if hasattr(logger, 'progress_tracker') and (i + 1) % 10 == 0:
                current_loss = total_loss / processed_graphs if processed_graphs > 0 else 0.0
                logger.progress_tracker.update_training(
                    epoch=epoch,
                    total_epochs=config.get('epochs', 1),
                    batch=i + 1,
                    total_batches=len(loader),
                    train_loss=current_loss,
                    message=f"Training Epoch {epoch}/{config.get('epochs', 1)}, Batch {i+1}/{len(loader)}, Loss: {current_loss:.4f}"
                )

        except torch.cuda.OutOfMemoryError as e:
            if checkpoint_dir:
                checkpoint_data = create_checkpoint_data(
                    epoch, i, model, optimizer, scaler, 
                    trainer.scheduler, trainer.best_val_loss, interrupted=True
                )
                oom_save_filename = f"oom_checkpoint_epoch_{epoch}_batch_{i}.pth.tar"
                save_checkpoint(checkpoint_data, checkpoint_dir, oom_save_filename, logger)
            raise e

        if trigger_file and checkpoint_dir and os.path.exists(trigger_file):
            checkpoint_data = create_checkpoint_data(
                epoch, i + 1, model, optimizer, scaler, 
                trainer.scheduler, trainer.best_val_loss, interrupted=True
            )
            manual_save_filename = f"manual_checkpoint_epoch_{epoch}_batch_{i}.pth.tar"
            save_checkpoint(checkpoint_data, checkpoint_dir, manual_save_filename, logger)
            
            try:
                os.remove(trigger_file)
            except OSError:
                pass
        
        del loss, output, data
        torch.cuda.empty_cache()

        if debug_mode and (i + 1) % 50 == 0:
            report_gpu_memory(f"After batch {i+1}", logger)

    if debug_mode:
        report_gpu_memory(f"End of Epoch {epoch}", logger)

    return total_loss / processed_graphs if processed_graphs > 0 else 0

def validate_epoch(model, loader, device, logger=None, current_stage=0):
    """
    Runs a single validation epoch.
    """
    model.eval()
    total_loss, processed_graphs = 0, 0
    # For validation, we typically want to evaluate on the original, unmodified data.
    # Therefore, we do not apply diffusion noise here.
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

            # Explicitly free memory to prevent accumulation during the validation loop
            del loss, output, data
            torch.cuda.empty_cache()

    return total_loss / processed_graphs if processed_graphs > 0 else 0
