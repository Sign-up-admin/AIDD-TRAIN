import os
import logging
import torch
import torch.nn.functional as F
from tqdm import tqdm
from torch.amp import autocast
from torch.profiler import record_function

from src.utils import save_checkpoint

def train(model, loader, optimizer, device, scaler, grad_accum_steps, epoch, best_val_loss, config, scheduler, training_state, start_batch=0):
    model.train()
    total_loss, processed_graphs = 0, 0
    optimizer.zero_grad()  # Reset gradients at the start of the epoch

    # --- Manual Save Trigger Setup ---
    trigger_file = config.get('manual_save_trigger_file')
    checkpoint_dir = config.get('checkpoint_dir')

    pbar = tqdm(enumerate(loader), desc=f"Training Epoch {epoch}", leave=False, total=len(loader), initial=start_batch)
    if start_batch > 0:
        pbar.set_description(f"Resuming Epoch {epoch} (from batch {start_batch})")

    for i, batch in pbar:
        if i < start_batch:
            continue

        # Skip batch if it is None (e.g., due to filtering in collate_fn)
        if batch is None:
            continue
        
        if i == start_batch and start_batch > 0:
            pbar.set_description(f"Training Epoch {epoch}")

        training_state['batch_idx'] = i

        data = batch.to(device)
        # Skip batch if it's empty to prevent NaN loss
        if data.num_graphs == 0:
            continue

        # record_function is used by the profiler to label code regions.
        with record_function("model_forward_pass"):
            # --- Automatic Mixed Precision (AMP) ---
            with autocast(device_type=device.type, dtype=torch.float16):
                output = model(data)
                loss = F.mse_loss(output, data.y.view(-1, 1)) / grad_accum_steps

            # --- NaN/Inf Debugging ---
            if torch.isnan(output).any() or torch.isinf(output).any() or torch.isnan(loss).any() or torch.isinf(loss).any():
                logging.warning(f"NaN/Inf detected at epoch {epoch}, batch {i}.")
                if hasattr(data, 'pdb_code'):
                    logging.warning(f"Problematic PDBs in batch: {data.pdb_code}")

                if config.get('debug_mode', False):
                    debug_dir = os.path.join(config.get('checkpoint_dir', '.'), 'debug_batches')
                    os.makedirs(debug_dir, exist_ok=True)
                    problem_batch_path = os.path.join(debug_dir, f"problem_batch_epoch_{epoch}_batch_{i}.pt")
                    torch.save(data.to('cpu'), problem_batch_path)
                    logging.error(f"NaN detected. Saved problematic batch to {problem_batch_path}. Stopping training.")
                    raise RuntimeError(f"NaN detected in batch {i} of epoch {epoch}. Batch saved for debugging.")
                else:
                    logging.warning("Skipping batch.")
                    continue

        with record_function("model_backward_pass"):
            scaler.scale(loss).backward()

        with record_function("optimizer_step"):
            if (i + 1) % grad_accum_steps == 0 or (i + 1) == len(loader):
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                scaler.step(optimizer)
                scaler.update()
                optimizer.zero_grad()
        
        # --- Manual Save Trigger Check ---
        if trigger_file and checkpoint_dir and os.path.exists(trigger_file):
            logging.info(f"--> Manual save triggered by file: {trigger_file}")
            
            checkpoint_data = {
                'epoch': epoch,
                'batch_idx': i,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scaler_state_dict': scaler.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'val_loss': best_val_loss,
                'interrupted': True,
            }
            
            manual_save_filename = f"manual_checkpoint_epoch_{epoch}_batch_{i}.pth.tar"
            save_checkpoint(checkpoint_data, checkpoint_dir, manual_save_filename)
            
            try:
                os.remove(trigger_file)
                logging.info(f"--> Removed trigger file: {trigger_file}")
            except OSError as e:
                logging.error(f"--> Error removing trigger file {trigger_file}: {e}")

        total_loss += loss.item() * grad_accum_steps * data.num_graphs
        processed_graphs += data.num_graphs

    return total_loss / processed_graphs if processed_graphs > 0 else 0

def test(model, loader, device):
    model.eval()
    total_loss, processed_graphs = 0, 0
    with torch.no_grad():
        for batch in tqdm(loader, desc="Validation", leave=False):
            # Skip batch if it is None (e.g., due to filtering in collate_fn)
            if batch is None:
                continue
            data = batch.to(device)
            # Skip batch if it's empty to prevent NaN loss
            if data.num_graphs == 0:
                continue
            with autocast(device_type=device.type, dtype=torch.float16):
                output = model(data)
                # --- NaN/Inf Check ---
                # Check for invalid values in output before calculating loss
                if torch.isnan(output).any() or torch.isinf(output).any():
                    logging.warning(f"NaN/Inf detected in validation output for a batch. Skipping loss calculation for this batch.")
                    if hasattr(data, 'pdb_code'):
                        logging.warning(f"Problematic PDBs in validation batch: {data.pdb_code}")
                    continue
                
                loss = F.mse_loss(output, data.y.view(-1, 1))
            
            # Ensure loss is valid before accumulating
            if not torch.isnan(loss):
                total_loss += loss.item() * data.num_graphs
                processed_graphs += data.num_graphs
    return total_loss / processed_graphs if processed_graphs > 0 else 0
