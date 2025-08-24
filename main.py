import os
import sys
import signal
import platform
from multiprocessing import set_start_method

import torch
from torch_geometric.data import Batch
from torch.utils.data import random_split
from torch_geometric.loader import DataLoader
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.amp import GradScaler
from torch.profiler import profile, ProfilerActivity

from config import CONFIG
from src.data_processing import get_pdb_info, get_data_paths
from src.dataset import PDBBindDataset
from src.model import ViSNetPDB
from src.training import train, test
from src.utils import save_checkpoint, load_checkpoint, set_seed, get_file_hash
from src.hardware_utils import get_hardware_recommendations
from src.logger import TrainingLogger


def collate_filter_none(batch):
    """Filters out None values from a batch and returns a new batch."""
    batch = list(filter(lambda x: x is not None, batch))
    if not batch:
        return None
    return Batch.from_data_list(batch)


def worker_init_fn(worker_id):
    """
    Prevents worker processes from catching KeyboardInterrupt.
    This is a common solution for multiprocessing data loading issues.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def main():
    config = CONFIG

    # --- Dynamic Directory Setup ---
    # Create run-specific directories for logs and checkpoints to avoid overwrites.
    run_name = config.get('run_name', 'default_run')
    log_dir = os.path.join('logs', run_name)
    checkpoint_dir = os.path.join('checkpoints', run_name)
    config['checkpoint_dir'] = checkpoint_dir
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.makedirs(config['processed_data_dir'], exist_ok=True)

    # --- Logger Setup ---
    logger = TrainingLogger(log_dir=log_dir)
    logger.log("--- Training Process Started ---")

    # --- Hardware & Config Logging ---
    get_hardware_recommendations(config, logger)
    set_seed(config.get('seed', 42), logger)

    logger.log("--- Configuration Loaded ---")
    logger.log(f"Run Name: {run_name}")
    logger.log(f"Development Mode: {config.get('development_mode', 'N/A')}")
    logger.log(f"Effective Batch Size: {config['batch_size'] * config['gradient_accumulation_steps']}")
    logger.log(f"Target Epochs: {config['epochs']}")
    logger.log("--------------------------")

    logger.log("Step 1: Parsing PDBbind index...")
    pdb_info = get_pdb_info(config['index_file'])
    all_data_paths = get_data_paths(pdb_info, config['dataset_path'])
    logger.log(f"-> Found {len(all_data_paths)} total pairs with existing data files.")

    logger.log("Step 2: Verifying data consistency and creating dataset...")

    data_processing_script_path = os.path.join(os.path.dirname(__file__), 'src', 'data_processing.py')
    data_processing_version = get_file_hash(data_processing_script_path)
    if data_processing_version is None:
        logger.log_error(f"Could not find data processing script at '{data_processing_script_path}' to generate version hash.")
        return

    version_file_path = os.path.join(config['processed_data_dir'], 'processing_version.txt')

    if any(os.path.exists(os.path.join(config['processed_data_dir'], item['year_dir'], f"{item['pdb_code']}.pt")) for item in all_data_paths):
        if os.path.exists(version_file_path):
            with open(version_file_path, 'r') as f:
                stored_version = f.read().strip()
            if stored_version != data_processing_version:
                logger.log_error("="*70)
                logger.log_error("Data processing logic has changed (code hash mismatch).")
                logger.log_error(f"  - Stored data version: {stored_version[:12]}...")
                logger.log_error(f"  - Current code version: {data_processing_version[:12]}...")
                logger.log_error("  - Please DELETE the processed data directory to allow reprocessing:")
                logger.log_error(f"    {config['processed_data_dir']}")
                logger.log_error("="*70)
                return
        else:
            logger.log_error("="*70)
            logger.log_error("Found old, unversioned processed data.")
            logger.log_error("  - The data processing logic has been updated for consistency.")
            logger.log_error("  - Please DELETE the processed data directory to allow reprocessing:")
            logger.log_error(f"    {config['processed_data_dir']}")
            logger.log_error("="*70)
            return

    dataset = PDBBindDataset(
        root=config['processed_data_dir'],
        data_paths=all_data_paths,
        num_workers=config['processing_num_workers']
    )

    if not os.path.exists(version_file_path):
        with open(version_file_path, 'w') as f:
            f.write(data_processing_version)

    valid_indices = [i for i, f in enumerate(dataset.processed_paths) if os.path.exists(f)]
    dataset = dataset.index_select(valid_indices)
    logger.log(f"-> Found {len(dataset)} processable data points.")

    if len(dataset) == 0:
        logger.log_error("No valid data points found after filtering. Cannot proceed.")
        return

    logger.log("Step 3: Splitting data and creating loaders...")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.log(f"-> Using device: {device}")

    train_size = int(config['train_split'] * len(dataset))
    val_size = len(dataset) - train_size

    generator = torch.Generator().manual_seed(config.get('seed', 42))
    train_subset, val_subset = random_split(dataset, [train_size, val_size], generator=generator)
    train_dataset = dataset.index_select(train_subset.indices)
    val_dataset = dataset.index_select(val_subset.indices)

    loader_num_workers = config['loader_num_workers']
    init_fn = worker_init_fn if loader_num_workers > 0 else None

    pin_memory = True if device.type == 'cuda' else False
    train_loader = DataLoader(train_dataset, batch_size=config['batch_size'], shuffle=True, num_workers=loader_num_workers, pin_memory=pin_memory, collate_fn=collate_filter_none, worker_init_fn=init_fn)
    val_loader = DataLoader(val_dataset, batch_size=config['batch_size'], shuffle=False, num_workers=loader_num_workers, pin_memory=pin_memory, collate_fn=collate_filter_none, worker_init_fn=init_fn)
    logger.log(f"-> Train: {len(train_dataset)} | Validation: {len(val_dataset)}")

    logger.log("Step 4: Setting up model, optimizer, and scheduler...")
    if len(train_dataset) == 0:
        logger.log_error("Training dataset is empty. Cannot proceed.")
        return

    model = ViSNetPDB(
        hidden_channels=config.get('visnet_hidden_channels', 128),
        num_layers=config.get('visnet_num_layers', 6),
        num_rbf=config.get('visnet_num_rbf', 64),
        cutoff=config.get('visnet_cutoff', 8.0),
        max_num_neighbors=config.get('max_num_neighbors', 32),
        lmax=config.get('visnet_lmax', 1),
        vecnorm_type=config.get('visnet_vecnorm_type', 'max_min')
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], eps=1e-7, weight_decay=1e-5)
    scaler = GradScaler()
    scheduler = ReduceLROnPlateau(optimizer, 'min', patience=3, factor=0.5)
    logger.log(f"-> Model initialized with {sum(p.numel() for p in model.parameters()):,} parameters.")

    start_epoch = 1
    start_batch = 0
    best_val_loss = float('inf')

    training_state = {'epoch': 0, 'batch_idx': 0}

    def graceful_exit_handler(sig=None, _frame=None):
        """Handles graceful shutdown on SIGTERM or KeyboardInterrupt."""
        if sig == signal.SIGTERM:
            logger.log_warning("\n\n--- SIGTERM received. Saving final state... ---")
        elif sig == signal.SIGINT:
            logger.log_warning("\n\n--- SIGINT (Ctrl+C / Stop button) received. Saving final state... ---")
        else:
            logger.log_warning("\n\n--- Training interrupted by user (KeyboardInterrupt). Saving final state... ---")

        epoch_to_save = training_state.get('epoch', 0)
        batch_idx_to_save = training_state.get('batch_idx', 0)

        if epoch_to_save > 0:
            interrupt_checkpoint_data = {
                'epoch': epoch_to_save,
                'batch_idx': batch_idx_to_save,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scaler_state_dict': scaler.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'val_loss': best_val_loss,
                'interrupted': True
            }
            save_checkpoint(interrupt_checkpoint_data, config['checkpoint_dir'], 'INTERRUPTED.pth.tar', logger)
            logger.log(f"--- Final state for epoch {epoch_to_save}, batch {batch_idx_to_save} saved to INTERRUPTED.pth.tar. Exiting. ---")
        else:
            logger.log_warning("--- Interruption occurred before training loop. No state to save. Exiting. ---")

        sys.exit(0)

    signal.signal(signal.SIGTERM, graceful_exit_handler)
    signal.signal(signal.SIGINT, graceful_exit_handler)

    checkpoint = load_checkpoint(config['checkpoint_dir'], device, logger)
    if checkpoint:
        try:
            model.load_state_dict(checkpoint['model_state_dict'])
            optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

            if checkpoint.get('interrupted', False):
                start_epoch = checkpoint['epoch']
                start_batch = checkpoint.get('batch_idx', 0)
                logger.log_warning(f"--> Resuming from an interrupted training run. Restarting epoch {start_epoch} from batch {start_batch}.")
            else:
                start_epoch = checkpoint['epoch'] + 1
                start_batch = 0

            best_val_loss = checkpoint.get('val_loss', float('inf'))
            if 'scaler_state_dict' in checkpoint:
                scaler.load_state_dict(checkpoint['scaler_state_dict'])
                logger.log("--> Resumed GradScaler state.")
            if 'scheduler_state_dict' in checkpoint:
                scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
                logger.log("--> Resumed LR Scheduler state.")

            logger.log(f"--> Checkpoint states (model, optimizer, etc.) loaded successfully.")
            logger.log(f"--> Will start/resume training at epoch {start_epoch}. Last best validation loss: {best_val_loss:.4f}")

        except (KeyError, RuntimeError) as e:
            logger.log_error(f"Could not load checkpoint due to an error: {e}. Starting from scratch.")
            start_epoch = 1
            start_batch = 0
            best_val_loss = float('inf')
    else:
        logger.log("--> No checkpoint found, starting from scratch.")

    logger.log("Step 5: Starting training...")
    for epoch in range(start_epoch, config['epochs'] + 1):
        training_state['epoch'] = epoch
        training_state['batch_idx'] = 0

        is_profiling_epoch = (epoch == 1 and start_batch == 0 and config.get('profile', False))

        if is_profiling_epoch:
            logger.log("--- Profiling enabled for the first epoch ---")
            with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA], record_shapes=True) as prof:
                train_loss = train(model, train_loader, optimizer, device, scaler, config['gradient_accumulation_steps'], epoch, best_val_loss, config, scheduler, training_state, start_batch=start_batch, logger=logger)

            logger.log("--- Profiler Results ---")
            logger.log(prof.key_averages().table(sort_by="cuda_time_total", row_limit=15))
        else:
            train_loss = train(model, train_loader, optimizer, device, scaler, config['gradient_accumulation_steps'], epoch, best_val_loss, config, scheduler, training_state, start_batch=start_batch, logger=logger)

        start_batch = 0

        val_loss = test(model, val_loader, device, logger=logger)

        scheduler.step(val_loss)

        logger.log(f"Epoch {epoch:02d}/{config['epochs']} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

        is_best = val_loss < best_val_loss
        if is_best:
            best_val_loss = val_loss

        checkpoint_data = {
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'scaler_state_dict': scaler.state_dict(),
            'scheduler_state_dict': scheduler.state_dict(),
            'val_loss': val_loss,
        }

        save_checkpoint(checkpoint_data, config['checkpoint_dir'], 'checkpoint.pth.tar', logger)
        if is_best:
            save_checkpoint(checkpoint_data, config['checkpoint_dir'], 'model_best.pth.tar', logger)
            logger.log(f"-> New best model saved with validation loss: {val_loss:.4f}")

    logger.log("--- Training Finished ---")


if __name__ == '__main__':
    try:
        set_start_method('spawn')
    except RuntimeError:
        pass

    main()
