import os
import sys
import signal
import platform
import hashlib
from multiprocessing import set_start_method

import torch
import torch_geometric
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
from src.utils import save_checkpoint, load_checkpoint, set_seed


def collate_filter_none(batch):
    """Filters out None values from a batch and returns a new batch."""
    batch = list(filter(lambda x: x is not None, batch))
    # If the batch is empty after filtering, return None
    if not batch:
        return None
    return torch_geometric.data.Batch.from_data_list(batch)


def get_file_hash(filepath):
    """Computes the SHA256 hash of a file to be used as a version identifier."""
    if not os.path.exists(filepath):
        return None
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def main():
    config = CONFIG

    # Dynamically set processing and loader workers based on CPU cores
    num_cpu_cores = os.cpu_count()
    if num_cpu_cores is not None:
        if num_cpu_cores < 4:
            config['processing_num_workers'] = 2
            config['loader_num_workers'] = 2
        else:
            config['processing_num_workers'] = num_cpu_cores
            config['loader_num_workers'] = num_cpu_cores
    else:
        print("Warning: Could not detect number of CPU cores. Using default worker settings from config.py.")

    # --- Reproducibility ---
    set_seed(config.get('seed', 42))

    print("--- Configuration Loaded ---")
    print(f"Processing Cores: {config['processing_num_workers']}")
    print(f"Loader Cores: {config['loader_num_workers']}")
    print(f"Batch Size: {config['batch_size']} (Effective: {config['batch_size'] * config['gradient_accumulation_steps']})")
    print(f"Epochs: {config['epochs']}")
    print("--------------------------")
    
    config['checkpoint_dir'] = os.path.join(os.path.dirname(config.get('processed_data_dir', '.')), 'checkpoints')
    os.makedirs(config['processed_data_dir'], exist_ok=True)
    os.makedirs(config['checkpoint_dir'], exist_ok=True)

    print("Step 1: Parsing PDBbind index...")
    pdb_info = get_pdb_info(config['index_file'])
    all_data_paths = get_data_paths(pdb_info, config['dataset_path'])
    print(f"-> Found {len(all_data_paths)} total pairs with existing data files.")

    print("Step 2: Verifying data consistency and creating dataset...")

    # Use the hash of the data processing script as the version
    data_processing_script_path = os.path.join(os.path.dirname(__file__), 'src', 'data_processing.py')
    data_processing_version = get_file_hash(data_processing_script_path)
    if data_processing_version is None:
        print(f"FATAL: Could not find data processing script at '{data_processing_script_path}' to generate version hash.")
        return
        
    version_file_path = os.path.join(config['processed_data_dir'], 'processing_version.txt')

    # Check for existing data and compare versions
    if any(os.path.exists(os.path.join(config['processed_data_dir'], item['year_dir'], f"{item['pdb_code']}.pt")) for item in all_data_paths):
        if os.path.exists(version_file_path):
            with open(version_file_path, 'r') as f:
                stored_version = f.read().strip()
            if stored_version != data_processing_version:
                print("\n" + "="*70)
                print("FATAL: Data processing logic has changed (code hash mismatch).")
                print(f"  - Stored data version: {stored_version[:12]}...")
                print(f"  - Current code version: {data_processing_version[:12]}...")
                print("  - Please DELETE the processed data directory to allow reprocessing:")
                print(f"    {config['processed_data_dir']}")
                print("="*70 + "\n")
                return
        else:
            print("\n" + "="*70)
            print("FATAL: Found old, unversioned processed data.")
            print("  - The data processing logic has been updated for consistency.")
            print("  - Please DELETE the processed data directory to allow reprocessing:")
            print(f"    {config['processed_data_dir']}")
            print("="*70 + "\n")
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
    print(f"-> Found {len(dataset)} processable data points.")

    if len(dataset) == 0:
        print("FATAL: No valid data points found after filtering. Cannot proceed.")
        return

    print("Step 3: Splitting data and creating loaders...")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"-> Using device: {device}")

    train_size = int(config['train_split'] * len(dataset))
    val_size = len(dataset) - train_size
    
    generator = torch.Generator().manual_seed(config.get('seed', 42))
    train_subset, val_subset = random_split(dataset, [train_size, val_size], generator=generator)
    train_dataset = dataset.index_select(train_subset.indices)
    val_dataset = dataset.index_select(val_subset.indices)
    
    # Force num_workers to 0 on Windows to prevent instability
    loader_num_workers = config['loader_num_workers']
    if platform.system() == 'Windows' and loader_num_workers > 0:
        print("Warning: Using multiple workers for DataLoader on Windows can be unstable. Forcing num_workers to 0.")
        loader_num_workers = 0

    pin_memory = True if device.type == 'cuda' else False
    train_loader = DataLoader(train_dataset, batch_size=config['batch_size'], shuffle=True, num_workers=loader_num_workers, pin_memory=pin_memory, collate_fn=collate_filter_none)
    val_loader = DataLoader(val_dataset, batch_size=config['batch_size'], shuffle=False, num_workers=loader_num_workers, pin_memory=pin_memory, collate_fn=collate_filter_none)
    print(f"-> Train: {len(train_dataset)} | Validation: {len(val_dataset)}")

    print("Step 4: Setting up model, optimizer, and scheduler...")
    if len(train_dataset) == 0:
        print("FATAL: Training dataset is empty. Cannot proceed.")
        return
        
    model = ViSNetPDB(
        hidden_channels=config.get('visnet_hidden_channels', 128),
        num_layers=config.get('visnet_num_layers', 6),
        num_rbf=config.get('visnet_num_rbf', 64),
        cutoff=config.get('visnet_cutoff', 8.0),
        max_num_neighbors=config.get('visnet_max_neighbors', 32)
    ).to(device)

    # Increased epsilon for Adam optimizer for better stability with float16/AMP
    optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], eps=1e-7)
    scaler = GradScaler()
    scheduler = ReduceLROnPlateau(optimizer, 'min', patience=3, factor=0.5)
    print(f"-> Model initialized with {sum(p.numel() for p in model.parameters()):,} parameters.")

    start_epoch = 1
    start_batch = 0
    best_val_loss = float('inf')
    
    training_state = {'epoch': 0, 'batch_idx': 0}

    def graceful_exit_handler(sig=None, _frame=None):
        """Handles graceful shutdown on SIGTERM or KeyboardInterrupt."""
        if sig == signal.SIGTERM:
            print("\n\n--- SIGTERM received. Saving final state... ---")
        elif sig == signal.SIGINT:
            print("\n\n--- SIGINT (Ctrl+C / Stop button) received. Saving final state... ---")
        else:
            print("\n\n--- Training interrupted by user (KeyboardInterrupt). Saving final state... ---")

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
            save_checkpoint(interrupt_checkpoint_data, config['checkpoint_dir'], 'INTERRUPTED.pth.tar')
            print(f"--- Final state for epoch {epoch_to_save}, batch {batch_idx_to_save} saved to INTERRUPTED.pth.tar. Exiting. ---")
        else:
            print("--- Interruption occurred before training loop. No state to save. Exiting. ---")
        
        sys.exit(0)

    signal.signal(signal.SIGTERM, graceful_exit_handler)
    signal.signal(signal.SIGINT, graceful_exit_handler)

    checkpoint = load_checkpoint(config['checkpoint_dir'], device)
    if checkpoint:
        try:
            model.load_state_dict(checkpoint['model_state_dict'])
            optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            
            if checkpoint.get('interrupted', False):
                start_epoch = checkpoint['epoch']
                start_batch = checkpoint.get('batch_idx', 0)
                print(f"--> NOTE: Resuming from an interrupted training run. Restarting epoch {start_epoch} from batch {start_batch}.")
            else:
                start_epoch = checkpoint['epoch'] + 1
                start_batch = 0

            best_val_loss = checkpoint.get('val_loss', float('inf'))
            if 'scaler_state_dict' in checkpoint:
                scaler.load_state_dict(checkpoint['scaler_state_dict'])
                print("--> Resumed GradScaler state.")
            if 'scheduler_state_dict' in checkpoint:
                scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
                print("--> Resumed LR Scheduler state.")
            
            print(f"--> Checkpoint states (model, optimizer, etc.) loaded successfully.")
            print(f"--> Will start/resume training at epoch {start_epoch}. Last best validation loss: {best_val_loss:.4f}")

        except (KeyError, RuntimeError) as e:
            print(f"Could not load checkpoint due to an error: {e}. Starting from scratch.")
            start_epoch = 1
            start_batch = 0
            best_val_loss = float('inf')
    else:
        print("--> No checkpoint found, starting from scratch.")

    print("Step 5: Starting training...")
    try:
        for epoch in range(start_epoch, config['epochs'] + 1):
            training_state['epoch'] = epoch
            training_state['batch_idx'] = 0

            is_profiling_epoch = (epoch == 1 and start_batch == 0 and config.get('profile', True))
            
            if is_profiling_epoch:
                print("--- Profiling enabled for the first epoch ---")
                with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA], record_shapes=True) as prof:
                    train_loss = train(model, train_loader, optimizer, device, scaler, config['gradient_accumulation_steps'], epoch, best_val_loss, config, scheduler, training_state, start_batch=start_batch)
                
                print("--- Profiler Results ---")
                print(prof.key_averages().table(sort_by="cuda_time_total", row_limit=15))
            else:
                train_loss = train(model, train_loader, optimizer, device, scaler, config['gradient_accumulation_steps'], epoch, best_val_loss, config, scheduler, training_state, start_batch=start_batch)

            start_batch = 0

            val_loss = test(model, val_loader, device)
            
            scheduler.step(val_loss)
            
            print(f"Epoch {epoch:02d}/{config['epochs']} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

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
            
            save_checkpoint(checkpoint_data, config['checkpoint_dir'], 'checkpoint.pth.tar')
            if is_best:
                save_checkpoint(checkpoint_data, config['checkpoint_dir'], 'model_best.pth.tar')
                print(f"-> New best model saved with validation loss: {val_loss:.4f}")

    except KeyboardInterrupt:
        graceful_exit_handler()

    print("--- Training Finished ---")

if __name__ == '__main__':
    try:
        set_start_method('spawn')
    except RuntimeError:
        pass
    
    main()
