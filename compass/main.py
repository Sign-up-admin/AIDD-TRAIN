import os
import signal
from torch_geometric.data import Batch

def collate_filter_none(batch):
    """
    Filters out None values from a batch and returns a new batch.

    Args:
        batch (list): A list of data samples.

    Returns:
        torch_geometric.data.Batch or None: A batch object containing the valid samples,
                                            or None if the batch is empty after filtering.
    """
    batch = list(filter(lambda x: x is not None, batch))
    if not batch:
        return None
    return Batch.from_data_list(batch)


def worker_init_fn(_):
    """
    Prevents worker processes from catching KeyboardInterrupt.
    This is a common solution for multiprocessing data loading issues.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def main():
    """
    Main function to run the training pipeline.

    This function sets up the configuration, data loaders, model, and trainer,
    and then starts the training process.
    """
    import torch
    from torch.utils.data import random_split
    from torch_geometric.loader import DataLoader

    from .config import get_config
    from .data.loader.paths import get_pdb_info, get_data_paths
    from .data.dataset import PDBBindDataset
    from .training.model import ViSNetPDB
    from .training.engine import Trainer
    from .utils import set_seed, get_file_hash
    from .logger import TrainingLogger
    
    config = get_config()

    # --- Dynamic Directory Setup ---
    log_dir = config['log_dir']
    checkpoint_dir = config['checkpoint_dir']
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.makedirs(config['processed_data_dir'], exist_ok=True)

    # --- Logger Setup ---
    logger = TrainingLogger(log_dir=log_dir)
    logger.log("--- Compass Training Process Started ---")

    # --- Hardware & Config Logging ---
    set_seed(config.get('seed', 42), logger)

    logger.log("--- Configuration Loaded ---")
    logger.log(f"Run Name: {config.get('run_name', 'default_run')}")
    logger.log(f"Execution Mode: {config.get('execution_mode', 'N/A')}")
    logger.log(f"Effective Batch Size: {config['batch_size'] * config.get('gradient_accumulation_steps', 1)}")
    logger.log(f"Target Epochs: {config['epochs']}")
    logger.log("--------------------------")

    logger.log("Step 1: Parsing PDBbind index...")
    pdb_info = get_pdb_info(config['index_file'])
    all_data_paths = get_data_paths(pdb_info, config['dataset_path'])
    logger.log(f"-> Found {len(all_data_paths)} total pairs with existing data files.")

    logger.log("Step 2: Verifying data consistency and creating dataset...")

    data_processing_script_path = os.path.join(os.path.dirname(__file__), 'data', 'processing.py')
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

    loader_num_workers = config.get('loader_num_workers', os.cpu_count())
    init_fn = worker_init_fn if loader_num_workers > 0 else None

    pin_memory = True if device.type == 'cuda' else False
    train_loader = DataLoader(train_dataset, batch_size=config['batch_size'], shuffle=True, num_workers=loader_num_workers, pin_memory=pin_memory, collate_fn=collate_filter_none, worker_init_fn=init_fn)
    val_loader = DataLoader(val_dataset, batch_size=config['batch_size'], shuffle=False, num_workers=loader_num_workers, pin_memory=pin_memory, collate_fn=collate_filter_none, worker_init_fn=init_fn)
    logger.log(f"-> Train: {len(train_dataset)} | Validation: {len(val_dataset)}")

    logger.log("Step 4: Setting up model...")
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

    trainer = Trainer(
        config=config,
        model=model,
        train_loader=train_loader,
        val_loader=val_loader,
        device=device,
        logger=logger
    )
    
    trainer.run()