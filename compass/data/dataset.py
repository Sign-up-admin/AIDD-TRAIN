import logging
import os
import shutil
from multiprocessing import Pool

import torch
from torch_geometric.data import Dataset
from tqdm import tqdm

from .processing import process_item
from ..config import CONFIG  # Directly import the configuration


def _process_and_save_helper(args):
    """Helper function for parallel processing to be used by the Pool."""
    item_info, dest_path, pre_transform = args
    try:
        data = process_item(item_info)
        if data is not None:
            if pre_transform: 
                data = pre_transform(data)
            torch.save(data, dest_path)
            return 1, None  # Return success and no error message
        # If process_item returns None, it means the item was intentionally skipped.
        return 0, f"PDB {item_info.get('pdb_code', 'N/A')} was skipped during processing."
    except Exception as e:
        # Catch any unexpected errors during processing of a single item
        pdb_code = item_info.get('pdb_code', 'N/A')
        logging.error(f"Critical error processing {pdb_code}: {e}", exc_info=True)
        return 0, f"Critical error processing {pdb_code}."

class PDBBindDataset(Dataset):
    """
    PyTorch Geometric Dataset for the PDBbind dataset.

    Handles processing of raw PDB files into graph-based data objects, 
    saving them to disk and loading them for training.
    """
    def __init__(self, root, data_paths, num_workers, transform=None, pre_transform=None):
        """
        Initializes the dataset.

        Args:
            root (str): The root directory where the dataset should be stored.
            data_paths (list): A list of dictionaries, each containing metadata for a single data item.
            num_workers (int): The number of worker processes to use for data processing.
            transform (callable, optional): A function/transform that takes in an 
                `torch_geometric.data.Data` object and returns a transformed version. 
                The data object will be transformed before every access. Defaults to None.
            pre_transform (callable, optional): A function/transform that takes in an 
                `torch_geometric.data.Data` object and returns a transformed version. 
                The data object will be transformed before being saved to disk. Defaults to None.
        """
        self.data_paths = data_paths
        self.num_workers = num_workers
        force_reprocess_flag = CONFIG.get('force_data_reprocessing', False)

        if force_reprocess_flag and os.path.exists(root):
            print(f"\n--- FORCE REPROCESSING ENABLED ---")
            print(f"Deleting existing root data directory: {root}")
            shutil.rmtree(root)
            print("--- Deletion complete. Starting fresh processing. ---\n")

        super(PDBBindDataset, self).__init__(root, transform, pre_transform)

    @property
    def raw_file_names(self):
        """The name of the files in the `raw_dir` to determine whether processing is needed."""
        return [item['pdb_code'] for item in self.data_paths]

    @property
    def processed_file_names(self):
        """The name of the files in the `processed_dir` which must be present to skip processing."""
        return [os.path.join(item['year_dir'], f"{item['pdb_code']}.pt") for item in self.data_paths]

    def process(self):
        """Processes the raw data and saves it to the `processed_dir`."""
        print("Starting data processing...")

        tasks = []
        for item, rel_path in zip(self.data_paths, self.processed_file_names):
            dest_path = os.path.join(self.processed_dir, rel_path)
            if not os.path.exists(dest_path):
                tasks.append((item, dest_path, self.pre_transform))

        for _, dest_path, _ in tasks:
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)

        if not tasks: 
            print("All data is up-to-date.")
            return

        print(f"Processing {len(tasks)} items...")
        success_count = 0
        if self.num_workers > 1:
            print(f"Processing data in parallel with {self.num_workers} cores...")
            with Pool(processes=self.num_workers) as pool:
                pbar = tqdm(pool.imap_unordered(_process_and_save_helper, tasks), total=len(tasks), desc="Processing Raw Data")
                for idx, (result, error_msg) in enumerate(pbar):
                    success_count += result
                    if error_msg: logging.warning(error_msg)
                    # Update progress if logger has progress tracker
                    if hasattr(self, '_logger') and hasattr(self._logger, 'progress_tracker'):
                        self._logger.progress_tracker.update_data_processing(
                            completed=idx + 1,
                            total=len(tasks),
                            message=f"Processing {idx + 1}/{len(tasks)} data points"
                        )
        else:
            print("Processing data sequentially...")
            for idx, task in enumerate(tqdm(tasks, desc="Processing Raw Data")):
                result, error_msg = _process_and_save_helper(task)
                success_count += result
                if error_msg: logging.warning(error_msg)
                # Update progress if logger has progress tracker
                if hasattr(self, '_logger') and hasattr(self._logger, 'progress_tracker'):
                    self._logger.progress_tracker.update_data_processing(
                        completed=idx + 1,
                        total=len(tasks),
                        message=f"Processing {idx + 1}/{len(tasks)} data points"
                    )

        print("--- Processing Finished ---")
        skipped_count = len(tasks) - success_count
        print(f"Successfully processed {success_count}/{len(tasks)} items ({skipped_count} skipped).")

    def len(self):
        """Returns the number of examples in the dataset."""
        return len(self.processed_file_names)

    def get(self, idx):
        """
        Gets the data object at the specified index.

        Args:
            idx (int): The index of the data object to retrieve.

        Returns:
            torch_geometric.data.Data: The data object at the specified index, or None if it fails to load.
        """
        try:
            path = os.path.join(self.processed_dir, self.processed_file_names[idx])
            if not os.path.exists(path): return None
            data = torch.load(path, weights_only=False)
            if data is not None:
                data.pdb_code = self.data_paths[idx].get('pdb_code', f'index_{idx}')
                
                if CONFIG['diffusion']['use_two_stage_diffusion']:
                    # Stage 1: Perturb atom positions slightly to diffuse bond lengths
                    if CONFIG['diffusion']['stage1']['enabled']:
                        noise_level = CONFIG['diffusion']['stage1']['noise_level']
                        noise = torch.randn_like(data.pos) * noise_level
                        data.pos += noise
                    
                    # Stage 2: Perturb atom types and coordinates
                    if CONFIG['diffusion']['stage2']['enabled']:
                        # Perturb coordinates more significantly
                        noise_level = CONFIG['diffusion']['stage2']['noise_level']
                        pos_noise = torch.randn_like(data.pos) * noise_level
                        data.pos += pos_noise
                        
                        # Perturb atom features (types).
                        x_noise = torch.randn_like(data.x) * noise_level
                        data.x += x_noise

            return data
        except (RuntimeError, EOFError, AttributeError, FileNotFoundError) as e:
            pdb_code = self.data_paths[idx].get('pdb_code', f'index {idx}')
            logging.warning(f"Skipping corrupt or incomplete data file for PDB {pdb_code} ({e}): {self.processed_file_names[idx]}")
            return None
