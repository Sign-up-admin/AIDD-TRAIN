import os
import logging
import hashlib
import shutil
from multiprocessing import Pool

import torch
from torch_geometric.data import Dataset
from tqdm import tqdm

from src.data_processing import process_item
from config import CONFIG # Directly import the configuration

def _process_and_save_helper(args):
    """Helper function for parallel processing to be used by the Pool."""
    item_info, dest_path, pre_transform = args
    try:
        data = process_item(item_info)
        if data is not None:
            if pre_transform: 
                data = pre_transform(data)
            torch.save(data, dest_path)
            return (1, None) # Return success and no error message
        # If process_item returns None, it means the item was intentionally skipped.
        return (0, f"PDB {item_info.get('pdb_code', 'N/A')} was skipped during processing.")
    except Exception as e:
        # Catch any unexpected errors during processing of a single item
        pdb_code = item_info.get('pdb_code', 'N/A')
        logging.error(f"Critical error processing {pdb_code}: {e}", exc_info=True)
        return (0, f"Critical error processing {pdb_code}.")

class PDBBindDataset(Dataset):
    def __init__(self, root, data_paths, num_workers, transform=None, pre_transform=None):
        self.data_paths = data_paths
        self.num_workers = num_workers
        force_reprocess_flag = CONFIG.get('force_data_reprocessing', False)

        # --- FINAL, CORRECT FIX ---
        # The check must happen on the `root` directory itself, *before* the parent
        # class __init__ is called. The parent class __init__ decides whether to call
        # self.process(). If we don't delete the entire root directory here, the parent
        # class will see the old files and incorrectly skip the processing step.
        if force_reprocess_flag and os.path.exists(root):
            print(f"\n--- FORCE REPROCESSING ENABLED ---")
            print(f"Deleting existing root data directory: {root}")
            shutil.rmtree(root)
            print("--- Deletion complete. Starting fresh processing. ---\n")

        super(PDBBindDataset, self).__init__(root, transform, pre_transform)

    @property
    def raw_file_names(self):
        return [item['pdb_code'] for item in self.data_paths]

    @property
    def processed_file_names(self):
        return [os.path.join(item['year_dir'], f"{item['pdb_code']}.pt") for item in self.data_paths]

    def process(self):
        print("Starting data processing...")

        tasks = []
        for item, rel_path in zip(self.data_paths, self.processed_file_names):
            dest_path = os.path.join(self.processed_dir, rel_path)
            # The check is simple now, as we know the folder is clean if it was forced.
            if not os.path.exists(dest_path):
                tasks.append((item, dest_path, self.pre_transform))

        # Create parent directories for the files that will be created.
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
                for result, error_msg in pbar:
                    success_count += result
                    if error_msg: logging.warning(error_msg)
        else:
            print("Processing data sequentially...")
            for task in tqdm(tasks, desc="Processing Raw Data"):
                result, error_msg = _process_and_save_helper(task)
                success_count += result
                if error_msg: logging.warning(error_msg)

        print("--- Processing Finished ---")
        skipped_count = len(tasks) - success_count
        print(f"Successfully processed {success_count}/{len(tasks)} items ({skipped_count} skipped).")

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        try:
            path = os.path.join(self.processed_dir, self.processed_file_names[idx])
            if not os.path.exists(path): return None
            data = torch.load(path, weights_only=False)
            if data is not None:
                data.pdb_code = self.data_paths[idx].get('pdb_code', f'index_{idx}')
            return data
        except (RuntimeError, EOFError, AttributeError, FileNotFoundError) as e:
            pdb_code = self.data_paths[idx].get('pdb_code', f'index {idx}')
            logging.warning(f"Skipping corrupt or incomplete data file for PDB {pdb_code} ({e}): {self.processed_file_names[idx]}")
            return None
