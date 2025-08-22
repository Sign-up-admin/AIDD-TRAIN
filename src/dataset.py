import os
import logging
from multiprocessing import Pool
import torch
from torch_geometric.data import Dataset
from tqdm import tqdm

from src.data_processing import process_item

def _process_and_save_helper(args):
    item_info, dest_path, pre_transform = args
    data = process_item(item_info)
    if data is not None:
        if pre_transform: data = pre_transform(data)
        torch.save(data, dest_path)
        return (1, None)
    return (0, f"PDB {item_info.get('pdb_code', 'N/A')} was skipped.")

class PDBBindDataset(Dataset):
    def __init__(self, root, data_paths, num_workers, transform=None, pre_transform=None):
        self.data_paths = data_paths
        self.num_workers = num_workers
        super(PDBBindDataset, self).__init__(root, transform, pre_transform)

    @property
    def raw_file_names(self):
        return [item['pdb_code'] for item in self.data_paths]

    @property
    def processed_file_names(self):
        return [os.path.join(item['year_dir'], f"{item['pdb_code']}.pt") for item in self.data_paths]

    def process(self):
        print("Starting data processing...")
        tasks = [(item, os.path.join(self.processed_dir, rel_path), self.pre_transform)
                 for item, rel_path in zip(self.data_paths, self.processed_file_names)
                 if not os.path.exists(os.path.join(self.processed_dir, rel_path))]

        for item, dest_path, _ in tasks:
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)

        if not tasks: 
            print("All data is already processed.")
            return

        print(f"Processing {len(tasks)} new items...")
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
        print(f"Successfully processed {success_count}/{len(tasks)} items ({skipped_count} skipped due to errors).")

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        try:
            path = os.path.join(self.processed_dir, self.processed_file_names[idx])
            if not os.path.exists(path): return None
            data = torch.load(path, weights_only=False)
            # Attach the PDB code to the data object for easier debugging
            if data is not None:
                data.pdb_code = self.data_paths[idx].get('pdb_code', f'index_{idx}')
            return data
        except (RuntimeError, EOFError, AttributeError):
            pdb_code = self.data_paths[idx].get('pdb_code', f'index {idx}')
            logging.warning(f"Skipping corrupt or incomplete data file for PDB {pdb_code}: {self.processed_file_names[idx]}")
            return None
