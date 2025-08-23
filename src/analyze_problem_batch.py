import torch
import argparse
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def analyze_tensor(tensor_name, tensor):
    """Analyzes a single tensor for NaNs, Infs, and other stats."""
    if not isinstance(tensor, torch.Tensor):
        logging.info(f"Skipping {tensor_name}: Not a tensor.")
        return

    logging.info(f"--- Analyzing Tensor: {tensor_name} ---")
    logging.info(f"Shape: {tensor.shape}")
    logging.info(f"Dtype: {tensor.dtype}")

    has_nan = torch.isnan(tensor).any()
    has_inf = torch.isinf(tensor).any()

    logging.info(f"Contains NaN: {has_nan}")
    logging.info(f"Contains Inf: {has_inf}")

    if has_nan or has_inf:
        logging.warning(f"!!! Problem detected in {tensor_name} !!!")

    if tensor.numel() == 0:
        logging.info("Tensor is empty.")
        return
        
    if tensor.dtype in [torch.long, torch.int, torch.short, torch.bool]:
        logging.info("Skipping stats for integer/boolean tensor.")
        return

    logging.info(f"Min value: {torch.min(tensor)}")
    logging.info(f"Max value: {torch.max(tensor)}")
    logging.info(f"Mean value: {torch.mean(tensor.float())}")
    if tensor.numel() > 1:
        logging.info(f"Std value: {torch.std(tensor.float())}")
    else:
        logging.info(f"Std value: 0.0 (single element)")
    logging.info("-" * (20 + len(tensor_name)))

def check_for_duplicate_positions(pos_tensor):
    """Checks for duplicate coordinates by rounding to a set precision."""
    if not isinstance(pos_tensor, torch.Tensor) or pos_tensor.dim() != 2 or pos_tensor.shape[1] != 3:
        return
    
    logging.info("--- Checking for Duplicate Atom Positions ---")
    seen_positions = {}
    duplicates_found = False
    for i, pos in enumerate(pos_tensor):
        pos_tuple = (round(pos[0].item(), 4), round(pos[1].item(), 4), round(pos[2].item(), 4))
        if pos_tuple in seen_positions:
            duplicates_found = True
            seen_positions[pos_tuple].append(i)
        else:
            seen_positions[pos_tuple] = [i]

    if duplicates_found:
        logging.warning("!!! Duplicate atom positions found! This is a likely cause of NaNs. !!!")
        for pos_tuple, indices in seen_positions.items():
            if len(indices) > 1:
                logging.warning(f"  - Position {list(pos_tuple)} appears {len(indices)} times at indices: {indices}")
    else:
        logging.info("No duplicate atom positions found.")
    logging.info("---------------------------------------------")

def check_for_unknown_atoms(x_tensor, pos_tensor):
    """Checks for atoms classified as UNK and prints their info."""
    if not isinstance(x_tensor, torch.Tensor) or x_tensor.dim() != 2:
        return
    
    logging.info("--- Checking for Unknown (UNK) Atoms ---")
    
    # The last column in the one-hot encoding corresponds to UNK.
    # Find indices where the last column is 1.
    unknown_atom_indices = (x_tensor[:, -1] == 1).nonzero(as_tuple=True)[0]
    
    if unknown_atom_indices.numel() > 0:
        logging.warning(f"!!! {unknown_atom_indices.numel()} Unknown (UNK) atoms found! This is a likely cause of NaNs. !!!")
        for idx in unknown_atom_indices:
            atom_pos = pos_tensor[idx].tolist()
            logging.warning(f"  - Atom at index {idx.item()} with position {atom_pos} is UNK.")
    else:
        logging.info("No unknown (UNK) atoms found.")
    logging.info("----------------------------------------")

def main():
    parser = argparse.ArgumentParser(description="Analyze a saved PyTorch Geometric data batch.")
    parser.add_argument("batch_path", type=str, help="Path to the saved batch file (.pt).")
    args = parser.parse_args()

    logging.info(f"Loading batch from: {args.batch_path}")
    try:
        # Load the data object. It's not a batch, but a single Data object.
        data = torch.load(args.batch_path, weights_only=False)
        logging.info("Data object loaded successfully.")
    except Exception as e:
        logging.error(f"Failed to load data file: {e}")
        return

    logging.info(f"Data type: {type(data)}")
    if hasattr(data, 'pdb_code'):
        logging.info(f"PDB code: {data.pdb_code}")

    # --- Detailed Analysis ---
    if hasattr(data, 'pos'):
        check_for_duplicate_positions(data.pos)

    if hasattr(data, 'x') and hasattr(data, 'pos'):
        check_for_unknown_atoms(data.x, data.pos)

    for key in data.keys():
        tensor = data[key]
        analyze_tensor(key, tensor)

    logging.info("Analysis complete.")

if __name__ == "__main__":
    main()
