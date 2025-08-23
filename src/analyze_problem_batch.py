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
    # Fix for std deviation on single-element tensors
    if tensor.numel() > 1:
        logging.info(f"Std value: {torch.std(tensor.float())}")
    else:
        logging.info(f"Std value: 0.0 (single element)")
    logging.info("-" * (20 + len(tensor_name)))

def check_for_duplicate_positions(pos_tensor):
    """Checks for duplicate coordinates in the position tensor."""
    if not isinstance(pos_tensor, torch.Tensor) or pos_tensor.dim() != 2 or pos_tensor.shape[1] != 3:
        return
    
    logging.info("--- Checking for Duplicate Atom Positions ---")
    unique_pos, counts = torch.unique(pos_tensor, dim=0, return_counts=True)
    duplicates = unique_pos[counts > 1]

    if duplicates.shape[0] > 0:
        logging.warning("!!! Duplicate atom positions found! This is a likely cause of NaNs. !!!")
        logging.warning(f"Number of unique duplicate positions: {duplicates.shape[0]}")
        for i in range(duplicates.shape[0]):
            duplicate_coord = duplicates[i]
            num_occurrences = counts[counts > 1][i]
            logging.warning(f"  - Position {duplicate_coord.tolist()} appears {num_occurrences} times.")
    else:
        logging.info("No duplicate atom positions found.")
    logging.info("---------------------------------------------")


def main():
    parser = argparse.ArgumentParser(description="Analyze a saved PyTorch Geometric data batch.")
    parser.add_argument("batch_path", type=str, help="Path to the saved batch file (.pt).")
    args = parser.parse_args()

    logging.info(f"Loading batch from: {args.batch_path}")
    try:
        data = torch.load(args.batch_path, weights_only=False)
        logging.info("Batch loaded successfully.")
    except Exception as e:
        logging.error(f"Failed to load batch file: {e}")
        return

    logging.info(f"Batch type: {type(data)}")
    if hasattr(data, 'pdb_code'):
        logging.info(f"PDB codes in batch: {data.pdb_code}")
    if hasattr(data, 'num_graphs'):
        logging.info(f"Number of graphs in batch: {data.num_graphs}")

    # --- Detailed Analysis ---
    if hasattr(data, 'pos'):
        check_for_duplicate_positions(data.pos)

    for key in data.keys():
        tensor = data[key]
        analyze_tensor(key, tensor)

    logging.info("Analysis complete.")


if __name__ == "__main__":
    main()
