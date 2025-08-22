import argparse
import os
import torch
from rdkit import Chem

from src.model import ViSNetPDB
from src.data_processing import get_protein_graph, get_ligand_graph, get_interaction_graph
from torch_geometric.data import Data
from config import CONFIG

def process_single_pair(protein_path, ligand_path, cutoff=8.0):
    """
    Processes a single protein-ligand pair into a PyG Data object.
    This function replicates the core logic from PDBBindDataset.process.
    """
    if not os.path.exists(protein_path) or not os.path.exists(ligand_path):
        raise FileNotFoundError(f"Could not find protein '{protein_path}' or ligand '{ligand_path}'")

    # 1. Load ligand
    ligand_mol = None
    if ligand_path.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(ligand_path, removeHs=False)
        if supplier:
            ligand_mol = supplier[0]
    elif ligand_path.endswith('.mol2'):
        ligand_mol = Chem.MolFromMol2File(ligand_path, removeHs=False)

    if ligand_mol is None:
        raise ValueError(f"RDKit could not process the ligand file: {ligand_path}")

    # 2. Create graphs
    ligand_graph = get_ligand_graph(ligand_mol)
    protein_graph = get_protein_graph(protein_path)
    interaction_graph = get_interaction_graph(ligand_mol, protein_path, cutoff)

    # 3. Combine into a single Data object
    data = Data()
    data.ligand_pos = ligand_graph['pos']
    data.ligand_v = ligand_graph['v']
    data.ligand_edge_index = ligand_graph['edge_index']
    data.ligand_edge_attr = ligand_graph['edge_attr']

    data.protein_pos = protein_graph['pos']
    data.protein_v = protein_graph['v']

    data.interaction_edge_index = interaction_graph['edge_index']
    data.interaction_edge_attr = interaction_graph['edge_attr']
    
    return data

def main():
    parser = argparse.ArgumentParser(description="Predict binding affinity for a given protein-ligand pair.")
    parser.add_argument('--protein_path', type=str, required=True, help='Path to the protein PDB file.')
    parser.add_argument('--ligand_path', type=str, required=True, help='Path to the ligand SDF or MOL2 file.')
    parser.add_argument('--checkpoint', type=str, default='checkpoints/model_best.pth.tar', help='Path to the model checkpoint file.')
    args = parser.parse_args()

    # --- 1. Setup --- 
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    if not os.path.exists(args.checkpoint):
        print(f"FATAL: Checkpoint file not found at '{args.checkpoint}'")
        print("Please run training first (python main.py) to generate a model checkpoint.")
        return

    # --- 2. Load Model ---
    print("Loading trained model...")
    model_params = {
        'hidden_channels': CONFIG.get('visnet_hidden_channels', 128),
        'num_layers': CONFIG.get('visnet_num_layers', 6),
        'num_rbf': CONFIG.get('visnet_num_rbf', 64),
        'cutoff': CONFIG.get('visnet_cutoff', 8.0),
        'max_num_neighbors': CONFIG.get('visnet_max_neighbors', 32)
    }
    model = ViSNetPDB(**model_params).to(device)

    checkpoint_data = torch.load(args.checkpoint, map_location=device)
    model.load_state_dict(checkpoint_data['model_state_dict'])
    model.eval() # Set model to evaluation mode
    print(f"Model loaded successfully from epoch {checkpoint_data['epoch']} with validation loss {checkpoint_data['val_loss']:.4f}")

    # --- 3. Preprocess Data ---
    print(f"Processing input pair: Protein='{args.protein_path}', Ligand='{args.ligand_path}'")
    try:
        data = process_single_pair(args.protein_path, args.ligand_path, cutoff=model_params['cutoff'])
        data = data.to(device)
    except (FileNotFoundError, ValueError) as e:
        print(f"FATAL: Error during data processing: {e}")
        return

    # --- 4. Predict ---
    print("Making prediction...")
    with torch.no_grad():
        # The model expects a batch, so we create a batch of size 1
        from torch_geometric.loader import DataLoader
        loader = DataLoader([data], batch_size=1, shuffle=False)
        batch = next(iter(loader))
        
        prediction = model(batch)
    
    print("\n--- Prediction Result ---")
    print(f"Predicted -log(Ki/Kd): {prediction.item():.4f}")
    print("-------------------------")

if __name__ == '__main__':
    main()
