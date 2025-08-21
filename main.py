
import os
import re
import math
import random
import platform
import logging
from io import StringIO
from multiprocessing import Pool, set_start_method
from tqdm import tqdm

# --- Dependency Imports ---
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch.utils.data import random_split
from torch_geometric.data import Dataset, Data, Batch
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool
from torch_geometric.utils import add_self_loops
from rdkit import Chem, rdBase
from Bio.PDB import PDBParser
from scipy.spatial import KDTree

# --- Local Imports ---
from config import CONFIG

# --- Custom Data Class for Protein-Ligand Complexes ---
class ProteinLigandData(Data):
    def __inc__(self, key, value, *args, **kwargs):
        if key == 'protein_edge_index':
            return self.protein_x.size(0)
        if key == 'ligand_edge_index':
            return self.ligand_x.size(0)
        return super().__inc__(key, value, *args, **kwargs)

# --- Logging Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Suppress RDKit logs to keep output clean and focused on critical errors.
rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')

# ==============================================================================
# PART 1: DATA PATH AND BINDING VALUE PARSING
# ==============================================================================

def get_pdb_info(index_file):
    pdb_info = {}
    with open(index_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    pdb_info[parts[0]] = {'year': int(parts[2]), 'binding_data': parts[3]}
                except (ValueError, IndexError):
                    pass # Ignore malformed lines
    return pdb_info

def get_data_paths(pdb_info, dataset_path):
    data_paths = []
    year_to_dir = {y: '1981-2000' for y in range(1981, 2001)}
    year_to_dir.update({y: '2001-2010' for y in range(2001, 2011)})
    year_to_dir.update({y: '2011-2019' for y in range(2011, 2020)})

    for code, info in tqdm(pdb_info.items(), desc="Verifying data paths"):
        year_dir = year_to_dir.get(info['year'], '2011-2019')
        assert isinstance(code, str)
        base_path = os.path.join(dataset_path, year_dir, code)
        protein_path = os.path.join(base_path, code + '_protein.pdb')
        
        ligand_path = None
        sdf_path = os.path.join(base_path, code + '_ligand.sdf')
        mol2_path = os.path.join(base_path, code + '_ligand.mol2')
        if os.path.exists(sdf_path):
            ligand_path = sdf_path
        elif os.path.exists(mol2_path):
            ligand_path = mol2_path

        if os.path.exists(protein_path) and ligand_path:
            data_paths.append({
                'pdb_code': code,
                'protein_path': protein_path,
                'ligand_path': ligand_path,
                'binding_data': info['binding_data'],
                'year_dir': year_dir  # Store the year directory for structured output
            })
    return data_paths

def parse_binding_data(binding_str):
    match = re.search(r"([0-9.]+)([a-zA-Z]+)", binding_str)
    if not match: return None
    
    value, unit = float(match.group(1)), match.group(2)
    unit_map = {'mM': 1e-3, 'uM': 1e-6, 'nM': 1e-9, 'pM': 1e-12, 'fM': 1e-15}
    molar_value = value * unit_map.get(unit, 1.0)
    return -math.log10(molar_value) if molar_value > 0 else 0.0

# ==============================================================================
# PART 2: MOLECULE-TO-GRAPH CONVERSION (ROBUST)
# ==============================================================================

# Define a shared set of elements for one-hot encoding to ensure consistency.
ELEMENTS = ['C', 'O', 'N', 'S', 'P', 'H', 'F', 'Cl', 'Br', 'I', 'UNK']
ELEMENT_MAP = {el: i for i, el in enumerate(ELEMENTS)}

def get_atom_features(atom):
    feat = [0] * len(ELEMENTS)
    element = atom.GetSymbol()
    feat[ELEMENT_MAP.get(element, len(ELEMENTS) - 1)] = 1
    return feat

def ligand_mol_to_graph(mol):
    if mol is None: return None
    try:
        if mol.GetNumConformers() == 0: return None
        atom_features = torch.tensor([get_atom_features(atom) for atom in mol.GetAtoms()], dtype=torch.float)
        rows, cols = [], []
        for bond in mol.GetBonds():
            start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            rows.extend([start, end])
            cols.extend([end, start])
        edge_index = torch.tensor([rows, cols], dtype=torch.long)
        edge_index, _ = add_self_loops(edge_index, num_nodes=atom_features.size(0))
        pos = torch.tensor(mol.GetConformer().GetPositions(), dtype=torch.float)
    except (ValueError, AttributeError, RuntimeError):
        return None
    return Data(x=atom_features, edge_index=edge_index, pos=pos)

def protein_structure_to_graph(structure, cutoff=8.0):
    atoms_to_include = [atom for residue in structure.get_residues() if residue.id[0] == ' ' for atom in residue.get_atoms()]
    if not atoms_to_include:
        return None

    atom_features_list, positions_list = [], []
    for atom in atoms_to_include:
        element = atom.element.strip().upper() if atom.element else ''
        feat = [0] * len(ELEMENTS)
        feat[ELEMENT_MAP.get(element, len(ELEMENTS) - 1)] = 1
        atom_features_list.append(feat)
        positions_list.append(atom.get_coord())

    atom_features = torch.tensor(atom_features_list, dtype=torch.float)
    positions_np = np.array(positions_list)

    # Use KDTree for memory-efficient neighbor search, avoiding O(N^2) memory
    kdtree = KDTree(positions_np)
    edge_list = kdtree.query_ball_tree(kdtree, r=cutoff, p=2.0)

    rows, cols = [], []
    for i, neighbors in enumerate(edge_list):
        for j in neighbors:
            if i < j: # Avoid self-loops and duplicate edges
                rows.append(i)
                cols.append(j)

    # Add edges in both directions
    edge_index = torch.tensor([rows + cols, cols + rows], dtype=torch.long)
    edge_index, _ = add_self_loops(edge_index, num_nodes=atom_features.size(0))
    pos = torch.from_numpy(positions_np).float()

    return Data(x=atom_features, edge_index=edge_index, pos=pos)

def count_pdb_atoms(pdb_path):
    try:
        with open(pdb_path, 'r') as f:
            return sum(1 for line in f if line.startswith('ATOM'))
    except IOError:
        return 0

def process_item(item):
    pdb_code = item.get('pdb_code', 'N/A')
    try:
        protein_path = item['protein_path']
        if os.path.getsize(protein_path) > 100 * 1024 * 1024: # 100MB limit
            logging.warning(f"Skipping {pdb_code}: Protein file is very large (>100MB).")
            return None

        # Stricter atom count limit for safety with KDTree
        if count_pdb_atoms(protein_path) > 15000:
            logging.warning(f"Skipping {pdb_code}: Protein has too many atoms (>15,000) for safe processing.")
            return None

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_code, protein_path)

        ligand = None
        ligand_path = item['ligand_path']
        if ligand_path.endswith('.mol2'):
            ligand = Chem.MolFromMol2File(ligand_path, removeHs=True, sanitize=True)
        elif ligand_path.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(ligand_path, removeHs=True, sanitize=True)
            if suppl: ligand = next(suppl, None)
        
        if ligand is None or ligand.GetNumAtoms() == 0:
            logging.warning(f"Skipping {pdb_code}: Failed to load ligand or ligand has no atoms.")
            return None

        protein_graph = protein_structure_to_graph(structure)
        if protein_graph is None: return None

        ligand_graph = ligand_mol_to_graph(ligand)
        if ligand_graph is None: return None

        y = parse_binding_data(item['binding_data'])
        if y is None: return None
        
        return ProteinLigandData(
            protein_x=protein_graph.x, protein_pos=protein_graph.pos, protein_edge_index=protein_graph.edge_index,
            ligand_x=ligand_graph.x, ligand_pos=ligand_graph.pos, ligand_edge_index=ligand_graph.edge_index,
            y=torch.tensor([y], dtype=torch.float),
            pdb_code=pdb_code
        )
    except Exception as e:
        logging.error(f"ERROR for PDB {pdb_code}: {e}")
        return None

# ==============================================================================
# PART 3: PYTORCH GEOMETRIC DATASET
# ==============================================================================

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
        
        print(f"--- Processing Finished ---")
        skipped_count = len(tasks) - success_count
        print(f"Successfully processed {success_count}/{len(tasks)} items ({skipped_count} skipped due to errors).")

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        try:
            path = os.path.join(self.processed_dir, self.processed_file_names[idx])
            if not os.path.exists(path): return None
            return torch.load(path, weights_only=False)
        except (RuntimeError, EOFError, AttributeError):
            pdb_code = self.data_paths[idx].get('pdb_code', f'index {idx}')
            logging.warning(f"Skipping corrupt or incomplete data file for PDB {pdb_code}: {self.processed_file_names[idx]}")
            return None

# ==============================================================================
# PART 4: GRAPH NEURAL NETWORK MODEL
# ==============================================================================

class Net(nn.Module):
    def __init__(self, protein_in_dim, ligand_in_dim, hidden_dim=64, dropout_rate=0.5):
        super(Net, self).__init__()
        self.protein_conv1 = GCNConv(protein_in_dim, hidden_dim)
        self.protein_conv2 = GCNConv(hidden_dim, hidden_dim * 2)
        self.protein_conv3 = GCNConv(hidden_dim * 2, hidden_dim * 4)
        self.ligand_conv1 = GCNConv(ligand_in_dim, hidden_dim)
        self.ligand_conv2 = GCNConv(hidden_dim, hidden_dim * 2)
        self.ligand_conv3 = GCNConv(hidden_dim * 2, hidden_dim * 4)
        self.fc1 = nn.Linear(hidden_dim * 4 * 2, 1024)
        self.fc2 = nn.Linear(1024, 512)
        self.out = nn.Linear(512, 1)
        self.dropout_rate = dropout_rate

    def forward(self, data):
        p_x, p_edge_index, p_batch = data.protein_x, data.protein_edge_index, data.protein_x_batch
        l_x, l_edge_index, l_batch = data.ligand_x, data.ligand_edge_index, data.ligand_x_batch
        
        p_x = F.relu(self.protein_conv1(p_x, p_edge_index))
        p_x = F.relu(self.protein_conv2(p_x, p_edge_index))
        p_x = F.relu(self.protein_conv3(p_x, p_edge_index))
        p_x = global_mean_pool(p_x, p_batch)
        
        l_x = F.relu(self.ligand_conv1(l_x, l_edge_index))
        l_x = F.relu(self.ligand_conv2(l_x, l_edge_index))
        l_x = F.relu(self.ligand_conv3(l_x, l_edge_index))
        l_x = global_mean_pool(l_x, l_batch)
        
        x = torch.cat([p_x, l_x], dim=1)
        x = F.relu(self.fc1(x))
        x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = F.relu(self.fc2(x))
        return self.out(x)

# ==============================================================================
# PART 5: TRAINING AND EVALUATION LOOPS
# ==============================================================================

def train(model, loader, optimizer, device):
    model.train()
    total_loss, processed_graphs = 0, 0
    for batch in tqdm(loader, desc="Training", leave=False):
        if batch is None: continue
        data = batch.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = F.mse_loss(output, data.y.view(-1, 1))
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * data.num_graphs
        processed_graphs += data.num_graphs
    return total_loss / processed_graphs if processed_graphs > 0 else 0

def test(model, loader, device):
    model.eval()
    total_loss, processed_graphs = 0, 0
    with torch.no_grad():
        for batch in tqdm(loader, desc="Validation", leave=False):
            if batch is None: continue
            data = batch.to(device)
            output = model(data)
            loss = F.mse_loss(output, data.y.view(-1, 1))
            total_loss += loss.item() * data.num_graphs
            processed_graphs += data.num_graphs
    return total_loss / processed_graphs if processed_graphs > 0 else 0

# ==============================================================================
# PART 6: MAIN EXECUTION SCRIPT
# ==============================================================================

def main():
    config = CONFIG
    print("--- Configuration Loaded ---")
    print(f"Processing Cores: {config['processing_num_workers']}")
    print(f"Loader Cores: {config['loader_num_workers']}")
    print(f"Batch Size: {config['batch_size']}")
    print(f"Epochs: {config['epochs']}")
    print("--------------------------")
    os.makedirs(config['processed_data_dir'], exist_ok=True)

    print("Step 1: Parsing PDBbind index...")
    pdb_info = get_pdb_info(config['index_file'])
    all_data_paths = get_data_paths(pdb_info, config['dataset_path'])
    print(f"-> Found {len(all_data_paths)} total pairs with existing data files.")

    print("Step 2: Verifying data consistency and creating dataset...")

    DATA_PROCESSING_VERSION = f"v7_kdtree_{len(ELEMENTS)}"
    version_file_path = os.path.join(config['processed_data_dir'], 'processing_version.txt')

    processed_files_exist = any(os.path.exists(os.path.join(config['processed_data_dir'], item['year_dir'], f"{item['pdb_code']}.pt")) for item in all_data_paths)

    if processed_files_exist:
        if os.path.exists(version_file_path):
            with open(version_file_path, 'r') as f:
                stored_version = f.read().strip()
            if stored_version != DATA_PROCESSING_VERSION:
                print("\n" + "="*70)
                print("FATAL: Data processing logic has changed.")
                print(f"  - Stored data version: {stored_version}")
                print(f"  - Current code version: {DATA_PROCESSING_VERSION}")
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
        with open(version_file_path, 'w') as f: f.write(DATA_PROCESSING_VERSION)

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
    train_dataset, val_dataset = random_split(dataset, [train_size, val_size])
    
    loader_num_workers = config['loader_num_workers']
    if platform.system() == 'Windows' and loader_num_workers > 0:
        print("Warning: Using multiple workers for DataLoader on Windows can be unstable.")

    pin_memory = True if device.type == 'cuda' else False
    train_loader = DataLoader(train_dataset, batch_size=config['batch_size'], shuffle=True, num_workers=loader_num_workers, pin_memory=pin_memory, follow_batch=['protein_x', 'ligand_x'])
    val_loader = DataLoader(val_dataset, batch_size=config['batch_size'], shuffle=False, num_workers=loader_num_workers, pin_memory=pin_memory, follow_batch=['protein_x', 'ligand_x'])
    print(f"-> Train: {len(train_dataset)} | Validation: {len(val_dataset)}")

    print("Step 4: Setting up model and optimizer...")
    if len(train_dataset) == 0:
        print("FATAL: Training dataset is empty. Cannot proceed.")
        return
        
    try:
        sample_data = next(iter(train_loader))
        protein_in_dim = sample_data.protein_x.shape[1]
        ligand_in_dim = sample_data.ligand_x.shape[1]
    except (StopIteration, IndexError, AttributeError):
        print("FATAL: Cannot determine model dimensions from the dataset. Check for empty validation set.")
        return

    model = Net(
        protein_in_dim=protein_in_dim, 
        ligand_in_dim=ligand_in_dim,
        dropout_rate=config['dropout_rate']
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'])
    print(f"-> Model initialized with {sum(p.numel() for p in model.parameters()):,} parameters.")

    print("Step 5: Starting training...")
    for epoch in range(1, config['epochs'] + 1):
        train_loss = train(model, train_loader, optimizer, device)
        val_loss = test(model, val_loader, device)
        print(f"Epoch {epoch:02d}/{config['epochs']} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

    print("--- Training Finished ---")

# ==============================================================================
# SCRIPT ENTRY POINT
# ==============================================================================

if __name__ == '__main__':
    try:
        set_start_method('spawn')
    except RuntimeError:
        pass
    
    main()
