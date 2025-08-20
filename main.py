"""# ==============================================================================
# INSTALLATION GUIDE
# ==============================================================================
#
# Before running, please ensure the following libraries are installed.
# We recommend using a Conda environment.
#
# 1. PyTorch:      conda install pytorch torchvision torchaudio -c pytorch
# 2. RDKit:        conda install -c conda-forge rdkit
# 3. PyG:          Visit https://pyg.org/for/torch/ and get the correct command.
# 4. TQDM:         conda install tqdm
# 5. Biopython:    conda install -c conda-forge biopython
#
# ==============================================================================
"""

import os
import re
import math
import platform
import tempfile
from multiprocessing import Pool, set_start_method, cpu_count
from tqdm import tqdm

# --- Dependency Imports ---
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Dataset, Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool
from rdkit import Chem, rdBase
from Bio.PDB import PDBParser, PDBIO

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
                'binding_data': info['binding_data']
            })
    return data_paths

def parse_binding_data(binding_str):
    match = re.match(r'([A-Za-z]+)=([0-9.]+)([a-zA-Z]+)', binding_str)
    if not match: return None
    
    value, unit = float(match.group(2)), match.group(3)
    unit_map = {'mM': 1e-3, 'uM': 1e-6, 'nM': 1e-9, 'pM': 1e-12, 'fM': 1e-15}
    molar_value = value * unit_map.get(unit, 1.0)
    return -math.log10(molar_value) if molar_value > 0 else 0.0

# ==============================================================================
# PART 2: MOLECULE-TO-GRAPH CONVERSION (ROBUST)
# ==============================================================================

def get_atom_features(atom):
    return [
        atom.GetAtomicNum(), atom.GetDegree(), atom.GetFormalCharge(),
        atom.GetHybridization(), atom.GetIsAromatic(), atom.GetTotalNumHs()
    ]

def mol_to_graph(mol):
    if mol is None: return None
    try:
        atom_features = torch.tensor([get_atom_features(atom) for atom in mol.GetAtoms()], dtype=torch.float)
        
        # Memory-efficient edge index creation to avoid large adjacency matrices.
        rows, cols = [], []
        for bond in mol.GetBonds():
            start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            rows.append(start)
            cols.append(end)
            rows.append(end) # Add reverse bond for undirected graph
            cols.append(start)
        
        edge_index = torch.tensor([rows, cols], dtype=torch.long)
        
        pos = torch.tensor(mol.GetConformer().GetPositions(), dtype=torch.float)
    except (ValueError, AttributeError, RuntimeError): # Added RuntimeError for safety
        return None
    return Data(x=atom_features, edge_index=edge_index, pos=pos)

def process_item(item):
    """Worker function to process a single protein-ligand pair with robust error handling."""
    pdb_code = item.get('pdb_code', 'N/A')
    try:
        # 1. Load Protein using Biopython for robustness
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_code, item['protein_path'])
        
        # RDKit still needs a PDB file format, so we write the parsed structure to a temporary in-memory file
        # This helps clean up many formatting issues that cause RDKit to fail.
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode='w') as tmp:
            io = PDBIO()
            io.set_structure(structure)
            io.save(tmp.name)
            protein = Chem.MolFromPDBFile(tmp.name, removeHs=False, sanitize=True)
        
        os.remove(tmp.name) # Clean up the temporary file

        if protein is None:
            # This is a common failure point. Log it and exit early.
            print(f"Warning: RDKit failed to load protein for {pdb_code} at {item['protein_path']} even after Biopython parsing. Skipping.")
            return None

        # 2. Load Ligand
        ligand = None
        ligand_path = item['ligand_path']
        if ligand_path.endswith('.mol2'):
            ligand = Chem.MolFromMol2File(ligand_path, removeHs=False, sanitize=True)
        elif ligand_path.endswith('.sdf'):
            # Use a supplier for SDF as it can contain multiple molecules.
            suppl = Chem.SDMolSupplier(ligand_path, removeHs=False, sanitize=True)
            if suppl: # Check if the supplier was created successfully
                ligand = next(suppl, None)

        if ligand is None:
            print(f"Warning: RDKit failed to load ligand for {pdb_code} at {ligand_path}. Skipping.")
            return None

        # 3. Convert to Graphs
        protein_graph = mol_to_graph(protein)
        if protein_graph is None:
            print(f"Warning: Failed to convert protein to graph for {pdb_code}. Skipping.")
            return None
        
        ligand_graph = mol_to_graph(ligand)
        if ligand_graph is None:
            print(f"Warning: Failed to convert ligand to graph for {pdb_code}. Skipping.")
            return None

        # 4. Parse Binding Data
        y = parse_binding_data(item['binding_data'])
        if y is None:
            print(f"Warning: Failed to parse binding data for {pdb_code} ('{item['binding_data']}'). Skipping.")
            return None

        # 5. Assemble Final Data Object
        return Data(
            protein_x=protein_graph.x, protein_pos=protein_graph.pos, protein_edge_index=protein_graph.edge_index,
            ligand_x=ligand_graph.x, ligand_pos=ligand_graph.pos, ligand_edge_index=ligand_graph.edge_index,
            y=torch.tensor([y], dtype=torch.float),
            pdb_code=pdb_code # Store PDB code for easier debugging later
        )
    except Exception as e:
        # This is a fallback for any other unexpected errors during processing.
        print(f"ERROR: An unexpected error occurred in process_item for {pdb_code}: {e}")
        return None

# ==============================================================================
# PART 3: PYTORCH GEOMETRIC DATASET
# ==============================================================================

def _process_and_save_helper(args):
    """
    Worker function for multiprocessing. Must be a top-level function to be pickleable.
    Catches all exceptions to prevent worker processes from crashing silently.

    Args:
        args (tuple): A tuple containing index, item_info, processed_dir, and pre_transform function.

    Returns:
        int: 1 if processing is successful, 0 otherwise.
    """
    idx, item_info, processed_dir, pre_transform = args
    pdb_code = item_info.get('pdb_code', 'N/A')
    try:
        data = process_item(item_info)
        if data is not None:
            if pre_transform:
                data = pre_transform(data)
            torch.save(data, os.path.join(processed_dir, f'data_{idx}.pt'))
            return 1
    except Exception as e:
        # Catching all exceptions here is crucial to prevent a single bad file
        # from crashing a worker process and causing the pool to hang.
        print(f"ERROR: Unhandled exception for PDB {pdb_code} (task index {idx}): {e}")
    return 0

class PDBBindDataset(Dataset):
    def __init__(self, root, data_paths, num_workers, transform=None, pre_transform=None):
        self.data_paths = data_paths
        self.num_workers = num_workers
        super(PDBBindDataset, self).__init__(root, transform, pre_transform)

    @property
    def raw_file_names(self):
        return []

    @property
    def processed_file_names(self):
        return [f'data_{i}.pt' for i in range(len(self.data_paths))]

    def process(self):
        print("Starting data processing...")
        
        tasks = []
        for i, item in enumerate(self.data_paths):
            if not os.path.exists(os.path.join(self.processed_dir, f'data_{i}.pt')):
                tasks.append((i, item, self.processed_dir, self.pre_transform))
        
        if not tasks:
            print("All data is already processed.")
            return

        print(f"Processing {len(tasks)} new items...")

        success_count = 0
        if self.num_workers > 1:
            print(f"Processing data in parallel with {self.num_workers} cores...")
            # Using imap_unordered is more memory-efficient for large task lists,
            # as it doesn't need to hold all results in memory before returning.
            with Pool(processes=self.num_workers) as pool:
                pbar = tqdm(pool.imap_unordered(_process_and_save_helper, tasks), total=len(tasks), desc="Processing Raw Data")
                for result in pbar:
                    success_count += result
        else:
            print("Processing data sequentially...")
            for task in tqdm(tasks, desc="Processing Raw Data"):
                success_count += _process_and_save_helper(task)
        
        print(f"Successfully processed {success_count}/{len(tasks)} new items.")
        if success_count == 0 and len(tasks) > 0:
            print("WARNING: No new data could be processed. Check input file formats or paths.")

    def len(self):
        # By convention, the length of the dataset is the number of processed files.
        # We assume that if processing was successful, all files exist.
        return len(self.processed_file_names)

    def get(self, idx):
        try:
            data = torch.load(os.path.join(self.processed_dir, f'data_{idx}.pt'))
            if self.transform:
                data = self.transform(data)
            return data
        except (FileNotFoundError, RuntimeError, EOFError) as e:
            print(f"Warning: Could not load data_{idx}.pt. Error: {e}. Skipping this item.")
            # Returning None is handled by the DataLoader which will filter it out.
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
    total_loss = 0
    for data in tqdm(loader, desc="Training", leave=False):
        if data is None:
            continue
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = F.mse_loss(output, data.y.view(-1, 1))
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * data.num_graphs
    return total_loss / len(loader.dataset)

def test(model, loader, device):
    model.eval()
    total_loss = 0
    with torch.no_grad():
        for data in tqdm(loader, desc="Validation", leave=False):
            if data is None:
                continue
            data = data.to(device)
            output = model(data)
            loss = F.mse_loss(output, data.y.view(-1, 1))
            total_loss += loss.item() * data.num_graphs
    return total_loss / len(loader.dataset)

# ==============================================================================
# PART 6: MAIN EXECUTION SCRIPT
# ==============================================================================

def get_optimal_workers():
    """
    Calculates a conservative number of workers for data preprocessing.
    This is based on the number of CPU cores but is capped to prevent
    excessive memory usage, which is critical on systems with around 16GB of RAM.
    Using too many workers can lead to swapping to virtual memory (SSD),
    slowing down the process significantly.
    """
    cores = cpu_count()
    if cores >= 32:
        return 8  # Cap at 8 for high-core-count machines
    elif cores >= 16:
        return 4
    elif cores >= 8:
        return 2
    else:
        return 1 # Default to 1 for low-core-count machines

def main():
    # --- Configuration ---
    config = {
        'index_file': r'E:/Qinchaojun/AIDD-TRAIN/index/INDEX_general_PL.2020R1.lst',
        'dataset_path': r'E:/Qinchaojun/PDBbind-2025.8.4/P-L/',
        'processed_data_dir': r'E:/Qinchaojun/AIDD-TRAIN/processed_data',
        'epochs': 20,
        'batch_size': 16, # Reduced batch size to lower memory usage during training
        'learning_rate': 0.001,
        'train_split': 0.8,
        'num_workers': get_optimal_workers(),
    }
    # --- End Configuration ---

    os.makedirs(config['processed_data_dir'], exist_ok=True)

    print("Step 1: Parsing PDBbind index...")
    pdb_info = get_pdb_info(config['index_file'])
    print(f"-> Found {len(pdb_info)} total entries in index.")
    all_data_paths = get_data_paths(pdb_info, config['dataset_path'])
    print(f"-> Found {len(all_data_paths)} pairs with existing data files.")

    print()
    print("Step 2: Creating/loading processed PyTorch Geometric dataset...")
    dataset = PDBBindDataset(
        root=config['processed_data_dir'], 
        data_paths=all_data_paths, 
        num_workers=config['num_workers']
    )
    print(f"-> Dataset loaded with {len(dataset)} samples.")
    
    if len(dataset) == 0:
        print("FATAL: Dataset is empty. Cannot proceed.")
        exit()

    print()
    print("Step 3: Splitting data and creating loaders...")
    dataset = dataset.shuffle()
    train_size = int(config['train_split'] * len(dataset))
    train_dataset = dataset[:train_size]
    val_dataset = dataset[train_size:]
    
    # On Windows, DataLoader num_workers must often be 0.
    loader_num_workers = 0 if platform.system() == 'Windows' else 2

    train_loader = DataLoader(train_dataset, batch_size=config['batch_size'], shuffle=True, num_workers=loader_num_workers)
    val_loader = DataLoader(val_dataset, batch_size=config['batch_size'], shuffle=False, num_workers=loader_num_workers)
    print(f"-> Train: {len(train_dataset)} | Validation: {len(val_dataset)}")

    print()
    print("Step 4: Setting up model and optimizer...")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"-> Using device: {device}")
        
    # Get a sample to determine model input dimensions.
    # This requires at least one processed data sample to exist.
    try:
        # Find the first valid sample to inspect dimensions
        first_sample = None
        for i in range(len(dataset)):
            sample = dataset.get(i)
            if sample is not None:
                first_sample = sample
                break
        
        if first_sample is None:
            raise FileNotFoundError("No valid processed samples found.")

        protein_in_dim = first_sample.protein_x.shape[1]
        ligand_in_dim = first_sample.ligand_x.shape[1]
    except (FileNotFoundError, IndexError):
        print("FATAL: Cannot determine model dimensions. No processed data found.")
        print("Please ensure your data paths are correct and at least one sample can be processed.")
        exit()


    model = Net(protein_in_dim=protein_in_dim, ligand_in_dim=ligand_in_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'])
    print(f"-> Model initialized with {sum(p.numel() for p in model.parameters()):,} parameters.")

    print()
    print("Step 5: Starting training...")
    for epoch in range(1, config['epochs'] + 1):
        train_loss = train(model, train_loader, optimizer, device)
        val_loss = test(model, val_loader, device)
        print(f"Epoch {epoch:02d}/{config['epochs']} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

    print()
    print("--- Training Finished ---")

# ==============================================================================
# SCRIPT ENTRY POINT
# ==============================================================================

if __name__ == '__main__':
    # This setup is crucial for multiprocessing to work correctly and safely 
    # on platforms like Windows and macOS. It must be inside this block.
    try:
        set_start_method('spawn')
    except RuntimeError:
        # This error is raised if the start method is already set, which is fine.
        pass
    
    main()
