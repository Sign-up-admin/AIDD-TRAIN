import os
import re
import math
import logging

import torch
import numpy as np
from torch_geometric.data import Data
from torch_geometric.utils import add_self_loops
from rdkit import Chem, RDLogger
from Bio.PDB import PDBParser
from scipy.spatial import KDTree
from tqdm import tqdm

# Suppress RDKit warnings about 2D/3D coordinates
RDLogger.logger().setLevel(RDLogger.CRITICAL)

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

def protein_structure_to_graph(structure, k=32, cutoff=10.0):
    atoms_to_include = [atom for residue in structure.get_residues() if residue.id[0] == ' ' for atom in residue.get_atoms()]
    
    num_atoms = len(atoms_to_include)
    if num_atoms < k + 1:
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

    kdtree = KDTree(positions_np)
    distances, indices = kdtree.query(positions_np, k=k+1, p=2.0)

    edge_set = set()
    for i in range(num_atoms):
        for j_idx in range(1, k + 1):
            neighbor_idx = indices[i, j_idx]
            if distances[i, j_idx] < cutoff:
                pair = tuple(sorted((i, neighbor_idx)))
                edge_set.add(pair)
    
    if not edge_set:
        return None

    rows, cols = zip(*edge_set)
    rows = list(rows)
    cols = list(cols)
    
    edge_index = torch.tensor([rows + cols, cols + rows], dtype=torch.long)
    
    pos = torch.from_numpy(positions_np).float()
    data = Data(x=atom_features, edge_index=edge_index, pos=pos)
    data.edge_index, _ = add_self_loops(data.edge_index, num_nodes=data.num_nodes)
    return data

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

        if count_pdb_atoms(protein_path) > 15000: # Relaxed limit
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
        
        # Combine protein and ligand graphs into a single Data object
        num_protein_nodes = protein_graph.x.size(0)
        
        x = torch.cat([protein_graph.x, ligand_graph.x], dim=0)
        pos = torch.cat([protein_graph.pos, ligand_graph.pos], dim=0)
        
        # Adjust ligand edge indices and combine
        ligand_edge_index_adjusted = ligand_graph.edge_index + num_protein_nodes
        edge_index = torch.cat([protein_graph.edge_index, ligand_edge_index_adjusted], dim=1)

        return Data(
            x=x,
            pos=pos,
            edge_index=edge_index,
            y=torch.tensor([y], dtype=torch.float),
            pdb_code=pdb_code
        )
    except Exception as e:
        logging.error(f"ERROR for PDB {pdb_code}: {e}")
        return None
