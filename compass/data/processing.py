import os
import re
import math
import logging

import torch
import numpy as np
from torch_geometric.data import Data
from rdkit import Chem, RDLogger
from tqdm import tqdm

from compass.config import CONFIG

# Suppress RDKit warnings
RDLogger.logger().setLevel(RDLogger.CRITICAL)

# =============================================================================
# PART 1: DATA PATH AND BINDING VALUE PARSING
# =============================================================================

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
                'binding_data': info['binding_data'],
                'year_dir': year_dir
            })
    return data_paths

def parse_binding_data(binding_str):
    match = re.search(r"([0-9.]+)([a-zA-Z]+)", binding_str)
    if not match: return None
    
    value, unit = float(match.group(1)), match.group(2)
    unit_map = {'mM': 1e-3, 'uM': 1e-6, 'nM': 1e-9, 'pM': 1e-12, 'fM': 1e-15}
    molar_value = value * unit_map.get(unit, 1.0)
    return -math.log10(molar_value) if molar_value > 0 else 0.0

# =============================================================================
# PART 2: MOLECULE-TO-GRAPH CONVERSION (WITH INTER-MOLECULE DUPLICATE REMOVAL)
# =============================================================================

ELEMENTS = ['C', 'O', 'N', 'S', 'P', 'H', 'F', 'Cl', 'Br', 'I', 'UNK']
ELEMENT_MAP = {el: i for i, el in enumerate(ELEMENTS)}

def get_atom_features(atom):
    feat = [0] * len(ELEMENTS)
    element = atom.GetSymbol()
    feat[ELEMENT_MAP.get(element, len(ELEMENTS) - 1)] = 1
    return feat

def get_ligand_graph(mol, pdb_code="N/A"):
    if mol is None or mol.GetNumConformers() == 0:
        return None

    conformer = mol.GetConformer()
    atom_features_list, positions_list = [], []
    seen_positions = set()
    old_to_new_idx = {}
    new_idx_counter = 0

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        pos = conformer.GetAtomPosition(atom_idx)
        pos_tuple = (round(pos.x, 4), round(pos.y, 4), round(pos.z, 4))

        if pos_tuple in seen_positions:
            logging.warning(f"Duplicate ligand atom position found and removed in {pdb_code}: {pos_tuple}")
            continue
        
        seen_positions.add(pos_tuple)
        old_to_new_idx[atom_idx] = new_idx_counter
        new_idx_counter += 1
        atom_features_list.append(get_atom_features(atom))
        positions_list.append(list(pos_tuple))

    if not atom_features_list:
        return None

    return {
        'v': torch.tensor(atom_features_list, dtype=torch.float),
        'pos': torch.tensor(positions_list, dtype=torch.float),
        'seen_pos': seen_positions
    }

def get_protein_graph(protein_path, ligand_positions_set):
    mol = Chem.MolFromPDBFile(protein_path, removeHs=False, sanitize=True)
    if mol is None or mol.GetNumConformers() == 0:
        return None

    atom_features_list, positions_list = [], []
    pdb_code = os.path.basename(protein_path).split('_')[0]
    protein_positions_seen = set()
    conformer = mol.GetConformer()

    for atom in mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info and pdb_info.GetIsHeteroAtom():
            continue

        pos = conformer.GetAtomPosition(atom.GetIdx())
        pos_tuple = (round(pos.x, 4), round(pos.y, 4), round(pos.z, 4))

        if pos_tuple in ligand_positions_set:
            logging.warning(f"Protein-Ligand Overlap: Protein atom at {pos_tuple} in {pdb_code} overlaps with ligand and will be removed.")
            continue
        
        if pos_tuple in protein_positions_seen:
            logging.warning(f"Intra-Protein Duplicate: Duplicate protein atom at {pos_tuple} in {pdb_code} will be removed.")
            continue
        
        protein_positions_seen.add(pos_tuple)
        element = atom.GetSymbol().strip().upper()
        feat = [0] * len(ELEMENTS)
        feat[ELEMENT_MAP.get(element, len(ELEMENTS) - 1)] = 1
        atom_features_list.append(feat)
        positions_list.append(list(pos_tuple))

    if not atom_features_list:
        return None

    return {
        'v': torch.tensor(atom_features_list, dtype=torch.float),
        'pos': torch.tensor(np.array(positions_list), dtype=torch.float)
    }

def count_pdb_atoms(pdb_path):
    try:
        with open(pdb_path, 'r') as f:
            return sum(1 for line in f if line.startswith('ATOM'))
    except IOError:
        return 0

def process_item(item):
    pdb_code = item.get('pdb_code', 'N/A')
    max_atoms = CONFIG.get('max_atoms', 10000)
    try:
        # --- Path Separator Fix for Cross-Platform Compatibility ---
        # On Windows, paths constructed with mixed separators (e.g., 'dir\\subdir/file')
        # can fail. This normalizes the separators to the OS-specific one.
        protein_path = os.path.normpath(item['protein_path'])
        ligand_path = os.path.normpath(item['ligand_path'])

        if os.path.getsize(protein_path) > 100 * 1024 * 1024: # 100MB limit
            logging.warning(f"Skipping {pdb_code}: Protein file is very large (>100MB).")
            return None

        atom_count = count_pdb_atoms(protein_path)
        if atom_count > max_atoms:
            logging.warning(f"Skipping {pdb_code}: Protein has too many atoms ({atom_count} > {max_atoms}).")
            return None

        ligand = None
        if ligand_path.endswith('.mol2'):
            ligand = Chem.MolFromMol2File(ligand_path, removeHs=False, sanitize=True)
        elif ligand_path.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(ligand_path, removeHs=False, sanitize=True)
            ligand = next(suppl, None)
        
        if ligand is None or ligand.GetNumAtoms() == 0:
            logging.warning(f"Skipping {pdb_code}: Could not load ligand or ligand has no atoms.")
            return None

        ligand_graph = get_ligand_graph(ligand, pdb_code)
        if ligand_graph is None:
            logging.warning(f"Skipping {pdb_code}: Ligand graph could not be generated.")
            return None

        protein_graph = get_protein_graph(protein_path, ligand_positions_set=ligand_graph['seen_pos'])
        if protein_graph is None:
            logging.warning(f"Skipping {pdb_code}: Protein graph could not be generated.")
            return None

        y = parse_binding_data(item['binding_data'])
        if y is None: 
            logging.warning(f"Skipping {pdb_code}: Could not parse binding data.")
            return None
        
        data = Data(y=torch.tensor([y], dtype=torch.float), pdb_code=pdb_code)
        data.x = torch.cat([ligand_graph['v'], protein_graph['v']], dim=0)
        data.pos = torch.cat([ligand_graph['pos'], protein_graph['pos']], dim=0)

        # Final check for duplicate positions before returning the data object.
        # This is a critical safeguard against data corruption.
        final_pos_tensor = data.pos
        
        # We use numpy's unique function here because it provides the indices of the unique elements,
        # which is essential for keeping the 'pos' and 'x' tensors synchronized.
        _, unique_indices = np.unique(final_pos_tensor.numpy(), axis=0, return_index=True)

        # If the number of unique indices is less than the total number of atoms, we have duplicates.
        if len(unique_indices) < final_pos_tensor.shape[0]:
            original_atom_count = final_pos_tensor.shape[0]
            logging.warning(f"Duplicate atom positions detected in {pdb_code}. Attempting automatic repair...")

            # --- Attempt Repair ---
            # We filter both the position and feature tensors using the indices of the unique atoms.
            # This ensures that the data remains consistent.
            data.pos = final_pos_tensor[unique_indices]
            data.x = data.x[unique_indices]
            
            # --- Final Verification ---
            # After attempting the repair, we perform one last check to ensure the data is clean.
            # This is the "trust but verify" principle in action.
            final_pos_tensor_after_fix = data.pos
            if torch.unique(final_pos_tensor_after_fix, dim=0).shape[0] < final_pos_tensor_after_fix.shape[0]:
                # This is our final, critical safeguard. If this check fails, it means our repair
                # logic has a flaw, and we must stop to avoid corrupting the dataset.
                logging.error(f"CRITICAL ERROR in {pdb_code}: Automatic repair of duplicate positions FAILED. "
                              f"The data for this entry will be discarded.")
                return None
            
            removed_count = original_atom_count - len(unique_indices)
            logging.warning(f"Successfully repaired {pdb_code}. Removed {removed_count} duplicate atom(s).")

        return data

    except Exception as e:
        logging.error(f"ERROR processing PDB {pdb_code}: {e}", exc_info=True)
        return None
