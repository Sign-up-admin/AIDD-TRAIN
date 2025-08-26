import os
import logging
import torch
from rdkit import Chem

from .features import get_atom_features


def get_ligand_graph(mol, pdb_code="N/A"):
    """
    Converts an RDKit molecule object into a graph representation for the ligand.
    Removes duplicate atoms based on their 3D coordinates.

    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule object for the ligand.
        pdb_code (str): The PDB code for logging purposes.

    Returns:
        dict or None: A dictionary containing feature and position tensors ('v', 'pos')
                      and a set of seen positions ('seen_pos'), or None if the graph
                      cannot be generated.
    """
    if mol is None or mol.GetNumConformers() == 0:
        return None

    conformer = mol.GetConformer()
    atom_features_list, positions_list = [], []
    seen_positions = set()

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        pos = conformer.GetAtomPosition(atom_idx)
        pos_tuple = (round(pos.x, 4), round(pos.y, 4), round(pos.z, 4))

        if pos_tuple in seen_positions:
            logging.warning(f"Duplicate ligand atom position found and removed in {pdb_code}: {pos_tuple}")
            continue
        
        seen_positions.add(pos_tuple)
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
    """
    Converts a protein PDB file into a graph representation.
    Removes heteroatoms, atoms overlapping with the ligand, and duplicate atoms.

    Args:
        protein_path (str): Path to the protein's PDB file.
        ligand_positions_set (set): A set of 3D coordinates of ligand atoms to check for overlaps.

    Returns:
        dict or None: A dictionary containing feature and position tensors ('v', 'pos'),
                      or None if the graph cannot be generated.
    """
    try:
        mol = Chem.MolFromPDBFile(protein_path, removeHs=False, sanitize=True)
    except Exception as e:
        logging.error(f"Failed to load or sanitize PDB file {protein_path}: {e}")
        return None

    if mol is None or mol.GetNumConformers() == 0:
        if mol is None:
            logging.warning(f"RDKit could not process PDB file: {protein_path}. Mol object is None.")
        else:
            logging.warning(f"No conformers found in PDB file: {protein_path}.")
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
        atom_features_list.append(get_atom_features(atom))
        positions_list.append(list(pos_tuple))

    if not atom_features_list:
        return None

    return {
        'v': torch.tensor(atom_features_list, dtype=torch.float),
        'pos': torch.tensor(positions_list, dtype=torch.float)
    }
