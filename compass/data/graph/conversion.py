import os
import logging
import torch
from rdkit import Chem

from .features import get_atom_features


def get_ligand_graph(mol, pdb_code="N/A"):
    """
    Converts an RDKit molecule object into a graph representation for the ligand.
    Removes duplicate atoms and checks for close contacts.

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
    
    # Step 1: Extract unique atoms and their features
    atom_info = {}
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        pos_tuple = (round(pos.x, 4), round(pos.y, 4), round(pos.z, 4))
        if pos_tuple not in atom_info:
            atom_info[pos_tuple] = get_atom_features(atom)
        else:
            logging.warning(f"Duplicate ligand atom position found and removed in {pdb_code}: {pos_tuple}")

    if not atom_info:
        return None

    positions_list = list(atom_info.keys())
    atom_features_list = list(atom_info.values())
    
    # Step 2: Atom proximity alert
    if len(positions_list) > 1:
        pos_tensor = torch.tensor(positions_list, dtype=torch.float)
        dists = torch.pdist(pos_tensor)
        if torch.any(dists < 0.5):
            min_dist = torch.min(dists)
            logging.warning(
                f"Atom proximity alert in {pdb_code} (ligand): "
                f"Minimum distance between atoms is {min_dist:.4f} Å, which is below the 0.5 Å threshold."
            )

    return {
        'v': torch.tensor(atom_features_list, dtype=torch.float),
        'pos': torch.tensor(positions_list, dtype=torch.float),
        'seen_pos': set(positions_list)
    }


def get_protein_graph(protein_path, ligand_positions_set):
    """
    Converts a protein PDB file into a graph representation.
    Removes heteroatoms, atoms overlapping with the ligand, and duplicate atoms.
    Checks for close contacts.

    Args:
        protein_path (str): Path to the protein's PDB file.
        ligand_positions_set (set): A set of 3D coordinates of ligand atoms to check for overlaps.

    Returns:
        dict or None: A dictionary containing feature and position tensors ('v', 'pos'),
                      or None if the graph cannot be generated.
    """
    pdb_code = os.path.basename(protein_path).split('_')[0]
    try:
        mol = Chem.MolFromPDBFile(protein_path, removeHs=False, sanitize=True)
    except Exception as e:
        logging.error(f"Failed to load or sanitize PDB file {protein_path} for {pdb_code}: {e}")
        return None

    if mol is None or mol.GetNumConformers() == 0:
        logging.warning(f"RDKit could not process or find conformers in PDB file: {protein_path}")
        return None

    conformer = mol.GetConformer()

    # Step 1: Extract unique atoms and their features, filtering heteroatoms and overlaps
    atom_info = {}
    for atom in mol.GetAtoms():
        if atom.GetPDBResidueInfo().GetIsHeteroAtom():
            continue

        pos = conformer.GetAtomPosition(atom.GetIdx())
        pos_tuple = (round(pos.x, 4), round(pos.y, 4), round(pos.z, 4))

        if pos_tuple in ligand_positions_set:
            logging.warning(f"Protein-Ligand Overlap: Protein atom at {pos_tuple} in {pdb_code} removed.")
            continue
        
        if pos_tuple not in atom_info:
            atom_info[pos_tuple] = get_atom_features(atom)
        else:
            logging.warning(f"Intra-Protein Duplicate: Protein atom at {pos_tuple} in {pdb_code} removed.")

    if not atom_info:
        return None
        
    positions_list = list(atom_info.keys())
    atom_features_list = list(atom_info.values())

    # Step 2: Atom proximity alert
    if len(positions_list) > 1:
        pos_tensor = torch.tensor(positions_list, dtype=torch.float)
        dists = torch.pdist(pos_tensor)
        if torch.any(dists < 0.5):
            min_dist = torch.min(dists)
            logging.warning(
                f"Atom proximity alert in {pdb_code} (protein): "
                f"Minimum distance between atoms is {min_dist:.4f} Å, which is below the 0.5 Å threshold."
            )

    return {
        'v': torch.tensor(atom_features_list, dtype=torch.float),
        'pos': torch.tensor(positions_list, dtype=torch.float)
    }


def check_graph_connectivity(pos_tensor, cutoff, pdb_code="N/A"):
    """
    Checks for isolated nodes in a graph.

    Args:
        pos_tensor (torch.Tensor): A tensor of atom positions (N, 3).
        cutoff (float): The cutoff radius for neighbor search.
        pdb_code (str): The PDB code for logging.
    """
    if pos_tensor.shape[0] < 2:
        return # Cannot check connectivity for a single node or empty graph

    # Compute pairwise distances
    dists = torch.pdist(pos_tensor)
    
    # Create adjacency matrix based on cutoff
    adj_matrix = torch.zeros((pos_tensor.shape[0], pos_tensor.shape[0]), dtype=torch.bool)
    indices = torch.triu_indices(pos_tensor.shape[0], pos_tensor.shape[0], offset=1)
    adj_matrix[indices[0], indices[1]] = dists < cutoff
    adj_matrix = adj_matrix | adj_matrix.T

    # Check for isolated nodes
    num_neighbors = adj_matrix.sum(dim=1)
    isolated_nodes = torch.where(num_neighbors == 0)[0]

    if len(isolated_nodes) > 0:
        logging.warning(
            f"Graph connectivity issue in {pdb_code}: Found {len(isolated_nodes)} isolated atom(s) "
            f"with no neighbors within the {cutoff} Å cutoff."
        )
