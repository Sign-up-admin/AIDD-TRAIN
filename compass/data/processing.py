import os
import logging
from pathlib import Path

import torch
import numpy as np
from torch_geometric.data import Data
from rdkit import Chem, RDLogger

from compass.config import CONFIG
from compass.data.graph.conversion import get_ligand_graph, get_protein_graph
from compass.data.loader.parser import parse_binding_data
from compass.data.utils.helpers import count_pdb_atoms

# Suppress RDKit warnings
RDLogger.logger().setLevel(RDLogger.CRITICAL)

# Configure logging to file - use relative path or environment variable
log_directory = os.getenv('COMPASS_LOG_DIR', 'logs')
log_directory = Path(log_directory).resolve()
log_directory.mkdir(parents=True, exist_ok=True)
log_file = log_directory / "processing.log"

# Clear existing handlers to ensure no console output for warnings
# and configure new handler for file logging.
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename=str(log_file),
    filemode='a'
)

def process_item(item):
    """
    Processes a single data item (PDB entry) to create a graph-based data object.

    This function performs the following steps:
    1. Validates file paths and sizes.
    2. Loads ligand and protein structures.
    3. Converts structures to graph representations, handling duplicates and overlaps.
    4. Parses binding affinity data.
    5. Assembles the final Data object for use in PyTorch Geometric.
    6. Performs a final validation and repair step to ensure data integrity.

    Args:
        item (dict): A dictionary containing information about the PDB entry, including
                     paths to protein and ligand files, and binding data.

    Returns:
        torch_geometric.data.Data or None: A processed Data object ready for model training,
                                           or None if processing fails at any stage.
    """
    pdb_code = item.get('pdb_code', 'N/A')
    max_atoms = CONFIG.get('max_atoms', 10000)
    try:
        protein_path = os.path.normpath(item['protein_path'])
        ligand_path = os.path.normpath(item['ligand_path'])

        if os.path.getsize(protein_path) > 100 * 1024 * 1024:  # 100MB limit
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
        # This step is crucial as duplicates can arise from rounding and merging.
        _, unique_indices = np.unique(data.pos.numpy(), axis=0, return_index=True)

        if len(unique_indices) < data.pos.shape[0]:
            original_atom_count = data.pos.shape[0]
            removed_count = original_atom_count - len(unique_indices)
            logging.warning(f"Found and removed {removed_count} duplicate atom(s) in {pdb_code}.")

            # Sort indices to maintain the original relative order of atoms as much as possible.
            unique_indices = np.sort(unique_indices)
            data.x = data.x[unique_indices]
            data.pos = data.pos[unique_indices]

        return data

    except Exception as e:
        logging.error(f"ERROR processing PDB {pdb_code}: {e}", exc_info=True)
        return None
