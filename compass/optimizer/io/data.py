import os
import logging

from compass.data.processing import process_item
from compass.config import CONFIG
from compass.optimizer.config import TEST_SAMPLE_CONFIG

logger = logging.getLogger("HardwareOptimizer")

def prepare_real_data(project_root_path):
    sample_name = TEST_SAMPLE_CONFIG['pdb_code']
    logger.info(f"--- Preparing Real-World Test Case: '{sample_name}' ---")

    # --- WORKAROUND for path issue ---
    protein_path = os.path.join(project_root_path, 'compass', 'optimizer', '1jmf', '1jmf_pocket.pdb')
    ligand_path = os.path.join(project_root_path, 'compass', 'optimizer', '1jmf', '1jmf_ligand.sdf')
    item = {
        'pdb_code': sample_name,
        'protein_path': protein_path,
        'ligand_path': ligand_path,
        'binding_data': TEST_SAMPLE_CONFIG['binding_data']
    }

    max_atoms_for_sample = TEST_SAMPLE_CONFIG['max_atoms']
    logger.info(f"Temporarily setting max_atoms to {max_atoms_for_sample} to process {sample_name}...")
    # TODO: Refactor `process_item` to accept a config dictionary directly,
    # which would remove this temporary and fragile patching of the global CONFIG.
    # For now, this maintains existing behavior.
    original_max_atoms = CONFIG.get('max_atoms')
    processed_test_data = None
    try:
        CONFIG['max_atoms'] = max_atoms_for_sample
        processed_test_data = process_item(item)
    finally:
        if original_max_atoms is not None:
            CONFIG['max_atoms'] = original_max_atoms
        elif 'max_atoms' in CONFIG:
            del CONFIG['max_atoms']

    if processed_test_data is None:
        logger.critical(f"Failed to process {sample_name} data. Aborting."); exit()

    logger.info(f"--- Successfully processed {sample_name}. Atom count: {processed_test_data.num_nodes} ---")
    return processed_test_data