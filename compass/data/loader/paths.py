import os
from tqdm import tqdm

def get_pdb_info(index_file):
    """
    Parses the PDBbind index file to extract PDB code, year, and binding data.

    Args:
        index_file (str): Path to the index file.

    Returns:
        dict: A dictionary mapping PDB codes to their year and binding data string.
    """
    pdb_info = {}
    with open(index_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    pdb_info[parts[0]] = {'year': int(parts[2]), 'binding_data': parts[3]}
                except (ValueError, IndexError):
                    pass  # Ignore malformed lines
    return pdb_info

def get_data_paths(pdb_info, dataset_path):
    """
    Constructs and verifies paths to protein and ligand files for each PDB entry.

    Args:
        pdb_info (dict): A dictionary with PDB metadata from get_pdb_info.
        dataset_path (str): The root path of the PDBbind dataset.

    Returns:
        list: A list of dictionaries, each containing paths and metadata for a valid entry.
    """
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
