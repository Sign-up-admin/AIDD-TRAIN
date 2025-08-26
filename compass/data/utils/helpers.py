def count_pdb_atoms(pdb_path):
    """
    Counts the number of ATOM records in a PDB file.

    Args:
        pdb_path (str): The path to the PDB file.

    Returns:
        int: The number of atoms, or 0 if the file cannot be read.
    """
    try:
        with open(pdb_path, 'r') as f:
            return sum(1 for line in f if line.startswith('ATOM'))
    except IOError:
        return 0
