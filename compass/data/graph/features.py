ELEMENTS = ['C', 'O', 'N', 'S', 'P', 'H', 'F', 'Cl', 'Br', 'I', 'UNK']
ELEMENT_MAP = {el: i for i, el in enumerate(ELEMENTS)}

def get_atom_features(atom):
    """
    Generates a one-hot encoded feature vector for a given atom.

    Args:
        atom (rdkit.Chem.Atom): The RDKit atom object.

    Returns:
        list: A list representing the one-hot encoded feature vector.
    """
    feat = [0] * len(ELEMENTS)
    element = atom.GetSymbol()
    feat[ELEMENT_MAP.get(element, len(ELEMENTS) - 1)] = 1
    return feat
