import re
import math

def parse_binding_data(binding_str):
    """
    Parses a string to extract and convert binding affinity data to pK_d/i value.

    Args:
        binding_str (str): The string containing binding data (e.g., 'Kd=12.3nM').

    Returns:
        float or None: The calculated -log10 molar value, or None if parsing fails.
    """
    match = re.search(r"([0-9.]+)([a-zA-Z]+)", binding_str)
    if not match:
        return None

    value, unit = float(match.group(1)), match.group(2)
    unit_map = {'mM': 1e-3, 'uM': 1e-6, 'nM': 1e-9, 'pM': 1e-12, 'fM': 1e-15}
    molar_value = value * unit_map.get(unit, 1.0)
    return -math.log10(molar_value) if molar_value > 0 else 0.0
