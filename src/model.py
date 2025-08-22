import torch
from torch_geometric.nn.models import ViSNet

# Mapping from element index to atomic number (0 for UNK)
# Based on this list from the data processing section:
# ELEMENTS = ['C', 'O', 'N', 'S', 'P', 'H', 'F', 'Cl', 'Br', 'I', 'UNK']
INDEX_TO_Z = torch.tensor([6, 8, 7, 16, 15, 1, 9, 17, 35, 53, 0], dtype=torch.long)

class ViSNetPDB(torch.nn.Module):
    """
    A wrapper for the ViSNet model to handle the combined protein-ligand graph data.
    """
    def __init__(self, hidden_channels=128, num_layers=6, num_rbf=64, cutoff=8.0, max_num_neighbors=32):
        super().__init__()
        self.visnet = ViSNet(
            hidden_channels=hidden_channels,
            num_layers=num_layers,
            num_rbf=num_rbf,
            cutoff=cutoff,
            max_num_neighbors=max_num_neighbors,
            lmax=1,
            vecnorm_type='rms', # Changed from 'max_min' for better stability
            trainable_vecnorm=True, # Enabled for adaptive normalization
            num_heads=8,
            trainable_rbf=False,
            max_z=100, # Max atomic number in periodic table
            reduce_op='add'
        )
        self.register_buffer('index_to_z', INDEX_TO_Z)

    def forward(self, data):
        # Convert one-hot encoded features to atomic numbers (z)
        z_idx = data.x.argmax(dim=-1)
        z = self.index_to_z[z_idx]
        
        # Run ViSNet forward pass
        output = self.visnet(z, data.pos, data.batch)
        if isinstance(output, tuple):
            return output[0]
        return output
