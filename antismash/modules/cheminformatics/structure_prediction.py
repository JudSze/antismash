from rdkit import Chem


class Structures:
    def __init__(self, smiles):
        self.smiles = str
        self.mol = Chem.MolFromSmiles(smiles)

def smiles_to_objects(smiles: str):
    mol_structure = Chem.MolFromSmiles(smiles)
    return mol_structure

