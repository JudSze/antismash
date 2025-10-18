from typing import List

from rdkit import Chem

from antismash.modules.nrps_pks.results import CandidateClusterPrediction

def smiles_to_objects(predicted_clusters=List[CandidateClusterPrediction]):
    mol_structures = {}

    for predicted_cluster in predicted_clusters:
        mol_structure = Chem.MolFromSmiles(predicted_cluster.smiles)
        mol_structures[predicted_cluster.polymer]= mol_structure

    return mol_structures
