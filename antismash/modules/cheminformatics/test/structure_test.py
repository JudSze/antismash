from antismash.main import get_all_modules

from antismash.common import secmet

from antismash.config import build_config
from antismash.modules import nrps_pks
from antismash.modules.nrps_pks import specific_analysis

from antismash.modules.cheminformatics import structure_prediction

genbank="/home/szenei/antismash/BGC0000034.gbk"

config = build_config(['--help'])
config

args = ["--enable-nrps-pks",
        "--genefinding-tool", "prodigal",
        "--cb-general",
        "--pfam2go",
        "--smcog-trees",
        "--asf",
        "--rre"]

options = build_config(args, isolated=True,
                       modules=get_all_modules())

gbk_record = secmet.record.Record.from_genbank(genbank)[0]
domain_res = nrps_pks.run_on_record(gbk_record, None, options)

chemical_structure = specific_analysis(gbk_record, domain_res, options)

for prediction in chemical_structure.region_predictions[1]:
    structure_prediction.smiles_to_objects(prediction.smiles)

structure_prediction.smiles_to_objects(chemical_structure.consensus)


for region in gbk_record.get_regions():
    for candidate_cluster in region._candidate_clusters:
        print(candidate_cluster.polymer)
