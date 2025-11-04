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

mols = []
for prediction in chemical_structure.region_predictions[1]:
    mol = structure_prediction.smiles_to_objects(prediction.smiles)
    mols.append(mol)

# Drawing molecules dynamically
# https://iwatobipen.wordpress.com/2020/01/17/draw-rdkit-mol-reaction-object-on-html-without-static-png-image-rdkit-memo/
from IPython.display import HTML, display

import io
import base64
from PIL import Image

from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

img = Draw.MolToImage(mols[0])
drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)

drawer.SetFontSize(1.0)
drawer.DrawMolecules(mols)
rdMolDraw2D.PrepareAndDrawMolecule(drawer, mols[0])
drawer.FinishDrawing()
bio = io.BytesIO()

text = drawer.GetDrawingText()
text
imtext = base64.b64encode(text).decode("utf8")

display(HTML(f"<div><img src='data:image/png;base64, {text}' alt='hoge'/></div>").data)
