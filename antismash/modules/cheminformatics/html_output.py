from rdkit.Chem import Draw

from antismash.common import path

from antismash.common.html_renderer import FileTemplate, HTMLSection

from antismash.modules.nrps_pks.results import NRPS_PKS_Results

def generate_html(results: NRPS_PKS_Results):
    html = HTMLSection("chemical-structure")
    structure_viz = Draw.MolsToImage(results.chemical_structure)

    details_template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    details = details_template.render(structure_viz)

    html.add_setail_section("Chemical Structure", details, class_name="chemical-structure")
    return structure_viz