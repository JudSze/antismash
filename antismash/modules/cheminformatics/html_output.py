"""
HTML output generation module for antiSMASH

This module generates HTML output for displaying chemical structures
in a custom tab alongside Tailoring and Gene overview tabs.
"""

from typing import List

from rdkit import Chem
from rdkit.Chem import Draw

from antismash.common.secmet import Record
from antismash.common import secmet

from antismash.modules import nrps_pks
from antismash.modules.nrps_pks import specific_analysis
from antismash.modules.nrps_pks.results import NRPS_PKS_Results
from antismash.modules.cheminformatics import structure_prediction

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RegionLayer


def will_handle(products: List[str]) -> bool:
    """
    Returns true if the chemical structure prediction is enabled
    a placeholder for now
    """
    return True


def generate_html(region_layer: RegionLayer, results) -> HTMLSections:
    """
    """
    html = HTMLSections("cheminformatics")

    # Generate the chemical structure visualization
    structure_html = generate_structure_section(region_layer, results)

    # Add as a sidepanel detail (appears in tabs next to Tailoring/Gene overview)
    html.add_sidepanel_section(
        "Structure",  # Tab name
        structure_html,
        class_name="structure-section"
    )

    return html


def generate_structure_section(region_layer: RegionLayer, results: NRPS_PKS_Results) -> Markup:
    """
    """
    # TODO: Extract structure data from results
    # structure_data = results.get_structure_data()

    # Generate SVG using RDKit
    svg_data = generate_structure_svg(results)

    # Create HTML template
    template = FileTemplate(path.get_full_path(__file__, "templates", "structure.html"))

    # TODO: Prepare template variables
    template_vars = {
        "region": region_layer,
        "structure_svg": svg_data,
    }

    return template.render(**template_vars)


def generate_structure_svg(nrps_pks_results: NRPS_PKS_Results) -> str:
    """
    """
    structures = []
    visualizations = []

    for prediction in nrps_pks_results.region_predictions[1]:
        structure = structure_prediction.smiles_to_objects(prediction.smiles)
        structures.append(structure)

    for molecule in structures:
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(400, 300)
        drawer.DrawMolecule(molecule)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        visualizations.append(svg)

    return visualizations


# Optional: Add a details section (appears below the cluster visualization)
def generate_details_div(region_layer: RegionLayer, results) -> Markup:
    """
    """
    template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))

    template_vars = {
        "region": region_layer,
        "results": results,
    }

    return template.render(**template_vars)