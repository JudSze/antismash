# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Replace this text with a description of the module.
    It can also include references for the method implemented.
"""

# start with standard library imports
import logging
from typing import Any, Dict, List, Optional

# then any imports from external modules, e.g. biopython, if relevant

# then any imports from antismash
from antismash.modules.nrps_pks.results import NRPS_PKS_Results
from antismash.common.secmet import Record
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs

# then any local file imports, e.g. from .somefile import..., if relevant

NAME = "chemical-structure"
SHORT_DESCRIPTION = "Chemical structure prediction"

def get_arguments() -> ModuleArgs:
    """ Builds any commandline argument constructs that may be required

        Returns:
            an empty or populated ModuleArgs instance
    """
    # construct the argument group, with section and prefix
    # the prefix will be enforced for all command line options for the module
    args = ModuleArgs('Additional analysis', 'chemicalstructure', enabled_by_default=True)

    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks that the provided options are compatible with each other

        Arguments:
            options: the current antismash config object

        Returns:
            a list of strings, each string being an issue with the given options
    """
    issues = []
    return issues


def check_prereqs(options: ConfigType) -> List[str]:
    """ Check that all prerequisites are present

        Arguments:
            options: the current antismash config object

        Returns:
            a list of strings, each string being an issue with prerequisites
    """
    # behaves similarly to check_options(), though checking for built databases,
    # that external programs are available, and so on
    # see antismash.detection.hmm_detection for an example of these

    # if there are no external prerequisites, this can just return an empty list
    return []


def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the module is enabled with the options provided
    """
    # the logic here depends on which command options you've created
    # using the example above, this is as simple as returning the toggle
    return options


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType):
    """ Regenerate the previous results from JSON format.

        Arguments:
            previous: the previous results as from JSON
            record: the Record these previous results were originally created from
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation,
            or None if the current options require the analysis to be rerun or cannot be regenerated
    """
    # if there isn't anything to work with, just return None
    if not previous:
        return None
    return NRPS_PKS_Results.from_json(previous, record)


def run_on_record(record: Record, results: NRPS_PKS_Results, options: ConfigType) -> NRPS_PKS_Results:
    """ Run the analysis, unless the previous results apply to the given record

        Arguments:
            record: the Record being analysed
            results: an existing instance of the module's ModuleResults implementation (or None)
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation
    """
    # after a safety check that the results are the correct ones for the record, return them
    if isinstance(results, NRPS_PKS_Results) and results.record_id == record.id:
        return results
    # otherwise run the actual analysis and generate a results instance with your analysis results
    results = NRPS_PKS_Results(record.id)
    # and return it
    return results