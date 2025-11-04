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
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs

# then any local file imports, e.g. from .somefile import..., if relevant

NAME = "the_name_of_your_module"
SHORT_DESCRIPTION = "a short description of the module"


# define a results class, this is important as adding information to the record
# during analysis will cause issues
# for detailed examples, see any of the other analysis modules' implementations
class TemplateAnalysisResults(ModuleResults):
    """ Example results class for the analysis module template """
    schema_version = 1  # when the data format in the results changes, this needs to be incremented

    # define whatever construction arguments are needed, record_id is required by the superclass
    # it's good to keep any command line option values here to know when they're changed for --reuse-results
    def __init__(self, record_id: str, cutoff: float) -> None:
        super().__init__(record_id)
        self.cutoff = cutoff
        self.some_other_information = []  # this could be added to during analysis

    # implement a conversion to a JSON-compatible dictionary
    # all elements must one of: str, int, float, list, or a dict of those types (this can recurse)
    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """

        return {
            "schema_version": self.schema_version,
            "cutoff": self.cutoff,
            "other": [str(item) for item in self.some_other_information],  # an example only
        }

    # once _all_ analysis modules have completed, their results are added with this method
    # adding to the record during the analysis will cause issues
    def add_to_record(self, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        # any results would be added here
        # for an example of new features, see antismash.modules.tta
        # for an example of qualifiers, see antismash.modules.t2pks
        # any new feature types or qualifiers would be implemented in antismash.common.secmet,
        #   and would need to be able to be converted to and from biopython's SeqFeature without loss
        raise NotImplementedError()  # remove this when completed

    # implement a conversion from the JSON-compatible data returned by to_json()
    # this allows --results-reuse to avoid running the module again if not neccessary
    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["TemplateAnalysisResults"]:
        """ Constructs a new results instance from a JSON format and the
            original record analysed.
        """
        # check that the previous data version is the same as current, if not, discard the results
        if json["schema_version"] != TemplateAnalysisResults.schema_version:
            return None

        # as an example, checking if the example cutoff option matches that of the previous run
        options = get_config()
        if options.template_cutoff != json["cutoff"]:
            # it's nice to log some decisions to the debug logger so that it's easier to follow
            logging.debug("TemplateAnalysis cutoff has changed, discarding previous results")
            return None

        # the exact reconstruction depends on what details are stored
        # to match the conversion to JSON that would be:
        results = TemplateAnalysisResults(json["record_id"], json["cutoff"])
        for other in json["other"]:
            results.some_other_information.append(other)

        return results


def get_arguments() -> ModuleArgs:
    """ Builds any commandline argument constructs that may be required

        Returns:
            an empty or populated ModuleArgs instance
    """
    # construct the argument group, with section and prefix
    # the prefix will be enforced for all command line options for the module
    args = ModuleArgs('Additional analysis', 'template')

    # an example toggle to turn on your analysis, if not set to always be enabled
    # can also be used to turn on/off extra features of your analysis
    args.add_analysis_toggle('analysis',  # the option as it would appear on the command line without prefix, e.g. "--{prefix}-analysis"
                             dest='analysis',  # the storage location in the antismash Config object, e.g. "{prefix}_analysis"
                             default=False,             # disabled by default
                             action='store_true',       # enabled if --template-analysis is given on the commandline
                             # and finally, text to show when the user runs with --help
                             help="Compare identified clusters against a "
                                  "database of antiSMASH-predicted clusters.")

    # an example option setting a specific value
    args.add_option('cutoff',     # the option as it appears on the command line
                    dest='cutoff',  # as it appears in the antismash Config object
                    type=float,              # the type of the option (int, str, float, ...)
                    default=0.65,            # the default value of the option
                    help="Lowest GC content to annotate TTA codons at (default: %(default)s).")
    # more complicated options are possible, for further information see antismash.config.args,
    # look at how other modules construct arguments, or ask for help
    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks that the provided options are compatible with each other

        Arguments:
            options: the current antismash config object

        Returns:
            a list of strings, each string being an issue with the given options
    """
    issues = []
    # test the options used in get_arguments here, if any
    # for example, enforcing that values are within a certain range
    if not 0 < options.template_cutoff < 1:
        issues.append("Supplied cutoff is outside the range of 0 to 1: %s" % options.template_cutoff)
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
    return options.template_analysis


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[TemplateAnalysisResults]:
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
    return TemplateAnalysisResults.from_json(previous, record)


def run_on_record(record: Record, results: TemplateAnalysisResults, options: ConfigType) -> TemplateAnalysisResults:
    """ Run the analysis, unless the previous results apply to the given record

        Arguments:
            record: the Record being analysed
            results: an existing instance of the module's ModuleResults implementation (or None)
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation
    """
    # after a safety check that the results are the correct ones for the record, return them
    if isinstance(results, TemplateAnalysisResults) and results.record_id == record.id:
        return results
    # otherwise run the actual analysis and generate a results instance with your analysis results
    results = TemplateAnalysisResults(record.id, options.cutoff)
    # and return it
    return results