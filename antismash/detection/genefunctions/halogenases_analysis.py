from collections import defaultdict
from typing import Any, Dict, Optional, Iterable

import logging

from antismash.common.secmet import Record, CDSFeature, Region
# from antismash.config import ConfigType
# from antismash.detection.nrps_pks_domains import ModularDomain

from antismash.common.secmet.locations import location_from_string
from antismash.common.module_results import ModuleResults
from antismash.common import (
    all_orfs,
    comparippson,
    module_results,
    secmet,
    subprocessing,
    path,
    fasta,
    utils,
)

from antismash.common.module_results import ModuleResults

from antismash.common.signature import HmmSignature
from antismash.modules.lanthipeptides.specific_analysis import (
    find_neighbours_in_range,
    contains_feature_with_single_domain)

from antismash.config import build_config

from typing import Union
from Bio.SeqFeature import FeatureLocation

# from antismash.common.fasta import get_fasta_from_features, get_fasta_from_record

# TEST
from Bio import SeqIO

# test = secmet.Record.from_genbank("/home/szenei/antismash/MAGic-MOLFUN/kutzneride_copy.gbk")
# for record in records:

build_config([])

test = secmet.Record.from_genbank("/home/szenei/antismash/MAGic-MOLFUN/kutzneride.gbk")

record = test[0]
test_translations = []
test_cds_features = []

features = record.get_all_features()
for gene in features:
    location = FeatureLocation(gene.location.start, gene.location.end, 1)
    secmet.features.cds_feature._verify_location(location)
    protein = fasta.get_fasta_from_features(record.get_cds_features_within_location(location))
    test_translations.append(protein)

    protein_features = protein.split("\n", 1)

    if ">" not in protein_features[1]:
        plus_cds = CDSFeature(location, protein_features[1], record, "EU074211", protein_features[0], "")
        test_cds_features.append(plus_cds)

    # hmmsearch
    # signATURES
    # return Halogenase

# DEV
TRP_5_SIGNATURE = [33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186, 187, 195, 272, 302, 310, 342, 350, 400, 446, 450, 452, 482 ]
TRP_6_SIGNATURE = [19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219, 221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505]

TRP_5_SIGNATURE_RESIDUES = "VSILIREPGLPRGVPRAVLPGEA"
TRP_6_SIGNATURE_RESIDUES = "TEGCAGFDAYHDRFGNADYGLSIIAKIL"

# rename, add quality
TRP_5_LOW_CUTOFF = 350
TRP_5_HIGH_CUTOFF = 850
TRP_6_7_CUTOFF = 770

class HalogenasesResults(ModuleResults):
    """ Example results class for the analysis module template """
    schema_version = 1  # when the data format in the results changes, this needs to be incremented

    # define whatever construction arguments are needed, record_id is required by the superclass
    # it's good to keep any command line option values here to know when they're changed for --reuse-results
    # should be default values, and every attribute should be set in init
    def __init__(self, record_id: str, gbk_id: str, family: str,  position: int, substrate: str, uncertainity: str = "", signature: str = ""):
        super().__init__(record_id)
        # change this
        self.family = family                   # type of halogenase family (FD, SAMD, HD, VD, I/KGD)
        self.gbk_id = gbk_id                            # protein id in genbank file
        self.substrate = substrate          
        self.position = position       # position number of halogenated atom
        self.uncertainity = uncertainity     # not clearly identifiable enzyme, if position is not supplied, it should be weak
        self.signature = signature                        # string of amino acid residues

    # __str__
    def get_properties(self):
        """Short description of enzyme"""
        return f"Family: {self.family} Substrate: {self.substrate} Position: {self.position} Signature: {self.signature}"

    # implement a conversion to a JSON-compatible dictionary
    # all elements must one of: str, int, float, list, or a dict of those types (this can recurse)
    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """

        return {
            "record_id": self.record_id,
            "schema_version": self.schema_version,
            "halogenase_status": ""
        }

    # depends on where the module is in the antiSMASH timeline
    def add_to_record(self, cds: CDSFeature, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")

        features = [self.family, self.gbk_id, self.position, self.record_id, self.substrate, self.uncertainity]
        for feature in features:
            record.add_feature(feature)

    # implement a conversion from the JSON-compatible data returned by to_json()
    # this allows --results-reuse to avoid running the module again if not neccessary
    @staticmethod
    def from_json(json: Dict[str, Any], record: secmet.Record) -> Optional["HalogenasesResults"]:      
        return None

# cluster_fasta will have to come from a Record object
# use dataclass instead of the dictionary
def run_halogenase_phmms(cluster_fasta:  str, txt: str):
    
    with open(path.get_full_path(__file__, "data", txt),"r", encoding = "utf-8") as handle:
        hmmdetails = [line.split("\t") for line in handle.read().splitlines() if line.count("\t") == 3]

    signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]

    halogenase_hmms_by_id: Dict[str, Any] = {}
    for sig in signature_profiles:
        sig.path = path.get_full_path(f'{sig.hmm_file}/{sig.name}')
        runresults = subprocessing.run_hmmsearch(sig.path, cluster_fasta)
        for runresult in runresults:
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    halogenase_hmms_by_id[hsp.query_id] = {"type": sig.name,
                                                         "bitscore": hsp.bitscore, 
                                                         "PROFILE_NAME": hsp.query_id, 
                                                         "PROFILE": sig.path, # maybe unpractitional
                                                         "protein": hsp.hit_id}

    return halogenase_hmms_by_id

# changed results[0].hsps index to results[1].hsps
# check the query output
def search_signatures(sequence, positions, halogenase_hmms_by_id, max_evalue: float = 0.1) \
        -> tuple[Optional[str], Optional[str]]:
    # get the signature residues from the pHMM for the search protein sequence
    args = ["-E", str(max_evalue)]
    results = subprocessing.hmmpfam.run_hmmpfam2(halogenase_hmms_by_id["PROFILE"], f">query\n{sequence}", extra_args=args)
    # if not (results and results[1].hsps):
    if not (results and results[0].hsps):
        logging.debug("no hits for query %s, is this a legitimate enzyme?")
        return None, None

    found = False
    hit = None  # will be set in the loop or we abort anyway, just to make pylint happy
    for hit in results[0].hsps:
        if hit.hit_id == halogenase_hmms_by_id["PROFILE_NAME"]:
            found = True
            break

    if not found:
        logging.debug(
            "no hits for the enzyme in %s, is this a legitimate sequence?", halogenase_hmms_by_id["PROFILE_NAME"])
        return None, None

    profile = hit.aln[1].seq
    query = hit.aln[0].seq
    offset = hit.hit_start

    sites = utils.extract_by_reference_positions(query, profile, [p - offset for p in positions if offset < p])

    return sites, query

# search for it in the nearby, otherwise let it go
def is_coenzyme_present(center: secmet.CDSFeature, candidates: Iterable[secmet.CDSFeature]):
    """Should check reductase pHMM"""
    candidate_neighbors = find_neighbours_in_range(center, candidates)
    if contains_feature_with_single_domain(candidate_neighbors, {"reductase"}):
        return 'reductase'

def classify_halogenase(cds: CDSFeature, halogenase_hmms_by_id: dict, coenzyme_present: Union[str, bool]) -> HalogenasesResults:
    """ Input: 
            query: string of protein sequence returned by search_signature
            halogenase_hmms_by_id: dictionary of protein with hit returned by run_halogenase_phmms
            
        Output:
            enzyme: instance of Halogenase class"""
    
    if not halogenase_hmms_by_id:
        logging.debug("Hmmsearch did not return any hit.")
        return None

    for signature1 in [TRP_5_SIGNATURE, TRP_6_SIGNATURE]:
        # it is possible it has both signatures
        # list of tuples and ranking system with all the signatures
        sig_search_res = search_signatures(cds.translation, signature1, halogenase_hmms_by_id)

        if sig_search_res[0] == TRP_5_SIGNATURE_RESIDUES:
            signature = sig_search_res[0]
        elif sig_search_res[0] == TRP_6_SIGNATURE_RESIDUES:
            signature = sig_search_res[0]
        else:
            signature = ""

    # fragment it into more functions
    # ranking systems for 5-6-7, make decision based on the ranking
    if halogenase_hmms_by_id['PROFILE_NAME'] == 'trp_5':
        substrate = "tryptophan"
        if halogenase_hmms_by_id['bitscore'] >= TRP_5_HIGH_CUTOFF and signature == TRP_5_SIGNATURE_RESIDUES:
            position = 5
        elif halogenase_hmms_by_id['bitscore'] >= TRP_5_LOW_CUTOFF and bool(coenzyme_present) and signature == TRP_5_SIGNATURE_RESIDUES:
            uncertainity = "might be a funky Trp-5 halogenase" # shouldn't be a text
        else:
            logging.info("Enzyme did not much any of the Trp-5 pHMM cutoff.")
            return None

    if halogenase_hmms_by_id['PROFILE_NAME'] == 'trp_6_7':
        substrate = "tryptophan"
        if signature == TRP_6_SIGNATURE_RESIDUES:
            position = 6
        else:
            position = 7

    enzyme = HalogenasesResults(record.id, cds.get_name(), halogenase_hmms_by_id['type'], position, substrate, signature=signature) # shouldn't be created if its empty
    return enzyme

def specific_analysis(cds: CDSFeature) -> HalogenasesResults:
    hmm_search_res = dict()
    potential_enzymes = list()

    hmm_hit = run_halogenase_phmms(f">{cds.product}\n{cds.translation}", "hmmdetails.txt")

    for index, res in hmm_hit.items():
        if cds.product not in hmm_search_res.keys():
            hmm_search_res[cds.product] = [res]
        else:
            hmm_search_res[cds.product].append(res)
    
    for protein, hit in hmm_hit.items():
        potential_enzymes.append(classify_halogenase(cds, hit, True))

    return potential_enzymes


#def sig(sequence, hmm_results):
    # handle gaps in hits
    # return sig
                
TRP_5 = "/home/szenei/antismash/MAGic-MOLFUN/antismash/modules/halogenases/data/trp_5_v2.hmm"
TRP_6_7 = "/home/szenei/antismash/MAGic-MOLFUN/antismash/modules/halogenases/data/trp_6_7_v2.hmm"

if __name__ == "__main__":
    try:
        for feature in test_cds_features[0:4]:
            x = specific_analysis(feature)
            for enzyme in x:
                if enzyme != None:
                    print(enzyme.get_properties())
    except RuntimeError:
        print("File or directory does not exist or is misformatted")
