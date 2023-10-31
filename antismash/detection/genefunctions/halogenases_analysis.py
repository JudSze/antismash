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
    contains_feature_with_domain as contains_feature_with_single_domain)

from antismash.config import build_config

from typing import Union
from Bio.SeqFeature import FeatureLocation

from dataclasses import dataclass, field

# DEV
TRP_5_SIGNATURE = [33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186, 187, 195, 272, 302, 310, 342, 350, 400, 446, 450, 452, 482 ]
TRP_6_SIGNATURE = [19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219, 221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505]

TRP_5_SIGNATURE_RESIDUES = "VSILIREPGLPRGVPRAVLPGEA"
TRP_6_SIGNATURE_RESIDUES = "TEGCAGFDAYHDRFGNADYGLSIIAKIL"

# rename, add quality
TRP_5_LOW_QUALITY_CUTOFF = 350
TRP_5_HIGH_QUALITY_CUTOFF = 850
TRP_6_7_CUTOFF = 770

@dataclass
class Match:
    profile: str
    position: int
    confidence: float
    signature: str

@dataclass
class HalogenasesResults(ModuleResults):
    """ Example results class for the analysis module template """
    schema_version = 1  # when the data format in the results changes, this needs to be incremented

    # define whatever construction arguments are needed, record_id is required by the superclass
    # it's good to keep any command line option values here to know when they're changed for --reuse-results
    # should be default values, and every attribute should be set in init
    family: str                       # type of halogenase family (FD, SAMD, HD, VD, I/KGD)
    gbk_id: str                       # protein id in genbank file
    substrate: str          
    position: Optional[int] = None    # position number of halogenated atom
    confidence: int = 0               # not clearly identifiable enzyme, if position is not supplied, it should be weak
    signature: str = ""               # string of amino acid residues
    coenzyme: bool = False
    potential_matches: list[Match] = field(default_factory = list)

    def add_potential_matches(self, match: Match) -> None:

        self.potential_matches.append(match)

    def get_best_match(self) -> list[Match]:
        best_match = []
        
        if self.potential_matches:
            if len(self.potential_matches) == 1:
                return [self.potential_matches[0]]
            
            highest_confidence = max([profile.confidence for profile in self.potential_matches])
            for profile in self.potential_matches:
                if profile.confidence == highest_confidence:
                    best_match.append(profile)
        
        return best_match

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

        features = [self.family, self.gbk_id, self.position, self.record_id, self.substrate, self.confidence]
        for feature in features:
            record.add_feature(feature)

    # implement a conversion from the JSON-compatible data returned by to_json()
    # this allows --results-reuse to avoid running the module again if not neccessary
    @staticmethod
    def from_json(json: Dict[str, Any], record: secmet.Record) -> Optional["HalogenasesResults"]:      
        return None

# cluster_fasta will have to come from a Record object

def open_hmm_files() -> list[HmmSignature]:
    with open(path.get_full_path(__file__, "data", "halogenases", 'hmmdetails.txt'),"r", encoding = "utf-8") as handle:
        hmmdetails = [line.split("\t") for line in handle.read().splitlines() if line.count("\t") == 3]

    signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]
    return signature_profiles

def run_halogenase_phmms(cluster_fasta:  str, signature_profiles: list[HmmSignature]) -> dict[str, dict[str, Union[float, str]]]:
    
    halogenase_hmms_by_id: dict[str, dict[str, Any]] = defaultdict(dict)
    for sig in signature_profiles:
        sig.path = path.get_full_path(f'{sig.hmm_file}/{sig.name}')
        print('Im about to run hmmsearch', sig.name)
        runresults = subprocessing.run_hmmsearch(sig.path, cluster_fasta)
        for runresult in runresults:
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    if hsp.hit_id not in halogenase_hmms_by_id:
                        halogenase_hmms_by_id[hsp.hit_id] = {}
                    halogenase_hmms_by_id[hsp.hit_id][hsp.query_id] = {  # TODO convert to dataclass
                            "type": sig.name,
                            "bitscore": hsp.bitscore, 
                            "PROFILE_NAME": hsp.query_id, 
                            "PROFILE": sig.path
                        }


    return halogenase_hmms_by_id

# changed results[0].hsps index to results[1].hsps
# check the query output
def search_signatures(sequence: str, positions: list[int], halogenase_hmms_by_id: dict[str, Any], max_evalue: float = 0.1) \
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
def is_coenzyme_present(center: secmet.CDSFeature, candidates: Iterable[secmet.CDSFeature]) -> bool:
    """Should check reductase pHMM"""
    candidate_neighbors = find_neighbours_in_range(center, candidates)
    return contains_feature_with_single_domain(candidate_neighbors, {"reductase"})

def get_signatures_from_profiles(translation: str, halogenase_hmms_by_id:  dict[str, Union[str, float]]) -> dict[str, Optional[str]]:
    signatures: dict[str, Optional[str]] = {}

    for signature in [TRP_5_SIGNATURE, TRP_6_SIGNATURE]:
    # it is possible it has both signatures
    # list of tuples and ranking system with all the signatures
        residue, alignment = search_signatures(translation, signature, halogenase_hmms_by_id)
        signatures[halogenase_hmms_by_id['PROFILE_NAME']] = residue
    
    return signatures

def check_for_FDHs(cds: CDSFeature,halogenase_hmms_by_id: dict) -> Optional[HalogenasesResults]:
    if halogenase_hmms_by_id['type'] == 'Flavin-dependent':
        family, substrate = 'Flavin-dependent', 'tryptophan'
        return HalogenasesResults(family, cds.product, substrate)

    return None


def check_for_halogenases(cds: CDSFeature, halogenase_hmms_by_id: dict[str, Union[float, str]]) -> Optional[HalogenasesResults]:
    """ Input: 
            query: string of protein sequence returned by search_signature
            halogenase_hmms_by_id: dictionary of protein with hit returned by run_halogenase_phmms
            
        Output:
            enzyme: instance of Halogenase class"""
    
    if not halogenase_hmms_by_id:
        logging.debug("Hmmsearch did not return any hit.")
        return None
    
    signatures = get_signatures_from_profiles(cds.translation, halogenase_hmms_by_id)
    if signatures:
        profile_name = halogenase_hmms_by_id['PROFILE_NAME']
        residues = signatures[profile_name]
        FDHs_match = check_for_FDHs(cds, halogenase_hmms_by_id)
        if FDHs_match:
            if profile_name == 'trp_5':
                if halogenase_hmms_by_id['bitscore'] >= TRP_5_HIGH_QUALITY_CUTOFF and residues == TRP_5_SIGNATURE_RESIDUES:
                    residues = TRP_5_SIGNATURE_RESIDUES
                    FDHs_match.add_potential_matches(Match('trp_5', 5, 1, residues))
                elif halogenase_hmms_by_id['bitscore'] >= TRP_5_LOW_QUALITY_CUTOFF and residues == TRP_5_SIGNATURE_RESIDUES:
                    residues = TRP_5_SIGNATURE_RESIDUES
                    FDHs_match.add_potential_matches(Match('trp_5', 5, 0.5, residues))

            if profile_name == 'trp_6_7':       
                if residues == TRP_6_SIGNATURE_RESIDUES:
                    position, confidence = 6, 1
                else:
                    position, confidence = 7, 0.8

                FDHs_match.add_potential_matches(Match('trp_6_7', position, confidence, residues))
    
    return FDHs_match


def specific_analysis(record: Record) -> list[HalogenasesResults]:
    potential_enzymes: list[HalogenasesResults] = []

   
    features = record.get_cds_features_within_regions()
    hmmsearch_fasta = fasta.get_fasta_from_features(features)

    signature_profiles = open_hmm_files()
    hmm_hits = run_halogenase_phmms(hmmsearch_fasta, signature_profiles)
    
    for protein, hit in hmm_hits.items():
        for cds in features:
            if cds.get_name() == protein:
                for profile, hit_feature in hit.items():
                    potential_enzyme = check_for_halogenases(cds, hit_feature)
                    if potential_enzyme:
                        potential_enzymes.append(potential_enzyme)

    for enzyme in potential_enzymes:
        best_matches = enzyme.get_best_match()
        assert isinstance(best_matches, list), best_matches
        if not best_matches:
            continue
        assert len(best_matches) == 1, best_matches
        best_match = best_matches[0]
        enzyme.position = best_match.position
        enzyme.confidence = best_match.confidence
        enzyme.signature = best_match.signature


    return potential_enzymes