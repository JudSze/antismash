from typing import Any, Dict, Optional, Iterable, Union
import logging
from dataclasses import dataclass, field

from antismash.common.secmet import Record, CDSFeature
from antismash.common.module_results import ModuleResults
from antismash.common.module_results import ModuleResults
from antismash.common.signature import HmmSignature
from antismash.modules.lanthipeptides.specific_analysis import (
    find_neighbours_in_range,
    contains_feature_with_domain)
from antismash.common import (
    secmet,
    subprocessing,
    path,
    fasta,
    utils,
)

# RESIDUE SIGNATURES
TRP_5_SIGNATURE = [33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186, 187, 195, 272, 302, 310, 342, 350, 400, 446, 450, 452, 482 ]
TRP_6_SIGNATURE = [19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219, 221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505]

TRP_5_SIGNATURE_RESIDUES = "VSILIREPGLPRGVPRAVLPGEA"
TRP_6_SIGNATURE_RESIDUES = "TEGCAGFDAYHDRFGNADYGLSIIAKIL"

# CUTOFFS FOR pHMMs
TRP_5_LOW_QUALITY_CUTOFF = 350
TRP_5_HIGH_QUALITY_CUTOFF = 850
TRP_6_7_CUTOFF = 770

@dataclass
class HmmResult:
    hit_id: str
    query_id: str
    type: str
    bitscore: float
    profile_name: str
    profile: str


@dataclass
class Match:
    profile: str
    position: int
    confidence: Union[int, float]
    signature: str

@dataclass
class HalogenasesResults:
    """ Result class of the categorized enzymes """
    schema_version = 1  # when the data format in the results changes, this needs to be incremented

    gbk_id: str                                                         # protein id in genbank file
    family: str = "Halogenase"                                          # type of halogenase family (FD, SAMD, HD, VD, I/KGD)  
    substrate: str = ""                                                 # name of the common substrate (e.g. tryptophan, pyrrole, aliphatic molecule)
    position: Optional[int] = None                                      # position number of halogenated atom
    confidence: int = 0                                                 # confidence of the categorization
    signature: str = ""                                                 # string of amino acid residues
    coenzyme: bool = False                                              # is coenzyme present or not
    potential_matches: list[Match] = field(default_factory = list)      # possible categorizations (e.g. if enzyme meets requirements for more groups)

    def add_potential_matches(self, match: Match) -> None:
        """ Adds the features of a group to which the enzyme belongs to a list"""
        self.potential_matches.append(match)

    def get_best_match(self) -> list[Match]:
        """ If enzyme meets the requirements for several groups, 
            it compares the confidences of the categorizations and returns the one with the highest confidence."""
        best_match = []
        
        if self.potential_matches:
            if len(self.potential_matches) == 1:
                return [self.potential_matches[0]]
            
            highest_confidence = max([profile.confidence for profile in self.potential_matches])
            for profile in self.potential_matches:
                if profile.confidence == highest_confidence:
                    best_match.append(profile)
        
        return best_match
    
    def specify_details(self, best_match: Match) -> None:
        self.position = best_match.position
        self.confidence = best_match.confidence
        self.signature = best_match.signature

    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """

        return {
            "family": self.family,
            "gbk_id": self.gbk_id,
            "substrate": self.substrate,
            "position": self.position,
            "confidence": self.confidence,
            "signature": self.signature,
            "coenzyme": self.coenzyme,
            "potential_matches": self.potential_matches
        }

    # depends on where the module is in the antiSMASH timeline
    # def add_to_record(self, cds: CDSFeature, record: Record) -> None:
    #     """ Adds the analysis results to the record """
    #     if record.id != self.record_id:
    #         raise ValueError("Record to store in and record analysed don't match")

    #     features = [self.family, self.gbk_id, self.position, self.record_id, self.substrate, self.confidence]
    #     for feature in features:
    #         record.add_feature(feature)

    # implement a conversion from the JSON-compatible data returned by to_json()
    # this allows --results-reuse to avoid running the module again if not neccessary
    @staticmethod
    def from_json(json: Dict[str, Any], record: secmet.Record) -> Optional["HalogenasesResults"]:
        """ Constructs the HalogenasesResults from the JSON representation """      
        return None

def get_profile_signature() -> list[HmmSignature]:
    """ Get the features of a certain profile from a txt file """
    with open(path.get_full_path(__file__, "data", "halogenases", 'hmmdetails.txt'),"r", encoding = "utf-8") as handle:
        hmmdetails = [line.split("\t") for line in handle.read().splitlines() if line.count("\t") == 3]

    signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]
    return signature_profiles

def run_halogenase_phmms(cluster_fasta: str, signature_profiles: list[HmmSignature]) -> dict[str, list[HmmResult]]:
    """ Check if protein sequences hit any pHMM
        Input: 
            cluster_fasta: string of protein sequences in a fasta format
            signature_profiles: features of a profile returned by get_profile_signature
        Output:
            instance of the HmmResult class
            If there are no hits it returns an empty HmmResult instance
    """
    halogenase_hmms_by_id: dict = {}
    for sig in signature_profiles:
        sig.path = path.get_full_path(f'{sig.hmm_file}/{sig.name}')

        runresults = subprocessing.run_hmmsearch(sig.path, cluster_fasta)
        for runresult in runresults:
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    hit = HmmResult(hsp.hit_id, hsp.query_id, sig.name, hsp.bitscore, hsp.query_id, sig.path)
                    if hsp.hit_id not in halogenase_hmms_by_id.keys():
                        halogenase_hmms_by_id[hsp.hit_id] = []
                    else:
                        halogenase_hmms_by_id[hsp.hit_id].append(hit)

    return halogenase_hmms_by_id

# changed results[0].hsps index to results[1].hsps

def search_signature_residues(sequence: str, positions: list[int], hmm_result: HmmResult, max_evalue: float = 0.1) \
        -> Optional[str]:
    """ Get the signature residues from the pHMM for the searched protein sequence
        Input:
            sequence: protein sequence as str in fasta format
            positions: list of position numbers in the pHMM
            hmm_result: instance of the HmmResult class
            max_evalue: set to 0.1
        Output:
            residues that are present in the given positions
    """
    args = ["-E", str(max_evalue)]
    results = subprocessing.hmmpfam.run_hmmpfam2(hmm_result.profile, f">query\n{sequence}", extra_args=args)

    if not (results and results[0].hsps):
        logging.debug("no hits for query %s, is this a legitimate enzyme?")
        return None

    found = False
    hit = None  
    for hit in results[0].hsps:
        if hit.hit_id == hmm_result.profile_name:
            found = True
            break

    if not found:
        logging.debug(
            "no hits for the enzyme in %s, is this a legitimate sequence?", hmm_result.profile_name)
        return None

    profile = hit.aln[1].seq
    query = hit.aln[0].seq
    offset = hit.hit_start

    sites = utils.extract_by_reference_positions(query, profile, [p - offset for p in positions if offset < p])

    return sites


def is_coenzyme_present(center: secmet.CDSFeature, candidates: Iterable[secmet.CDSFeature]) -> bool:
    """ Check if coenzyme is present in the cluster """
    candidate_neighbors = find_neighbours_in_range(center, candidates)
    return contains_feature_with_domain(candidate_neighbors, {"reductase"})

def get_residues_from_profiles(translation: str, hmm_result: HmmResult) -> dict[str, Optional[str]]:
    """ Get signature residues for an enzyme from each pHMM
        Input:
            translation: string of protein sequence formatted in fasta
            hmm_result: instance of HmmResult class
        Output:
            dictionary of the profile name as key and the signature residues as value
    """
    signature_residues: dict[str, Optional[str]] = {}

    for signature in [TRP_5_SIGNATURE, TRP_6_SIGNATURE]:

        residue = search_signature_residues(translation, signature, hmm_result)
        signature_residues[hmm_result.profile_name] = residue
    
    return signature_residues

def check_for_halogenases(cds: CDSFeature, halogenase_hmms_by_id: list[HmmResult]) -> Optional[HalogenasesResults]:
    """ Categorizes enzymes based on wether they hit the pHMM and have the required signature residues
        Input: 
            cds: CDSFeature object
            hmm_result: HmmResult object
        Output:
            HalogenasesResults object
    """
    if not halogenase_hmms_by_id:
        logging.debug("Hmmsearch did not return any hit.")
        return None
    
    halogenase_match = HalogenasesResults(cds.get_name())
    for hit in halogenase_hmms_by_id:
        signature_residues = get_residues_from_profiles(cds.translation, hit)
        if signature_residues:
            profile_name = hit.profile_name
            residues = signature_residues[profile_name]
            if halogenase_match:
                halogenase_match.family = "Flavin-dependent halogenase"
                if profile_name == 'trp_5':
                    if hit.bitscore >= TRP_5_HIGH_QUALITY_CUTOFF and residues == TRP_5_SIGNATURE_RESIDUES:
                        halogenase_match.add_potential_matches(Match('trp_5', 5, 1, residues))
                    elif hit.bitscore >= TRP_5_LOW_QUALITY_CUTOFF and residues == TRP_5_SIGNATURE_RESIDUES:
                        halogenase_match.add_potential_matches(Match('trp_5', 5, 0.5, residues))

                if profile_name == 'trp_6_7':       
                    if residues == TRP_6_SIGNATURE_RESIDUES:
                        halogenase_match.add_potential_matches(Match('trp_6_7', 6, 1, residues))
                    else:
                        halogenase_match.add_potential_matches(Match('trp_6_7', 7, 0.8, residues))
    
    return halogenase_match


def specific_analysis(record: Record) -> list[HalogenasesResults]:
    """ Final categorization of enzyme, where the best match is taken to define the
        substrate, position and signature
        Input: record object
        Output: list of HalogenasesResults objects
    """
    potential_enzymes: list[HalogenasesResults] = []
    final_enzymes: list[HalogenasesResults] = []

   
    features = record.get_cds_features_within_regions()
    hmmsearch_fasta = fasta.get_fasta_from_features(features)

    signature_profiles = get_profile_signature()
    hmm_hits = run_halogenase_phmms(hmmsearch_fasta, signature_profiles)

    for protein, hits in hmm_hits.items():
        cds = record.get_cds_by_name(protein)
        potential_enzyme = check_for_halogenases(cds, hits)
        if potential_enzyme:
            potential_enzymes.append(potential_enzyme)

    for enzyme in potential_enzymes:
        best_matches = enzyme.get_best_match()
        assert isinstance(best_matches, list), best_matches
        if not best_matches:
            continue
        assert len(best_matches) == 1, best_matches
        best_match = best_matches[0]
        enzyme.specify_details(best_match)

    return potential_enzymes