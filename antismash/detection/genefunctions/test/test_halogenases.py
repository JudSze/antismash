from typing import Dict

import random
import unittest
from unittest.mock import patch

from antismash.common import utils
from antismash.common import secmet
from antismash.common import subprocessing
from antismash.common.fasta import read_fasta
from antismash.common.signature import HmmSignature
# fake classes for testing
from antismash.common.secmet.test.helpers import (
    DummyRecord,
    DummyFeature,
    DummyCDS)
from antismash.common.test.helpers import FakeHSPHit, FakeHit
#functions and classes to test
from antismash.detection.genefunctions import halogenases_analysis
from antismash.detection.genefunctions.halogenases_analysis import (
    HalogenasesResults,
    HmmResult,
    Match,
    get_profile_signature,
    run_halogenase_phmms,
    search_signature_residues,
    is_coenzyme_present,
    get_fdh_residues,
    check_for_fdh,
    check_for_halogenases,
    specific_analysis
)

TRP_5_SIGNATURE_RESIDUES = "VSILIREPGLPRGVPRAVLPGEA"
TRP_6_SIGNATURE_RESIDUES = "TEGCAGFDAYHDRFGNADYGLSIIAKIL"

class TestHalogenasesAnalysis(unittest.TestCase):
    def setUp(self):
        self.test_match = Match("trp_6_7", 6, 1, "")
        self.positive_test_hmm_result = HmmResult(hit_id='foo', query_id='trp_6_7', type='Flavin-dependent',
                            bitscore=1000, profile_name='trp_6_7',
                            profile='/home/szenei/antismash_tailoring_enzymes/antismash/antismash/detection/genefunctions/data/halogenases/trp_6_7_v2.hmm')
        self.positive_test_halogenase = HalogenasesResults("ktzR", potential_matches=[self.test_match])

        self.negtive_test_halogenase = HalogenasesResults("ktzR")

    def mock_get_cds_name(self):
        return ""

    def test_get_best_match(self):
        negative_test_best_match = self.negtive_test_halogenase.get_best_match()
        assert bool(negative_test_best_match) is False

        positive_test_best_match = self.positive_test_halogenase.get_best_match()
        assert len(positive_test_best_match) > 0 or isinstance(positive_test_best_match, Match)

        self.positive_test_halogenase.add_potential_matches(self.test_match)
        assert len(self.positive_test_halogenase.potential_matches) >= 2

        multiple_matches = self.positive_test_halogenase.get_best_match()
        assert len(multiple_matches) >= 2 or isinstance(positive_test_best_match, Match)

    def test_details_and_conversion_methods(self):
        self.positive_test_halogenase.specify_details(self.test_match)
        assert (bool(self.positive_test_halogenase.position)
                or bool(self.positive_test_halogenase.confidence)
                or bool(self.positive_test_halogenase.signature))
        
        converted_to_json = self.positive_test_halogenase.to_json()
        assert isinstance(converted_to_json, Dict)

        converted_from_json = self.positive_test_halogenase.from_json(converted_to_json)
        assert isinstance(converted_from_json, HalogenasesResults)


    def test_get_profile_signature(self):
        profile_signatures = get_profile_signature()
        for result in profile_signatures:
            assert isinstance(result, HmmSignature)

        with patch.object(halogenases_analysis, "HmmSignature", return_value=None):
            negative_profile_signatures = get_profile_signature()
            assert negative_profile_signatures is None

    @patch.object(subprocessing, "run_hmmsearch", return_value=[FakeHit("start", "end", 1000, "foo")])
    def test_run_halogenase_phmms(self, run_hmmsearch):
            for value in run_hmmsearch.return_value:
                value.hsps = [FakeHSPHit("foo", "foo", bitscore=250)]

            hmm_signature = get_profile_signature()
            negative_test_halogenase_hmms_by_id = run_halogenase_phmms("", list(hmm_signature))
            assert not negative_test_halogenase_hmms_by_id
            
            for value in run_hmmsearch.return_value:
                value.hsps = [FakeHSPHit("foo", "foo", bitscore=1000)]

            positive_test_halogenase_hmms_by_id = run_halogenase_phmms("", list(hmm_signature))
            for hit in positive_test_halogenase_hmms_by_id["foo"]:
                assert isinstance(hit, HmmResult)

    @patch.object(utils, "extract_by_reference_positions", return_value = "")
    def test_search_signature_residues(self, reference_positions):
        test_fasta = read_fasta("/home/szenei/antismash_tailoring_enzymes/test.fasta")
        positions = [random.randrange(512, 600) for number in range(10)]

        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2", return_value = []) as run_hmmpfam2:
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("foo", hit_id="foo")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = "ghjvghjkbln"
                    hit_profile.seq = "xfdhcgkbjlnkæml"
                    hit.aln = [hit_profile, hit_query]

            signature_residues = search_signature_residues(list(test_fasta.values())[0], positions, self.positive_test_hmm_result)
            assert signature_residues == None


        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2", return_value = [FakeHit("start", "end", 1000, "foo")]) as run_hmmpfam2:
            # checking for hit_id works
            for value in run_hmmpfam2.return_value:
                value.hsps = [FakeHSPHit("foo", hit_id=self.positive_test_hmm_result.profile_name)]
                for hit in value.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = "ghjvghjkbln"
                    hit_profile.seq = "xfdhcgkbjlnkæml"
                    hit.aln = [hit_profile, hit_query]

            signature_residues = search_signature_residues(list(test_fasta.values())[0], positions, self.positive_test_hmm_result)
            assert signature_residues == ""

            # checking if hit_id is bad it breaks
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("foo", hit_id="foo")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = "ghjvghjkbln"
                    hit_profile.seq = "xfdhcgkbjlnkæml"
                    hit.aln = [hit_profile, hit_query]

            signature_residues = search_signature_residues(list(test_fasta.values())[0], positions, self.positive_test_hmm_result)
            assert signature_residues == None


    def test_is_coenzyme_present(self):
        with self.assertRaises(NotImplementedError):
            is_coenzyme_present(DummyFeature(), [DummyFeature()])

    def test_get_fdh_residues(self):
        test_fasta = read_fasta("/home/szenei/antismash_tailoring_enzymes/test.fasta")
        fdh_residues = get_fdh_residues(list(test_fasta.values())[0],  self.positive_test_hmm_result)
        assert isinstance(fdh_residues, dict)

    def test_check_for_fdh(self):
        trp_7_fdh = check_for_fdh(DummyCDS(), self.positive_test_halogenase,
                                        "trp_6_7",
                                        self.positive_test_hmm_result)
        assert isinstance(trp_7_fdh, HalogenasesResults)

        with patch.object(halogenases_analysis, "get_fdh_residues", return_value = {"trp_6_7": TRP_6_SIGNATURE_RESIDUES}):
            trp_6 = check_for_fdh(DummyCDS(), self.positive_test_halogenase,
                                            "trp_6_7",
                                            self.positive_test_hmm_result)
            assert isinstance(trp_6, HalogenasesResults)

        with patch.object(halogenases_analysis, "get_fdh_residues", return_value = {"trp_5": TRP_5_SIGNATURE_RESIDUES}):
            high_quality_trp_5 = check_for_fdh(DummyCDS(), self.positive_test_halogenase,
                                            "trp_5",
                                            self.positive_test_hmm_result)
            assert isinstance(high_quality_trp_5, HalogenasesResults)
        
        with patch.object(halogenases_analysis, "get_fdh_residues", return_value = {"trp_5": TRP_5_SIGNATURE_RESIDUES}):
            low_quality_hit = HmmResult("trp_5", "trp_5_v2", "foo", 400, "trp_5_v2", "foo")
            low_quality_trp_5 = check_for_fdh(DummyCDS(), self.positive_test_halogenase, "trp_5", low_quality_hit)
            assert isinstance(low_quality_trp_5, HalogenasesResults)

    def test_check_for_halogenases(self):
        negative_checked_halogenases = check_for_halogenases(DummyFeature(), [])
        assert negative_checked_halogenases == None

        DummyFeature.get_name = TestHalogenasesAnalysis.mock_get_cds_name
        mock_cds_feature = DummyFeature()
        mock_cds_feature.translation = ""
        positive_checked_halogenases = check_for_halogenases(mock_cds_feature, [self.positive_test_hmm_result])
        assert isinstance(positive_checked_halogenases, HalogenasesResults)

    @patch.object(secmet.Record, "get_cds_by_name", return_value = "KtzR")
    @patch.object(subprocessing, "run_hmmsearch", return_value=[FakeHit("start", "end", 1000, "foo")])
    @patch.object(halogenases_analysis, "run_halogenase_phmms")
    @patch.object(halogenases_analysis, "check_for_halogenases")
    def test_specific_analysis(self, check_for_halogenase, run_halogenase_phmms, run_hmmsearch, get_cds_by_name):
        run_halogenase_phmms.return_value = {"trp_5": [self.positive_test_hmm_result]}
        for value in run_hmmsearch.return_value:
            value.hsps = [FakeHSPHit("foo", "foo", bitscore=250)]
        
        check_for_halogenase.return_value = self.positive_test_halogenase
        record = DummyRecord("", seq='>x\ngggipw')
        positive_test = specific_analysis(record)
        assert positive_test

        with patch.object(halogenases_analysis.HalogenasesResults, "get_best_match", return_value = []):
            check_for_halogenase.return_value = HalogenasesResults("id", "family")
            record = DummyRecord("", seq='>x\ngggipw')
            positive_test = specific_analysis(record)
            assert check_for_halogenase.return_value.position is None
