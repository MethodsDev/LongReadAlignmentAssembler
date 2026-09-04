#!/usr/bin/env python3

from io import StringIO

import LRAA_Globals
from Quantify import Quantify
from gene_symbol_utils import combine_gene_name_and_id


class _FakeMultipath:
    def __init__(self, mp_id, read_names):
        self._mp_id = mp_id
        self._read_names = set(read_names)

    def get_id(self):
        return self._mp_id

    def get_simple_path(self):
        """Node path, as the real MultiPath exposes.

        report_quant_results reads this to decide whether the read's intron chain is
        exactly the isoform's (is_FSM / uniq_FSM_reads). Returns a chain that matches
        no transcript double's path, so every row here is non-FSM -- which is correct
        for a fixture that asserts TPM renormalization and says nothing about FSM.
        """
        return ["E:mp5p", "I:mp", "E:mp3p"]

    def get_read_names(self):
        return self._read_names


class _FakeTranscript:
    def __init__(
        self,
        gene_id,
        transcript_id,
        read_count,
        multipath,
        multipath_weight=1.0,
        num_exons=1,
        gene_name=None,
    ):
        self._gene_id = gene_id
        self._transcript_id = transcript_id
        self._read_count = read_count
        self._multipath = multipath
        self._multipath_weight = multipath_weight
        self._num_exons = num_exons
        self._gene_name = gene_name

    def get_gene_id(self):
        return self._gene_id

    def get_output_gene_id(self):
        return combine_gene_name_and_id(self._gene_name, self._gene_id)

    def get_transcript_id(self):
        return self._transcript_id

    def get_read_counts_assigned(self):
        return self._read_count

    def get_isoform_fraction(self):
        return 1.0

    def get_multipaths_evidence_assigned(self):
        return [self._multipath]

    def get_multipath_weight(self, multipath):
        assert multipath == self._multipath
        return self._multipath_weight

    def get_TPM(self):
        return self._read_count / LRAA_Globals.config["num_total_reads"] * 1e6

    def get_num_exon_segments(self):
        return self._num_exons

    def get_exons_string(self):
        return "1-100"

    def get_simple_path(self):
        """Node path, as the real Transcript exposes.

        report_quant_results compares a read's intron chain to the isoform's to set the
        tracking file's is_FSM and the quant.expr uniq_FSM_reads column. This test
        asserts TPM renormalization, not FSM, so the path only has to be well formed
        and distinct per transcript.
        """
        return ["E:5p:" + str(self.get_transcript_id()),
                "I:" + str(self.get_transcript_id()),
                "E:3p:" + str(self.get_transcript_id())]

    def get_introns_string(self):
        return "151-199" if self._num_exons > 1 else ""


def test_report_quant_results_tpm_renormalizes_over_reported_transcripts():
    old_num_total_reads = LRAA_Globals.config["num_total_reads"]
    old_weight_reads = LRAA_Globals.config["weight_reads_by_3prime_agreement"]
    old_debug = LRAA_Globals.DEBUG
    try:
        LRAA_Globals.config["num_total_reads"] = 1000
        LRAA_Globals.config["weight_reads_by_3prime_agreement"] = True
        LRAA_Globals.DEBUG = False

        mp1 = _FakeMultipath("mp1", ["r1", "r2"])
        mp2 = _FakeMultipath("mp2", ["r3"])
        transcripts = [
            _FakeTranscript("gene1", "tx1", 100.0, mp1, 0.250, gene_name="GENE1"),
            _FakeTranscript("gene2", "tx2", 300.0, mp2, 1.000, num_exons=2),
        ]
        frac_assignments = {
            "tx1": {mp1: 1.0},
            "tx2": {mp2: 1.0},
        }

        quant_out = StringIO()
        tracking_out = StringIO()
        Quantify(False, 1).report_quant_results(
            transcripts, frac_assignments, quant_out, tracking_out
        )

        rows = [line.split("\t") for line in quant_out.getvalue().strip().splitlines()]
        tpm_sum = sum(float(row[6]) for row in rows)
        quant_by_transcript = {row[1]: row for row in rows}
        assert quant_by_transcript["tx1"][0] == "GENE1^gene1"
        assert quant_by_transcript["tx2"][0] == "gene2"
        # RPM_total_reads is the SECOND-to-last column: uniq_FSM_reads is appended
        # after it, deliberately last so that every column before it keeps its index.
        # row[-1] would read the FSM count here.
        rpm_total_reads_sum = sum(float(row[-2]) for row in rows)
        uniq_FSM_sum = sum(int(row[-1]) for row in rows)

        assert round(tpm_sum, 3) == 1000000.0
        assert round(rpm_total_reads_sum, 3) == 400000.0
        # The multipath double's chain matches no transcript double's, so no row is a
        # full splice match. Asserted so the column's presence is covered here too.
        assert uniq_FSM_sum == 0

        tracking_rows = [
            line.split("\t") for line in tracking_out.getvalue().strip().splitlines()
        ]
        tracking_by_transcript = {row[1]: row for row in tracking_rows}
        assert tracking_by_transcript["tx1"][0] == "GENE1^gene1"
        assert tracking_by_transcript["tx2"][0] == "gene2"
        # 8 original columns plus is_unique and is_FSM.
        assert all(len(row) == 10 for row in tracking_rows)
        # Both are 0/1 flags. NOT asserted: that is_FSM implies is_unique -- it does
        # not. is_FSM is a property of the (read, isoform) pair, so a read whose chain
        # matches several endpoint variants of that chain is FSM for each of them while
        # being unique to none. Only the AGGREGATE applies the conjunction, which is
        # why uniq_FSM_reads counts rows with both flags rather than the is_FSM rows.
        for row in tracking_rows:
            assert row[8] in ("0", "1"), "is_unique must be a 0/1 flag"
            assert row[9] in ("0", "1"), "is_FSM must be a 0/1 flag"
        assert {row[3] for row in tracking_rows if row[1] == "tx1"} == {"1"}
        assert {row[3] for row in tracking_rows if row[1] == "tx2"} == {"2"}
        assert {row[7] for row in tracking_rows if row[1] == "tx1"} == {"0.250000"}
    finally:
        LRAA_Globals.config["num_total_reads"] = old_num_total_reads
        LRAA_Globals.config["weight_reads_by_3prime_agreement"] = old_weight_reads
        LRAA_Globals.DEBUG = old_debug
