#!/usr/bin/env python3

from io import StringIO

import LRAA_Globals
from Quantify import Quantify


class _FakeMultipath:
    def __init__(self, mp_id, read_names):
        self._mp_id = mp_id
        self._read_names = set(read_names)

    def get_id(self):
        return self._mp_id

    def get_read_names(self):
        return self._read_names


class _FakeTranscript:
    def __init__(self, gene_id, transcript_id, read_count, multipath):
        self._gene_id = gene_id
        self._transcript_id = transcript_id
        self._read_count = read_count
        self._multipath = multipath

    def get_gene_id(self):
        return self._gene_id

    def get_transcript_id(self):
        return self._transcript_id

    def get_read_counts_assigned(self):
        return self._read_count

    def get_isoform_fraction(self):
        return 1.0

    def get_multipaths_evidence_assigned(self):
        return [self._multipath]

    def get_TPM(self):
        return self._read_count / LRAA_Globals.config["num_total_reads"] * 1e6

    def get_num_exon_segments(self):
        return 1

    def get_exons_string(self):
        return "1-100"

    def get_introns_string(self):
        return ""


def test_report_quant_results_tpm_renormalizes_over_reported_transcripts():
    old_num_total_reads = LRAA_Globals.config["num_total_reads"]
    old_debug = LRAA_Globals.DEBUG
    try:
        LRAA_Globals.config["num_total_reads"] = 1000
        LRAA_Globals.DEBUG = False

        mp1 = _FakeMultipath("mp1", ["r1", "r2"])
        mp2 = _FakeMultipath("mp2", ["r3"])
        transcripts = [
            _FakeTranscript("gene1", "tx1", 100.0, mp1),
            _FakeTranscript("gene2", "tx2", 300.0, mp2),
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
        rpm_total_reads_sum = sum(float(row[-1]) for row in rows)

        assert round(tpm_sum, 3) == 1000000.0
        assert round(rpm_total_reads_sum, 3) == 400000.0
    finally:
        LRAA_Globals.config["num_total_reads"] = old_num_total_reads
        LRAA_Globals.DEBUG = old_debug
