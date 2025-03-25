#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import logging
import string
import pysam
from collections import defaultdict


logger = logging.getLogger(__name__)


class GenomeFeature:

    def __init__(self, contig_acc, lend, rend, orient):

        assert orient in ("+", "-"), "Error, orient must be set to + or - "

        self._contig_acc = contig_acc
        self._lend = lend
        self._rend = rend
        self._orient = orient
        self._id = "__id_not_set__"
        self._read_types = set()
        self._weight = (
            0  # for ordering polyA and TSS with shared coordinates with exons
        )

        return

    def get_contig_acc(self):
        return self._contig_acc

    def get_coords(self):
        return (self._lend, self._rend)

    def get_read_support(self):
        # implement your own!
        return 0

    def get_feature_length(self):
        return self._rend - self._lend + 1

    def get_bed_row(self, pad=0):
        return "\t".join(
            [
                str(x)
                for x in [
                    self._contig_acc,
                    self._lend - pad - 1,
                    self._rend + pad,
                    self._id
                    + "[{}]".format(self._orient)
                    + ":{}".format(self.get_read_support()),
                    self.get_read_support(),
                    self._orient,
                ]
            ]
        )

    def get_id(self):
        return self._id

    def add_read_type(self, read_type: str):
        self._read_types.add(read_type)

    def add_read_types(self, read_types: list):
        for read_type in read_types:
            self.add_read_type(read_type)

    def has_read_type(self, read_type):
        return read_type in self._read_types

    def get_read_types(self):
        return list(self._read_types)

    def get_orient(self):
        return self._orient

    def get_coord_string(self):
        lend, rend = self.get_coords()
        return "{}:({})[({}, {}]".format(
            self.get_contig_acc(), self.get_strand(), lend, rend
        )


class Intron(GenomeFeature):

    intron_id_counter = 0

    splice_dinucs_top_strand = {"GTAG", "GCAG", "ATAC"}
    splice_dinucs_bottom_strand = {
        "CTAC",
        "CTGC",
        "GTAT",
    }  # revcomp of top strand dinucs

    def __init__(self, contig_acc, lend, rend, orient, count):
        super().__init__(contig_acc, lend, rend, orient)

        self._count = count

        Intron.intron_id_counter += 1
        self._id = "I:{}".format(Intron.intron_id_counter)

        return

    def get_read_support(self):
        return self._count

    def __repr__(self):
        return "Intron: {} {}:{}-{} [{}] count:{} rtypes:{}".format(
            self._id,
            self._contig_acc,
            self._lend,
            self._rend,
            self._orient,
            self._count,
            self.get_read_types(),
        )

    # static methods
    def check_canonical_splicing(intron_lend, intron_rend, contig_seq):

        contig_seq_len = len(contig_seq)
        # print(f"intron-lend: {intron_lend}, intron-rend: {intron_rend}, contig_seq_len: {contig_seq_len}")

        if (intron_lend >= contig_seq_len - 2) or (intron_rend >= contig_seq_len - 2):
            # out of bounds
            return None

        dinuc_left = contig_seq[intron_lend - 1] + contig_seq[intron_lend - 1 + 1]
        dinuc_right = contig_seq[intron_rend - 1 - 1] + contig_seq[intron_rend - 1]
        dinuc_combo = dinuc_left + dinuc_right

        # print("INTRON: {}-{}:{}".format(intron_lend, intron_rend, dinuc_combo))

        if dinuc_combo in Intron.splice_dinucs_top_strand:
            return "+"
        elif dinuc_combo in Intron.splice_dinucs_bottom_strand:
            return "-"
        else:
            return None


class Exon(GenomeFeature):

    exon_id_counter = 0

    def __init__(self, contig_acc, lend, rend, orient, mean_coverage):
        super().__init__(contig_acc, lend, rend, orient)

        self._mean_coverage = mean_coverage

        Exon.exon_id_counter += 1
        self._id = "E:{}".format(Exon.exon_id_counter)

        return

    def get_read_support(self):
        return self._mean_coverage

    def __repr__(self):
        return "Exon: {} {}:{}-{} [{}] mean_cov:{}".format(
            self._id,
            self._contig_acc,
            self._lend,
            self._rend,
            self._orient,
            self._mean_coverage,
        )

    @classmethod
    def reset_counter(cls):
        cls.exon_id_counter = 0


class TSS(GenomeFeature):

    TSS_id_counter = 0

    def __init__(self, contig_acc, lend, rend, orient, count):
        super().__init__(contig_acc, lend, rend, orient)

        self._count = count
        self._weight = -1 if orient == "+" else 1

        TSS.TSS_id_counter += 1
        self._id = "TSS:{}".format(TSS.TSS_id_counter)

        return

    def get_read_support(self):
        # implement your own!
        return self._count

    def __repr__(self):
        return "TSS: {} {}:{}-{} [{}] read_count:{}".format(
            self._id,
            self._contig_acc,
            self._lend,
            self._rend,
            self._orient,
            self._count,
        )


class PolyAsite(GenomeFeature):

    polyA_counter = 0

    def __init__(self, contig_acc, lend, rend, orient, count):
        super().__init__(contig_acc, lend, rend, orient)

        self._count = count
        self._weight = 1 if orient == "+" else -1

        PolyAsite.polyA_counter += 1
        self._id = "POLYA:{}".format(PolyAsite.polyA_counter)

        return

    def get_read_support(self):
        # implement your own!
        return self._count

    def __repr__(self):
        return "PolyAsite: {} {}:{}-{} [{}] read_count:{}".format(
            self._id,
            self._contig_acc,
            self._lend,
            self._rend,
            self._orient,
            self._count,
        )
