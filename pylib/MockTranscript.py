#!/usr/bin/env python3

import pytest
from unittest.mock import MagicMock, patch


# --- Minimal Transcript/Feature class for testing --- #
class MockTranscript:
    def __init__(self, transcript_id, chrom, strand, exons):
        self._id = transcript_id
        self._chrom = chrom
        self._strand = strand
        self._exons = exons  # list of (start, end) tuples, sorted
        self._introns = []
        if len(exons) > 1:
            self._introns = [
                (exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)
            ]

    def get_transcript_id(self):
        return self._id

    def get_coords(self):
        return (self._exons[0][0], self._exons[-1][1])

    def get_exon_segments(self):
        return self._exons

    def get_num_exon_segments(self):
        return len(self._exons)

    def has_introns(self):
        return len(self._introns) > 0

    def get_introns(self):
        return self._introns

    def get_introns_string(self, chrom=None):
        return f"{self._chrom}:{self._strand}:{self._introns}"

    def get_exons_string(self):
        return "_".join([f"{e[0]}-{e[1]}" for e in self._exons])

    def get_pretty_alignment_segments(self):
        return self._exons

    def get_pretty_alignment_string(self, chrom=None):
        return self.get_exons_string()
