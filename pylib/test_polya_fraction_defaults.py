#!/usr/bin/env python3

import LRAA_Globals


def test_min_polya_fraction_defaults_to_min_isoform_fraction():
    assert LRAA_Globals.resolve_min_polya_iso_fraction(0.01) == 0.01


def test_explicit_min_polya_fraction_override_is_preserved():
    assert LRAA_Globals.resolve_min_polya_iso_fraction(0.01, 0.05) == 0.05
