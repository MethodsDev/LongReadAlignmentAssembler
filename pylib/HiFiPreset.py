"""The --HiFi preset: the config values --HiFi sets, in one place.

LRAA applies it in main(), and every chunk worker re-derives it from its own
forwarded --HiFi. merge_LRAA_GTFs.py does not take --HiFi: none of these values
can act in a merge (see the note in its main()). Each function takes the config
dict and the parsed arguments; a caller without LRAA's full argument set passes a
namespace with the fields read here (HiFi, no_infer_TSS, no_infer_PolyA,
aggregate_adjacent_splice_boundaries, user_specified_dests, and
min_alt_splice_freq only when the user set it).
"""

import LRAA_Globals


def _apply_hifi_splice_default(config, args):
    if args.HiFi and not hasattr(args, "min_alt_splice_freq"):
        config["min_alt_splice_freq"] = 0.005


def _apply_hifi_isoform_fraction_defaults(config, args):
    """The other two thirds of the HiFi threshold set.

    These three gates were measured jointly, not separately: each is close to
    inert on its own, and only the combination recovers isoforms. On SIRVs, where
    the truth set is complete within multi-exon chains and every false positive
    is therefore visible, the set gains 4.6 points of sensitivity over 12 samples
    (p=0.016) for a precision change of -0.2 points that does not reach
    significance -- 6.6 true chains gained per false one. Evidence is PacBio HiFi
    throughout, which is why this is a HiFi preset and the library defaults for
    other platforms are left alone.

    min_isoform_fraction carries an argparse default, so explicitness comes from
    user_specified_dests (read from argv, set in main() before any derivation).
    min_frac_gene_unique_reads has no CLI flag at all and is written
    unconditionally; --config_update is applied after the preset, so a user file
    still wins.
    """

    if not args.HiFi:
        return
    specified = getattr(args, "user_specified_dests", None) or set()
    if "min_isoform_fraction" not in specified:
        config["min_isoform_fraction"] = 0.005
    config["min_frac_gene_unique_reads"] = 0.005


def _apply_hifi_config_overrides(config, args):
    """The HiFi preset, as a function so both run paths apply the identical set.

    Chunked mode exits before main()'s inline copy of this ran, so the outer
    driver's config stayed at library defaults while every worker re-derived the
    preset from its own forwarded --HiFi. That worked only because the preset is
    a pure function of args; the moment anything on the chunked path reads a
    HiFi-derived value from config, the two disagree. One function, both callers.
    """

    if not args.HiFi:
        return
    config["HiFi"] = True
    # Enable TSS/PolyA by default in HiFi, but respect user-specified disabling flags
    config["infer_TSS"] = False if args.no_infer_TSS else True
    config["infer_PolyA"] = False if args.no_infer_PolyA else True
    config["min_per_id"] = LRAA_Globals.HIFI_MIN_PER_ID
    config["rescue_unassigned_max_indel_length"] = 10
    config["max_exon_spur_length"] = 13
    config["fracture_splice_graph_at_input_transcript_bounds"] = True
    config["aggregate_adjacent_splice_boundaries"] = False
    if getattr(args, "aggregate_adjacent_splice_boundaries", False):
        config["aggregate_adjacent_splice_boundaries"] = True
    _apply_hifi_splice_default(config, args)
    _apply_hifi_isoform_fraction_defaults(config, args)
