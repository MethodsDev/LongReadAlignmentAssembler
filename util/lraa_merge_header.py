#!/usr/bin/env python3

"""Build the comment header for a file merged from several LRAA outputs.

A sharded run (LRAA.wdl scatters one LRAA invocation per contig) produces no
single command, so the merged output cannot carry one input's "# LRAA CMD:"
line as if it were its own: every shard names its own split genome, split BAM
and --contig.  This reduces the input commands to what they agree on, marking
the arguments that varied with <varies>, and records the contigs covered.

Nothing here is invented.  The emitted command is the input commands with their
differences masked, so it is a template rather than a runnable line -- and it is
labelled as one.  When the inputs disagree structurally, every distinct command
is emitted verbatim instead.
"""

import argparse
import gzip
import shlex
import sys


CMD_PREFIX = "# LRAA CMD: "
VARIES_PLACEHOLDER = "<varies>"

# Substring identifying a tracking file whose rows cover only the reads coverage
# normalization retained. Such a file understates any count derived by summing
# frac_assigned, so the warning has to reach the merged output: this header is
# reconstructed rather than copied, and a dropped warning leaves a file that is
# incomplete with nothing about it to say so. One incomplete input is enough, since
# the merge inherits its gaps.
INCOMPLETE_TRACKING_MARKER = "use_XW_read_weights_for_quant"
# Detection is scoped to a line starting with this prefix, not a blind substring scan
# of the whole leading comment block: every input also carries a "# LRAA CMD: ..."
# echo of its own command line, which contains INCOMPLETE_TRACKING_MARKER whenever
# --use_XW_read_weights_for_quant was passed at all -- including a --stream_reads run,
# whose tracking file is complete and does not carry the marker line. Matching the CMD
# echo would mark that file incomplete on the strength of its own invocation text.
INCOMPLETE_TRACKING_MARKER_LINE_PREFIX = "# WARNING:"
INCOMPLETE_TRACKING_COMMENT = (
    "# WARNING: produced with --use_XW_read_weights_for_quant. Rows cover only reads "
    "retained by coverage normalization, and frac_assigned is NOT weighted by XW. "
    "Counts derived by summing frac_assigned understate the library, unevenly, and most "
    "at high-coverage loci. Do not derive counts from this file."
)


def _open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def leading_comments(path):
    """The unbroken block of '#' lines at the top of a file, newlines stripped."""

    comments = []
    with _open_text(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            comments.append(line.rstrip("\r\n"))
    return comments


def input_argv(path):
    """The argv recorded by an input file, or None when it recorded no command."""

    for comment in leading_comments(path):
        if comment.startswith(CMD_PREFIX):
            return shlex.split(comment[len(CMD_PREFIX) :])
    return None


def _mask_varying_args(argvs):
    """One argv with the values that differ replaced by <varies>, else None.

    Masking is only safe when the commands line up token for token and every
    disagreement falls on an argument value.  A disagreement on a flag means the
    inputs were parameterized differently, not just pointed at different files,
    and collapsing them would hide that.
    """

    if len({len(argv) for argv in argvs}) != 1:
        return None

    masked = []
    for tokens in zip(*argvs):
        if len(set(tokens)) == 1:
            masked.append(shlex.quote(tokens[0]))
        elif any(token.startswith("-") for token in tokens):
            return None
        else:
            masked.append(VARIES_PLACEHOLDER)
    return " ".join(masked)


def _contig_values(argvs):
    contigs = []
    for argv in argvs:
        for i, token in enumerate(argv[:-1]):
            if token == "--contig":
                value = argv[i + 1]
                if value not in contigs:
                    contigs.append(value)
                break
    return contigs


def _marks_incomplete_tracking(path):
    """Whether this input declares that its rows cover only retained reads."""
    try:
        return any(
            line.startswith(INCOMPLETE_TRACKING_MARKER_LINE_PREFIX)
            and INCOMPLETE_TRACKING_MARKER in line
            for line in leading_comments(path)
        )
    except OSError:
        return False


def merge_header_lines(version_comment, input_files):
    lines = [version_comment]

    argvs = [argv for argv in (input_argv(path) for path in input_files) if argv]

    distinct = []
    for argv in argvs:
        if argv not in distinct:
            distinct.append(argv)

    notes = []
    if len(distinct) == 1:
        lines.append(CMD_PREFIX + shlex.join(distinct[0]))
    elif len(distinct) > 1:
        masked = _mask_varying_args(distinct)
        if masked is not None:
            lines.append(CMD_PREFIX + masked)
            notes.append(
                "{} marks arguments that differed between inputs".format(
                    VARIES_PLACEHOLDER
                )
            )
        else:
            for argv in distinct:
                lines.append(CMD_PREFIX + shlex.join(argv))
            notes.append("{} distinct commands".format(len(distinct)))

    contigs = _contig_values(distinct)
    if len(contigs) > 1:
        notes.append("--contig values: " + " ".join(contigs))

    if len(input_files) > 1:
        merge_note = "# LRAA merge: {} inputs".format(len(input_files))
        for note in notes:
            merge_note += "; " + note
        lines.append(merge_note)

    # Carried explicitly because this header is rebuilt, not copied: an input's own
    # warning would otherwise be dropped and the merged file would understate counts
    # with nothing about it to say so.
    if any(_marks_incomplete_tracking(path) for path in input_files):
        lines.append(INCOMPLETE_TRACKING_COMMENT)

    return lines


def main():
    parser = argparse.ArgumentParser(
        description="Emit the comment header for a file merged from LRAA outputs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version_comment",
        type=str,
        required=True,
        help="the '# LRAA version ...' line to lead with",
    )
    parser.add_argument(
        "--inputs",
        type=str,
        required=True,
        nargs="+",
        help="merge inputs carrying LRAA comment headers (plain or gzipped)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="write here instead of stdout",
    )

    args = parser.parse_args()

    lines = merge_header_lines(args.version_comment, args.inputs)

    ofh = open(args.output, "wt") if args.output else sys.stdout
    try:
        for line in lines:
            print(line, file=ofh)
    finally:
        if args.output:
            ofh.close()


if __name__ == "__main__":
    main()
