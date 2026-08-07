import gzip
import shlex
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).with_name("lraa_merge_header.py")

VERSION_COMMENT = "# LRAA version v0.18.0"

SHARD_BODY = 'SIRV1\tLRAA\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1";\n'


def write_shard(path, cmd=None, body=SHARD_BODY, gzipped=False):
    lines = [VERSION_COMMENT]
    if cmd is not None:
        lines.append("# LRAA CMD: " + cmd)
    text = "".join(line + "\n" for line in lines) + body
    opener = gzip.open if gzipped else open
    with opener(path, "wt") as ofh:
        ofh.write(text)
    return str(path)


def run_header(shards):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--version_comment", VERSION_COMMENT, "--inputs"]
        + shards,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.splitlines()


def contig_cmd(contig):
    return (
        "LRAA --genome /shard/{c}.fa --bam /shard/{c}.bam "
        "--output_prefix S.{c} --contig {c} --HiFi".format(c=contig)
    )


def test_shards_differing_only_in_values_collapse_to_one_masked_command(tmp_path):
    shards = [
        write_shard(tmp_path / "chr1.gtf", contig_cmd("chr1")),
        write_shard(tmp_path / "chr2.gtf", contig_cmd("chr2")),
        write_shard(tmp_path / "chrX.gtf", contig_cmd("chrX")),
    ]

    lines = run_header(shards)

    assert lines[0] == VERSION_COMMENT
    assert lines[1] == (
        "# LRAA CMD: LRAA --genome <varies> --bam <varies> "
        "--output_prefix <varies> --contig <varies> --HiFi"
    )
    assert lines[2] == (
        "# LRAA merge: 3 inputs; <varies> marks arguments that differed "
        "between inputs; --contig values: chr1 chr2 chrX"
    )
    assert len(lines) == 3


def test_identical_commands_are_reported_verbatim_and_stay_runnable(tmp_path):
    cmd = "LRAA --genome '/data/my ref/hg38.fa' --bam /data/s.bam --HiFi"
    shards = [
        write_shard(tmp_path / "a.expr", cmd),
        write_shard(tmp_path / "b.expr", cmd),
    ]

    lines = run_header(shards)

    assert lines[1] == "# LRAA CMD: " + cmd
    assert shlex.split(lines[1][len("# LRAA CMD: ") :]) == shlex.split(cmd)
    assert lines[2] == "# LRAA merge: 2 inputs"


def test_a_flag_that_only_some_inputs_carry_blocks_masking(tmp_path):
    """Masking may not hide a parameterization difference behind <varies>."""

    shards = [
        write_shard(tmp_path / "a.gtf", "LRAA --bam /a.bam --contig chr1 --HiFi"),
        write_shard(tmp_path / "b.gtf", "LRAA --bam /b.bam --contig chr2 --no_EM"),
    ]

    lines = run_header(shards)

    assert lines[1] == "# LRAA CMD: LRAA --bam /a.bam --contig chr1 --HiFi"
    assert lines[2] == "# LRAA CMD: LRAA --bam /b.bam --contig chr2 --no_EM"
    assert lines[3] == (
        "# LRAA merge: 2 inputs; 2 distinct commands; --contig values: chr1 chr2"
    )


def test_differing_argument_counts_block_masking(tmp_path):
    shards = [
        write_shard(tmp_path / "a.gtf", "LRAA --bam /a.bam --HiFi"),
        write_shard(tmp_path / "b.gtf", "LRAA --bam /b.bam"),
    ]

    lines = run_header(shards)

    assert lines[1:3] == [
        "# LRAA CMD: LRAA --bam /a.bam --HiFi",
        "# LRAA CMD: LRAA --bam /b.bam",
    ]


def test_input_without_a_command_is_counted_but_contributes_none(tmp_path):
    """An empty shard is stamped with the version line alone; it must not be lost."""

    shards = [
        write_shard(tmp_path / "chr1.gtf", contig_cmd("chr1")),
        write_shard(tmp_path / "chrM.gtf", cmd=None, body=""),
    ]

    lines = run_header(shards)

    assert lines == [
        VERSION_COMMENT,
        "# LRAA CMD: " + contig_cmd("chr1"),
        "# LRAA merge: 2 inputs",
    ]


def test_gzipped_inputs_are_read(tmp_path):
    shards = [
        write_shard(tmp_path / "a.tracking.gz", contig_cmd("chr1"), gzipped=True),
        write_shard(tmp_path / "b.tracking.gz", contig_cmd("chr2"), gzipped=True),
    ]

    lines = run_header(shards)

    assert "--contig <varies>" in lines[1]
    assert lines[2].endswith("--contig values: chr1 chr2")


def test_a_hash_line_below_the_header_block_is_not_read_as_a_command(tmp_path):
    body = SHARD_BODY + "# LRAA CMD: not-the-header\n"
    shards = [
        write_shard(tmp_path / "a.gtf", contig_cmd("chr1"), body=body),
        write_shard(tmp_path / "b.gtf", contig_cmd("chr2"), body=body),
    ]

    lines = run_header(shards)

    assert "not-the-header" not in "\n".join(lines)


def test_single_input_carries_its_command_and_no_merge_note(tmp_path):
    shards = [write_shard(tmp_path / "a.gtf", contig_cmd("chr1"))]

    lines = run_header(shards)

    assert lines == [VERSION_COMMENT, "# LRAA CMD: " + contig_cmd("chr1")]
