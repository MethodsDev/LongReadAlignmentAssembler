#!/usr/bin/env python3

"""Chunk first, split the orientations inside the chunk.

What these defend, in the order they would hurt:

ORDER. ``separate_bam_by_strand`` rewrites ``read.is_reverse`` when the
orientation it infers disagrees with the aligner, and the extractor's strand
filter reads the raw flag. Extract-then-split is the only order that gives the
right answer, and the wrong order gives a WRONG one that looks entirely normal
-- no crash, no empty file, just reads in the wrong chunk. So the refusal is
tested, not the happy path alone.

ACCOUNTING. The parity control is the whole bam MINUS the reads cut selection
named, so those names have to be exactly the reads extraction dropped. Under
strandless cuts one position severs reads on both orientations and the set is a
union rather than a per-orientation set, so the identity is checked on every run
and the check itself is tested here from both sides.

PARTITION. A strandless chunk feeds two quant units. Every record and every
transcript model it holds must reach exactly one of them: a lost record shrinks
the chunked arm against its control, a lost model drops a row from the merged
table, and a model reaching BOTH units doubles the table -- which it did, until
the annotation was partitioned alongside the bam.

MODES. One output directory serves one mode. The per-stage sentinels differ
between the modes by construction, which is exactly what would let a second mode
overwrite the first's merged output while reporting every stage as new work.
"""

import os
import sys

import pysam
import pytest

REPO = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
)
sys.path.insert(0, os.path.join(REPO, "pylib"))

import ChunkedRun  # noqa: E402  (path insert must precede the import)


# --------------------------------------------------------------------- fixtures


def write_bam(path, num_forward, num_reverse):
    """A bam of ``num_forward`` forward and ``num_reverse`` reverse records."""

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chunk", "LN": 100000}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as ofh:
        for i in range(num_forward + num_reverse):
            aln = pysam.AlignedSegment(ofh.header)
            aln.query_name = "read{}".format(i)
            aln.reference_id = 0
            aln.reference_start = 10 * i
            aln.mapping_quality = 60
            aln.cigarstring = "50M"
            aln.query_sequence = "A" * 50
            aln.flag = 0 if i < num_forward else 16
            ofh.write(aln)
    return path


def gtf_line(chrom, feature, lend, rend, strand, gene_id, transcript_id=None):
    attrs = 'gene_id "{}";'.format(gene_id)
    if transcript_id:
        attrs += ' transcript_id "{}";'.format(transcript_id)
    return "\t".join(
        [chrom, "test", feature, str(lend), str(rend), ".", strand, ".", attrs]
    ) + "\n"


def make_chunk(tmp_path, strand="", emitted=(0, 0), transcripts=(("+", 1), ("-", 1))):
    """A chunk record as stage 3 builds one, with its mini GTF on disk.

    ``emitted`` is (forward, reverse) as the extractor would have counted them.
    """

    cdir = tmp_path / "chunks" / "chrT_00"
    cdir.mkdir(parents=True, exist_ok=True)
    prefix = str(cdir / "chunk")
    forward, reverse = emitted

    with open(prefix + ".gtf", "wt") as ofh:
        index = 0
        for orientation, count in transcripts:
            for _ in range(count):
                index += 1
                gene = "G{}".format(index)
                tid = "T{}".format(index)
                lend, rend = 100 * index, 100 * index + 50
                ofh.write(gtf_line("chunk", "gene", lend, rend, orientation, gene))
                ofh.write(
                    gtf_line(
                        "chunk", "transcript", lend, rend, orientation, gene, tid
                    )
                )
                ofh.write(
                    gtf_line("chunk", "exon", lend, rend, orientation, gene, tid)
                )
    num_transcripts = sum(count for _, count in transcripts)

    manifest = {
        "strand": strand or None,
        "strand_split_required": not strand,
        "offset": 0,
        "window_origin": 0,
        "partition_lend": 1,
        "partition_rend": 100000,
        "dropped_read_names": [],
        "counts": {
            "alignments_emitted": forward + reverse,
            "alignments_emitted_forward": forward,
            "alignments_emitted_reverse": reverse,
            "alignments_dropped_overhang": 0,
            "gtf_transcripts_emitted": num_transcripts,
        },
    }
    return {
        "chunk_id": "chrT_00",
        "chrom": "chrT",
        "strand": strand,
        "strandless": not strand,
        "region": "chrT:1-100000",
        "index": 0,
        "order": 0,
        "dir": str(cdir),
        "prefix": prefix,
        "log": str(cdir / "chunk.log"),
        "manifest": manifest,
        "offset": 0,
        "window_origin": 0,
        "upstream_token": "stage3.up_aaaaaaaaaaaa",
        "units": ChunkedRun.chunk_quant_units(
            "chrT_00", str(cdir), prefix, strand, 0, 0
        ),
    }


def selection(chrom, *spans):
    return {
        "chrom": chrom,
        "segments": [{"lend": lend, "rend": rend} for lend, rend in spans],
    }


# ---------------------------------------------------------------- the ordering


def test_a_chunk_extracted_for_one_orientation_is_never_split(tmp_path):
    """The mistake that produces a plausible wrong answer, refused.

    Its bam was filtered on the RAW flag, so any read the split would have
    flipped is already in the wrong chunk and splitting cannot move it back.
    """

    chunk = make_chunk(tmp_path, strand="+", emitted=(3, 0))
    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.assert_extracted_strandlessly(chunk)
    assert "REFUSING" in str(err.value)


def test_a_chunk_already_strand_separated_is_never_split_again(tmp_path):
    """The extractor states it directly; splitting twice empties an orientation."""

    chunk = make_chunk(tmp_path, emitted=(3, 2))
    chunk["manifest"]["strand_split_required"] = False
    with pytest.raises(ChunkedRun.PipelineError):
        ChunkedRun.assert_extracted_strandlessly(chunk)


def test_a_strandless_chunk_is_accepted_for_splitting(tmp_path):
    chunk = make_chunk(tmp_path, emitted=(3, 2))
    ChunkedRun.assert_extracted_strandlessly(chunk)  # must not raise


def test_a_manifest_without_the_key_is_judged_on_its_strand(tmp_path):
    """Absence is not evidence of the wrong order; the strand says the same thing."""

    chunk = make_chunk(tmp_path, emitted=(3, 2))
    del chunk["manifest"]["strand_split_required"]
    ChunkedRun.assert_extracted_strandlessly(chunk)


# ------------------------------------------------------------ split accounting


def test_a_split_that_partitions_exactly_is_accepted(tmp_path):
    chunk = make_chunk(tmp_path, emitted=(4, 3))
    prefix = chunk["prefix"] + ".strand"
    write_bam(prefix + ".+.bam", 4, 0)
    write_bam(prefix + ".-.bam", 0, 3)

    counts = ChunkedRun.verify_chunk_split(chunk, prefix)
    assert counts["records_plus"] == 4
    assert counts["records_minus"] == 3
    assert counts["records_total"] == 7
    assert counts["records_lost"] == 0


def test_a_split_that_loses_a_record_fails_the_chunk(tmp_path):
    """Extraction already applied the split's own filter, so a loss here is a
    filter divergence -- and it would shrink the chunked arm against its
    control without any other symptom."""

    chunk = make_chunk(tmp_path, emitted=(4, 3))
    prefix = chunk["prefix"] + ".strand"
    write_bam(prefix + ".+.bam", 4, 0)
    write_bam(prefix + ".-.bam", 0, 2)

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.verify_chunk_split(chunk, prefix)
    assert "lost 1 record" in str(err.value)


def test_records_moving_between_orientations_fails_the_chunk(tmp_path):
    """Totals can agree while the partition is wrong. Only inference moves a
    read between the two bams, and this pipeline never asks for inference."""

    chunk = make_chunk(tmp_path, emitted=(4, 3))
    prefix = chunk["prefix"] + ".strand"
    write_bam(prefix + ".+.bam", 5, 0)
    write_bam(prefix + ".-.bam", 0, 2)

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.verify_chunk_split(chunk, prefix)
    assert "which records are which" in str(err.value)


def test_a_missing_orientation_bam_fails_the_chunk(tmp_path):
    chunk = make_chunk(tmp_path, emitted=(1, 1))
    prefix = chunk["prefix"] + ".strand"
    write_bam(prefix + ".+.bam", 1, 0)

    with pytest.raises(ChunkedRun.PipelineError):
        ChunkedRun.verify_chunk_split(chunk, prefix)


# ------------------------------------------------------- annotation partition


def test_the_annotation_is_partitioned_by_orientation(tmp_path):
    """Every model reaches exactly one unit. Both units sharing the chunk's GTF
    is what doubled the merged table: 1,110 rows where 555 were right."""

    chunk = make_chunk(tmp_path, emitted=(1, 1), transcripts=(("+", 3), ("-", 2)))
    prefix = chunk["prefix"] + ".strand"
    counts = ChunkedRun.split_chunk_gtf_by_strand(chunk, prefix)

    assert counts["gtf_transcripts_plus"] == 3
    assert counts["gtf_transcripts_minus"] == 2

    for orientation, expected in (("+", 3), ("-", 2)):
        with open("{}.{}.gtf".format(prefix, orientation)) as fh:
            lines = [line.split("\t") for line in fh if line.strip()]
        assert lines, orientation
        assert {line[6] for line in lines} == {orientation}
        assert sum(1 for line in lines if line[2] == "transcript") == expected


def test_an_orientation_with_no_models_gets_an_empty_annotation(tmp_path):
    """Normal, not exceptional: the strand-first arm already produces chunks
    with nothing annotated on one orientation."""

    chunk = make_chunk(tmp_path, emitted=(1, 0), transcripts=(("+", 2),))
    prefix = chunk["prefix"] + ".strand"
    counts = ChunkedRun.split_chunk_gtf_by_strand(chunk, prefix)

    assert counts["gtf_transcripts_minus"] == 0
    assert os.path.exists(prefix + ".-.gtf")
    assert os.path.getsize(prefix + ".-.gtf") == 0


def test_a_model_lost_in_the_partition_fails_the_chunk(tmp_path):
    """A dropped row in the merged table is not an error anyone would notice."""

    chunk = make_chunk(tmp_path, emitted=(1, 1), transcripts=(("+", 2), ("-", 1)))
    chunk["manifest"]["counts"]["gtf_transcripts_emitted"] = 4
    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.split_chunk_gtf_by_strand(chunk, chunk["prefix"] + ".strand")
    assert "annotation split accounting" in str(err.value)


# ------------------------------------------------------------- severed reads


def test_severed_accounting_counts_reads_not_mentions(tmp_path):
    """A severed read is listed by BOTH neighbouring chunks, so the per-chunk
    sum is about twice the read count. The control subtracts READS."""

    cut_dir = tmp_path / "cuts"
    cut_dir.mkdir()
    (cut_dir / "strandless.dropped_reads.txt").write_text("readA\nreadB\n")

    left = make_chunk(tmp_path, emitted=(1, 1))
    right = dict(left, chunk_id="chrT_01")
    left["manifest"] = dict(left["manifest"], dropped_read_names=["readA", "readB"])
    right["manifest"] = dict(right["manifest"], dropped_read_names=["readA", "readB"])

    counts = ChunkedRun.verify_severed_accounting(str(cut_dir), [left, right])
    assert counts["severed_reads"] == 2
    assert counts["per_chunk_drop_mentions"] == 4
    assert counts["sets_identical"] is True


def test_a_read_dropped_but_not_named_fails_the_run(tmp_path):
    """The control would keep a record the chunks never saw, and every
    difference downstream would be confounded by exactly that read."""

    cut_dir = tmp_path / "cuts"
    cut_dir.mkdir()
    (cut_dir / "strandless.dropped_reads.txt").write_text("readA\n")

    chunk = make_chunk(tmp_path, emitted=(1, 1))
    chunk["manifest"] = dict(
        chunk["manifest"], dropped_read_names=["readA", "readB"]
    )

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.verify_severed_accounting(str(cut_dir), [chunk])
    assert "readB" in str(err.value)


def test_a_read_named_but_not_dropped_fails_the_run(tmp_path):
    cut_dir = tmp_path / "cuts"
    cut_dir.mkdir()
    (cut_dir / "strandless.dropped_reads.txt").write_text("readA\nreadB\n")

    chunk = make_chunk(tmp_path, emitted=(1, 1))
    chunk["manifest"] = dict(chunk["manifest"], dropped_read_names=["readA"])

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.verify_severed_accounting(str(cut_dir), [chunk])
    assert "readB" in str(err.value)


# ----------------------------------------------------------------- the units


def test_a_strandless_chunk_feeds_two_units_with_their_own_inputs(tmp_path):
    units = ChunkedRun.chunk_quant_units("chr1_00", "/d", "/d/chunk", "", 7, 3)

    assert [u["unit_id"] for u in units] == ["chr1_00_plus", "chr1_00_minus"]
    assert [u["strand"] for u in units] == ["+", "-"]
    assert len({u["bam"] for u in units}) == 2
    assert len({u["gtf"] for u in units}) == 2
    assert len({u["quant_prefix"] for u in units}) == 2
    assert all(u["offset"] == 7 and u["order"] == 3 for u in units)
    # both read what stage 3b writes, never the unsplit chunk bam
    assert all(".strand." in u["bam"] and ".strand." in u["gtf"] for u in units)


def test_a_strand_first_chunk_keeps_the_paths_it_has_always_had(tmp_path):
    """A moved path or sentinel would silently re-run every existing directory."""

    units = ChunkedRun.chunk_quant_units("chr1_plus_00", "/d", "/d/chunk", "+", 0, 0)

    assert len(units) == 1
    unit = units[0]
    assert unit["unit_id"] == "chr1_plus_00"
    assert unit["sentinel_id"] == "chr1_plus_00"
    assert unit["bam"] == "/d/chunk.bam"
    assert unit["gtf"] == "/d/chunk.gtf"
    assert unit["norm_bam"] == "/d/chunk.norm.bam"
    assert unit["quant_name"] == "chunk_quant"


def test_strandless_sentinels_cannot_be_read_as_strand_first_ones():
    strandless = ChunkedRun.chunk_quant_units("chr1_00", "/d", "/d/chunk", "", 0, 0)
    assert [u["sentinel_id"] for u in strandless] == [
        "strandless_chr1_00_plus",
        "strandless_chr1_00_minus",
    ]


def test_units_merge_in_the_strand_first_order():
    """Stage 6 concatenates + then -, so the two modes' tables line up."""

    chunks = []
    for order, chunk_id in enumerate(("chr1_00", "chr1_01")):
        chunks.append(
            {
                "units": ChunkedRun.chunk_quant_units(
                    chunk_id, "/d", "/d/chunk", "", 0, order
                )
            }
        )
    assert [u["unit_id"] for u in ChunkedRun.ordered_units(chunks)] == [
        "chr1_00_plus",
        "chr1_01_plus",
        "chr1_00_minus",
        "chr1_01_minus",
    ]


# ------------------------------------------------------- cut sources and plan


def test_strandless_cut_selection_reads_the_raw_bam(tmp_path):
    args = ChunkedRun.default_args(bam=str(tmp_path / "raw.bam"), strandless_chunks=True)
    sources = ChunkedRun.cut_sources(args, None, "inputs.up_1", "split.up_1")

    assert len(sources) == 1
    key, tag, bam, parent = sources[0]
    assert key == ""
    assert tag == ChunkedRun.STRANDLESS_TAG
    assert bam == os.path.abspath(str(tmp_path / "raw.bam"))
    # chained on the inputs, because no split has happened
    assert parent == "inputs.up_1"


def test_strand_first_cut_selection_reads_the_split_bams():
    args = ChunkedRun.default_args(bam="raw.bam")
    strand_bams = {"+": "p.bam", "-": "m.bam"}
    sources = ChunkedRun.cut_sources(args, strand_bams, "inputs.up_1", "split.up_1")

    assert [s[0] for s in sources] == ["+", "-"]
    assert [s[2] for s in sources] == ["p.bam", "m.bam"]
    assert {s[3] for s in sources} == {"split.up_1"}


def test_a_strandless_region_carries_no_orientation():
    """The extractor keeps both orientations for a strandless region, and
    refuses a strand-suffixed one over a bam that still holds both."""

    sources = [("", ChunkedRun.STRANDLESS_TAG, "raw.bam", "tok")]
    planned = list(
        ChunkedRun.planned_chunks(sources, {"": [selection("chr1", (1, 10), (11, 20))]})
    )

    assert [p["chunk_id"] for p in planned] == ["chr1_00", "chr1_01"]
    assert [p["region"] for p in planned] == ["chr1:1-10", "chr1:11-20"]


def test_a_strand_first_region_still_carries_its_orientation():
    sources = [("+", "plus", "p.bam", "tok"), ("-", "minus", "m.bam", "tok")]
    selections = {"+": [selection("chr1", (1, 10))], "-": [selection("chr1", (1, 20))]}
    planned = list(ChunkedRun.planned_chunks(sources, selections))

    assert [p["chunk_id"] for p in planned] == ["chr1_plus_00", "chr1_minus_00"]
    assert [p["region"] for p in planned] == ["chr1+:1-10", "chr1-:1-20"]
    assert [p["order"] for p in planned] == [0, 1]


def test_the_plan_reports_intervals_and_the_split_inside_each():
    args = ChunkedRun.default_args(strandless_chunks=True, dry_run=True)
    sources = [("", ChunkedRun.STRANDLESS_TAG, "raw.bam", "tok")]
    selections = {"": [selection("chr1", (1, 1000000), (1000001, 2000000))]}

    plan = ChunkedRun.format_plan(args, ChunkedRun.STRANDLESS_MODE, sources, selections)

    assert "SKIPPED" in plan  # stage 1
    assert "strand split IN CHUNK" in plan
    assert "2 intervals, 2 extraction(s), 4 quant unit(s)" in plan
    assert "chr1_00_plus, chr1_00_minus" in plan


def test_the_strand_first_plan_does_not_claim_an_in_chunk_split():
    args = ChunkedRun.default_args()
    sources = [("+", "plus", "p.bam", "tok"), ("-", "minus", "m.bam", "tok")]
    selections = {"+": [selection("chr1", (1, 10))], "-": [selection("chr1", (1, 10))]}

    plan = ChunkedRun.format_plan(args, ChunkedRun.STRAND_FIRST_MODE, sources, selections)

    assert "1 run over the whole bam" in plan
    assert "not run -- the whole bam was split at stage 1" in plan
    assert "2 contig-strand chunks, 2 extraction(s), 2 quant unit(s)" in plan


# ----------------------------------------------------------------- the modes


def test_one_output_directory_serves_one_mode(tmp_path):
    """The sentinels differ between the modes, which is what would make a
    collision quiet: every stage reruns and overwrites a merged table that
    looked finished."""

    ckpt = ChunkedRun.Checkpoints(str(tmp_path / "__ckpt"))
    ChunkedRun.claim_pipeline_mode(ckpt, ChunkedRun.STRANDLESS_MODE)

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.claim_pipeline_mode(ckpt, ChunkedRun.STRAND_FIRST_MODE)
    assert "already served the strandless pipeline" in str(err.value)


def test_re_running_the_same_mode_is_fine(tmp_path):
    ckpt = ChunkedRun.Checkpoints(str(tmp_path / "__ckpt"))
    ChunkedRun.claim_pipeline_mode(ckpt, ChunkedRun.STRAND_FIRST_MODE)
    ChunkedRun.claim_pipeline_mode(ckpt, ChunkedRun.STRAND_FIRST_MODE)


def test_a_strandless_chunked_run_will_not_invent_its_denominator(
    tmp_path, monkeypatch
):
    """-N is the one number stage 6 does not rebase. Stage 1 is not run, so
    there is nothing to default to and a guess would read as a quantification
    difference against a strand-first run of the same substrate.

    ``run`` marks the environment so a chunk worker cannot re-enter the
    pipeline, and that marker outlives the call. Declared through monkeypatch so
    pytest takes it back off: left set, it makes every later test in this
    process -- and every LRAA it launches -- believe it is a chunk worker.
    """

    monkeypatch.setenv(ChunkedRun.WORKER_ENV, "1")
    for name in ("in.bam", "in.fa", "in.gtf"):
        (tmp_path / name).write_text(name)
    args = ChunkedRun.default_args(
        bam=str(tmp_path / "in.bam"),
        genome_fa=str(tmp_path / "in.fa"),
        gtf=str(tmp_path / "in.gtf"),
        output_dir=str(tmp_path / "out"),
        strandless_chunks=True,
        cpu_budget=1,
    )
    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.run(args)
    assert "-N is required" in str(err.value)


# ------------------------------------------------------------ the worker wiring


class Recorder:
    """Stands in for run_step, and writes what the real tool would have.

    The worker's own checks have to keep biting -- a stub that produced nothing
    would make every assertion below vacuous.
    """

    def __init__(self, emitted):
        self.calls = []
        self.emitted = emitted

    def __call__(self, name, cmd, log, cwd, rss_interval, append=True):
        self.calls.append((name, list(cmd)))
        if name.startswith("stage3b"):
            prefix = cmd[cmd.index("--output_prefix") + 1]
            write_bam(prefix + ".+.bam", self.emitted[0], 0)
            write_bam(prefix + ".-.bam", 0, self.emitted[1])
        elif name.startswith("stage4"):
            out = cmd[cmd.index("--output_bam") + 1]
            write_bam(out, 1, 1)
        elif name.startswith("stage5"):
            out = os.path.join(cwd, cmd[cmd.index("--output_prefix") + 1])
            base = out + "." + ChunkedRun.LRAA_QUANT_ONLY_SUFFIX
            open(base + ".quant.expr", "wt").close()
            open(base + ".quant.tracking.gz", "wt").close()
        return {"step": name, "wall_s": 0.0, "peak_tree_rss_kb": 0}


def run_worker(tmp_path, monkeypatch, chunk):
    recorder = Recorder(
        (
            chunk["manifest"]["counts"]["alignments_emitted_forward"],
            chunk["manifest"]["counts"]["alignments_emitted_reverse"],
        )
    )
    monkeypatch.setattr(ChunkedRun, "run_step", recorder)
    args = ChunkedRun.default_args(strandless_chunks=not chunk["strand"])
    ChunkedRun.chunk_worker(
        args,
        ChunkedRun.Checkpoints(str(tmp_path / "__ckpt")),
        str(tmp_path),
        chunk,
        1000,
        0.5,
        1,
    )
    return recorder


def test_the_split_runs_after_extraction_and_before_everything_else(
    tmp_path, monkeypatch
):
    chunk = make_chunk(tmp_path, emitted=(4, 3))
    recorder = run_worker(tmp_path, monkeypatch, chunk)
    # step names carry the unit they ran for; the stage is the leading token
    steps = [name.split("_", 2)[0] for name, _ in recorder.calls]

    assert steps == ["stage3b", "stage4", "stage5", "stage4", "stage5"]
    assert recorder.calls[0][0] == "stage3b_strand_split_chrT_00"
    assert [name for name, _ in recorder.calls[1:]] == [
        "stage4_normalize_strandless_chrT_00_plus",
        "stage5_quant_strandless_chrT_00_plus",
        "stage4_normalize_strandless_chrT_00_minus",
        "stage5_quant_strandless_chrT_00_minus",
    ]


def test_the_in_chunk_split_is_the_command_stage_1_runs(tmp_path, monkeypatch):
    """Same tool, same filters. A per-chunk difference in filtering would be a
    difference in the record set, not in scheduling -- and inference would
    reassign reads across the two bams, which nothing downstream expects."""

    chunk = make_chunk(tmp_path, emitted=(4, 3))
    recorder = run_worker(tmp_path, monkeypatch, chunk)
    _, cmd = recorder.calls[0]

    assert ChunkedRun.SEPARATE_BAM in cmd
    assert "--infer_read_orient" not in cmd
    assert cmd[cmd.index("--bam") + 1] == chunk["prefix"] + ".bam"
    assert "--max_intron_length" in cmd


def test_the_two_units_share_one_mini_contig_and_split_the_rest(
    tmp_path, monkeypatch
):
    """The saving: sequence has no orientation, so it is extracted once. The
    annotation does, so it is not."""

    chunk = make_chunk(tmp_path, emitted=(4, 3), transcripts=(("+", 2), ("-", 1)))
    recorder = run_worker(tmp_path, monkeypatch, chunk)
    quant_cmds = [cmd for name, cmd in recorder.calls if name.startswith("stage5")]

    assert len(quant_cmds) == 2
    genomes = {cmd[cmd.index("--genome") + 1] for cmd in quant_cmds}
    gtfs = {cmd[cmd.index("--gtf") + 1] for cmd in quant_cmds}
    bams = {cmd[cmd.index("--bam") + 1] for cmd in quant_cmds}

    assert genomes == {chunk["prefix"] + ".fa"}
    assert len(gtfs) == 2
    assert len(bams) == 2


def test_a_strand_first_chunk_still_runs_two_steps_and_no_split(
    tmp_path, monkeypatch
):
    chunk = make_chunk(tmp_path, strand="+", emitted=(4, 0))
    recorder = run_worker(tmp_path, monkeypatch, chunk)

    assert [name.rsplit("_chrT", 1)[0] for name, _ in recorder.calls] == [
        "stage4_normalize",
        "stage5_quant",
    ]
    quant_cmd = recorder.calls[-1][1]
    assert quant_cmd[quant_cmd.index("--gtf") + 1] == chunk["prefix"] + ".gtf"


def test_the_split_accounting_reaches_the_run_record(tmp_path, monkeypatch):
    chunk = make_chunk(tmp_path, emitted=(4, 3))
    run_worker(tmp_path, monkeypatch, chunk)

    rolled = ChunkedRun.roll_up_split_accounting([chunk])
    assert rolled == {
        "intervals_split": 1,
        "alignments_emitted": 7,
        "records_plus": 4,
        "records_minus": 3,
        "records_quantified": 7,
        "records_lost_in_split": 0,
    }
