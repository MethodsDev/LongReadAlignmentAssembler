import math

from filter_bam_to_secondary_rescue import (
    AlignmentMetrics,
    _filter_group,
    _metrics_are_complete,
    should_rescue_secondary,
)


def metrics(as_score=100, ms=100, nm_rate=0.01, aligned_query_len=1000, max_intron=0):
    return AlignmentMetrics(
        as_score=float(as_score),
        ms=float(ms),
        nm_rate=float(nm_rate),
        aligned_query_len=int(aligned_query_len),
        max_intron=int(max_intron),
    )


def test_branch1_passes_for_top_secondary_with_tied_ms_and_no_worse_nm_rate():
    primary = metrics(ms=100, nm_rate=0.01, max_intron=100)
    secondary = metrics(ms=99.6, nm_rate=0.01, max_intron=100)
    keep, branch = should_rescue_secondary(primary, secondary, 1)
    assert keep
    assert branch == "branch1"


def test_branch1_fails_for_non_top_secondary():
    primary = metrics(ms=100, nm_rate=0.01)
    secondary = metrics(ms=99.6, nm_rate=0.01)
    keep, branch = should_rescue_secondary(primary, secondary, 2)
    assert not keep
    assert branch == "failed_heuristic"


def test_branch2_passes_for_better_nm_rate_and_larger_secondary_intron():
    primary = metrics(ms=100, nm_rate=0.02, max_intron=100)
    secondary = metrics(ms=90, nm_rate=0.014, max_intron=103)
    keep, branch = should_rescue_secondary(primary, secondary, 3)
    assert keep
    assert branch == "branch2"


def test_branch2_fails_when_intron_ratio_is_not_large_enough():
    primary = metrics(ms=100, nm_rate=0.02, max_intron=100)
    secondary = metrics(ms=90, nm_rate=0.014, max_intron=101)
    keep, branch = should_rescue_secondary(primary, secondary, 3)
    assert not keep
    assert branch == "failed_heuristic"


def test_missing_ms_fails_conservatively():
    primary = metrics(ms=100, nm_rate=0.01)
    secondary = metrics(ms=math.nan, nm_rate=0.01)
    keep, branch = should_rescue_secondary(primary, secondary, 1)
    assert not keep
    assert branch == "missing_metrics"


def test_invalid_aligned_length_fails_completeness_check():
    assert not _metrics_are_complete(metrics(aligned_query_len=0))


class FakeRead:
    def __init__(
        self,
        name,
        *,
        is_secondary=False,
        is_supplementary=False,
        as_score=100,
        ms=100,
        nm=10,
        aligned_query_len=1000,
        max_intron=0,
        ref="chr1",
        start=100,
        flag=0,
    ):
        self.query_name = name
        self.is_unmapped = False
        self.is_supplementary = is_supplementary
        self.is_secondary = is_secondary
        self.reference_name = ref
        self.reference_start = start
        self.flag = flag
        self.cigarstring = f"{aligned_query_len}M"
        self.query_alignment_length = aligned_query_len
        self.cigartuples = [(0, aligned_query_len)]
        if max_intron:
            self.cigarstring = f"{aligned_query_len // 2}M{max_intron}N{aligned_query_len // 2}M"
            self.cigartuples = [
                (0, aligned_query_len // 2),
                (3, max_intron),
                (0, aligned_query_len // 2),
            ]
        self._tags = {"AS": as_score, "ms": ms, "NM": nm}

    def has_tag(self, tag_name):
        return tag_name in self._tags

    def get_tag(self, tag_name):
        return self._tags[tag_name]


class FakeWriter:
    def __init__(self):
        self.records = []

    def write(self, read):
        self.records.append(read)


def test_filter_group_always_retains_primary_and_only_passing_secondary():
    primary = FakeRead("read1", as_score=100, ms=100, nm=10, start=100)
    passing_secondary = FakeRead(
        "read1",
        is_secondary=True,
        as_score=95,
        ms=99.7,
        nm=10,
        start=200,
    )
    failing_secondary = FakeRead(
        "read1",
        is_secondary=True,
        as_score=90,
        ms=99.7,
        nm=12,
        start=300,
    )
    writer = FakeWriter()
    stats = {
        "read_groups": 0,
        "input_records": 0,
        "non_secondary_records": 0,
        "kept_non_secondary": 0,
        "secondary_candidates": 0,
        "kept_secondary_heuristic": 0,
        "kept_secondary_branch1": 0,
        "kept_secondary_branch2": 0,
        "dropped_secondary": 0,
        "dropped_secondary_missing_metrics": 0,
        "dropped_secondary_failed_delta_nm_rate": 0,
        "dropped_secondary_failed_heuristic": 0,
        "dropped_unmapped": 0,
        "dropped_supplementary": 0,
        "groups_with_secondaries_without_primary": 0,
    }

    _filter_group(
        "read1",
        [primary, passing_secondary, failing_secondary],
        writer,
        stats,
    )

    assert writer.records == [primary, passing_secondary]
    assert stats["kept_non_secondary"] == 1
    assert stats["kept_secondary_branch1"] == 1
    assert stats["dropped_secondary_failed_delta_nm_rate"] == 1
