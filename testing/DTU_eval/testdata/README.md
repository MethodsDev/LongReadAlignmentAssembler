# DTU_eval test inputs

A 999-transcript slice of a single-cell run over 12 clusters, paired across all
three files by `transcript_id`. Nothing regenerates these; they are checked in.

## `num_exons` in `gene_transcript_splicehashcode.tsv.wAnnotIDs`

The column was added after the fact and does **not** carry the original
per-transcript exon counts, which this fixture's source never recorded. It
carries the one distinction that is exactly recoverable, and the only one any
consumer reads.

`Quantify.py` writes a transcript's splice hash code as

```python
transcript_splice_hash_code = (
    Util_funcs.get_hash_code(transcript.get_introns_string())
    if num_exons > 1
    else transcript_id
)
```

so `transcript_splice_hash_code == transcript_id` iff the transcript is
monoexonic — an identity, not an inference, and it holds on both the original
and the annotated ID columns for all 999 rows. Monoexonic rows therefore carry a
true `1`; spliced rows carry `2`, the smallest count consistent with a hashed
splice pattern.

`--ignore_unspliced` is the only consumer, and it branches on `ne(1)`
(`sc_pseudobulk_test_isoform_DiffUsage.py:219`), so the floor is exact for the
behaviour under test: 689 rows are dropped and 310 retained. **Do not read the
`2`s as real exon counts**, and do not add a test that depends on their
magnitude. `TestDTUFixtureIntegrity` in `pylib/test_diff_iso_tests.py` pins the
column to the hash-code identity so the two cannot drift apart.
