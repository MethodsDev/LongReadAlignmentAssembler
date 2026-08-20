test: test_chunked_single_cell_forwarding
	../../../LRAA --genome ../data/chr19.genome.fa --bam ../data/chr19.bam --HiFi --verbose --no_cleanup --gtf ../data/chr19.gtf  --quant_only --no_chunk --no_stream_reads
	../../../util/sc/singlecell_tracking_to_sparse_matrix.py --tracking LRAA.LRAA.quant-only.quant.tracking.gz --output_prefix test
	../../../util/sc/incorporate_gene_symbols_in_sc_features.py  --ref_gtf ../data/chr19.gtf --id_mappings test.gene_transcript_splicehashcode.tsv   --sparseM_dirs 'test^gene-sparseM' 'test^isoform-sparseM' 'test^splice_pattern-sparseM'


# C2: --chunk forwards single-cell settings through to its chunk workers
# instead of silently dropping them (CodeBacklog.chunking-single-cell.md, "the
# missing forward"). A cut geometry no contig in this fixture can exceed keeps
# this a validation recipe in the C1 sense -- zero cuts, so a difference has
# nothing to be attributed to but the forward itself.
#
# Default --strandless_chunks, not --chunk_by_strand: this GTF used to trip an
# unrelated chunk-worker bug under the strandless split (annotation split
# accounting mismatch, 393+256 != 654 transcript lines, reproducible with no
# --cell_list at all) -- fixed in c8c09f2. Running the production-default mode
# here now also exercises that path rather than dodging it.
#
# The roster is HALF the bam's own distinct barcodes, built here rather than
# committed, so it can never drift from the bam it filters. A roster of ALL
# barcodes would prove nothing about --cell_list actually restricting anything;
# a roster built independently of the bam could name barcodes that were never
# going to appear in the output either way.
test_chunked_single_cell_forwarding:
	rm -rf __sc_unchunked* __sc_chunked* __cell_roster.txt
	python3 -c 'import pysam; bam = pysam.AlignmentFile("../data/chr19.bam"); cbs = sorted({r.get_tag("CB") for r in bam.fetch(until_eof=True) if r.has_tag("CB") and r.has_tag("XM")}); half = cbs[: len(cbs) // 2]; assert half, "no barcodes found"; open("__cell_roster.txt", "w").write("\n".join(half) + "\n"); print("roster: %d of %d distinct barcodes" % (len(half), len(cbs)))'
	../../../LRAA --genome ../data/chr19.genome.fa --bam ../data/chr19.bam --gtf ../data/chr19.gtf \
	    --quant_only --HiFi --cell_list __cell_roster.txt --no_chunk --no_stream_reads --output_prefix __sc_unchunked
	../../../LRAA --chunk --approx_MB_per_cut 100 \
	    --genome ../data/chr19.genome.fa --bam ../data/chr19.bam --gtf ../data/chr19.gtf \
	    --quant_only --HiFi --cell_list __cell_roster.txt --no_stream_reads --output_prefix __sc_chunked
	python3 -c 'import json; d = json.load(open("__sc_chunked.chunked_work/timing.json")); cp = d["cut_placement"]; assert cp["cuts_placed"] == 0, "cuts were placed: %r" % (cp,); assert cp["alignments_severed"] == 0, "alignments were severed: %r" % (cp,); print("structural invariants hold: %d chunks, 0 cuts, 0 severed" % cp["chunks"])'
	python3 -c 'a = {(f[0], f[1]): f[2:] for f in (ln.rstrip("\n").split("\t") for ln in open("__sc_unchunked.LRAA.quant-only.quant.expr") if not ln.startswith("#")) if f[0] != "gene_id"}; b = {(f[0], f[1]): f[2:] for f in (ln.rstrip("\n").split("\t") for ln in open("__sc_chunked.LRAA.quant-only.quant.expr") if not ln.startswith("#")) if f[0] != "gene_id"}; only_a = set(a) - set(b); only_b = set(b) - set(a); diffs = [k for k in a if k in b and a[k] != b[k]]; assert not only_a and not only_b and not diffs, "quant.expr rows differ: only_a=%r only_b=%r differing=%d" % (only_a, only_b, len(diffs)); print("quant.expr identical: %d rows, 0 unique keys, 0 differing values" % len(a))'
	python3 -c 'import gzip; names_of = lambda p: {ln.rstrip("\n").split("\t")[5] for ln in gzip.open(p, "rt") if not ln.startswith("#") and not ln.startswith("gene_id")}; a = names_of("__sc_unchunked.LRAA.quant-only.quant.tracking.gz"); b = names_of("__sc_chunked.LRAA.quant-only.quant.tracking.gz"); assert a == b, "tracking read-name sets differ: only_unchunked=%d only_chunked=%d" % (len(a - b), len(b - a)); tagged = sum(1 for n in a if "^" in n); assert tagged == len(a), "%d of %d read names lack the barcode^umi^ encoding -- the tags did not reach the worker" % (len(a) - tagged, len(a)); print("tracking read names identical: %d rows, all %d carry the single-cell encoding" % (len(a), tagged))'


test_wdls:


clean:
	rm -rf ./chr19.bam.SS* ./LRAA* ./chr19.bam.norm*
	rm -rf ./test* ./__chckpts ./__alignment_cache
	rm -rf ./__*  ./*.count


