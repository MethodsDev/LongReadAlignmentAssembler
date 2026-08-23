
# The equivalence gate for strandless chunking, on the bundled minigenome.
#
# Strandless chunking moves the orientation split out of a whole-BAM serial phase
# and into the per-chunk parallel phase.  That is a performance change only if it
# computes the same answer, so this directory answers whether it does.
#
# The corpus is ../single_contig's: 50 gene loci concatenated into one 3.4 Mb
# pseudo-contig.  It has no overlapping gene neighbours and therefore no read
# shared between two genes, so it exercises the plumbing and the severed-read
# accounting and says nothing about read-sharing components -- exactly as
# docs/chunked_quant_parity_evaluation.md records for the chunked-vs-whole
# comparison it is the counterpart of.  3.4 Mb yields no cuts at the default
# spacing, so the cut size is forced down here as it is there.
#
# The real-chromosome form is this same `gate` invocation pointed at a real
# chromosome and is deliberately not routine.

GATE = python3 ./run_strandless_gate.py
INPUTS = --bam ../single_contig/pacbio.PBLR.bam \
         --genome_fa ../single_contig/ref_genome.fa \
         --gtf ../single_contig/ref_annot.gtf
GEOMETRY = --approx_MB_per_cut 0.5 --approx_MB_per_cut_wiggle_window 0.1

.PHONY: test test_invariant test_ordering_cost test_gate clean


test: test_invariant test_ordering_cost test_gate


# THE CLAIM THE DESIGN RESTS ON, and the cheapest thing here to check: strand
# assignment is per-read local, so a read's strand must not depend on whether the
# split saw the whole contig or only the chunk holding it.  Checked in both modes.
#
# The aligner-flag mode is what the pipeline runs today and cannot fail for an
# interesting reason.  The inferred mode is where a read's strand can disagree
# with its alignment flag at all -- on this corpus 1,076 of 2,266 records have
# their flag REWRITTEN -- so it is the mode in which the invariant is falsifiable,
# and it is checked even though nothing turns it on yet.  A design whose soundness
# argument holds only for the settings currently in use is not sound.
test_invariant:
	rm -rf __invariant
	$(GATE) gate $(INPUTS) $(GEOMETRY) --contig minigenome --HiFi \
	    --output_dir __invariant --skip_arms --report __invariant.json


# The mistake the design must not make, executed on purpose so its cost is a
# number.  extract_contig_region_inputs filters on the RAW alignment flag and
# separate_bam_by_strand rewrites that flag, so extracting `minigenome+:...` from
# a raw BAM keeps every read the aligner called forward and then treats all of
# them as forward -- including those whose splice dinucleotides say otherwise,
# which no later stage revisits.  The output looks fine.  This target prints how
# many reads that would silently misassign, by name.
test_ordering_cost:
	rm -rf __ordering_cost
	$(GATE) ordering-cost $(INPUTS) --region 'minigenome+:1014001-1502500' \
	    --work_dir __ordering_cost --report __ordering_cost.json


# THE GATE. Both arms over the same substrate, then diffed.
#
# The pins are on the STRAND-FIRST arm, and they are what makes a strandless
# figure interpretable: without a control reproducing known values, "strandless
# gives 553 rows" does not say which arm moved.  555/544/2266/6 is the healthy
# strand-first run recorded in docs/chunked_quant_parity_evaluation.md.
#
# tracking_rows was re-pinned 541 -> 544 (and the figure above 5 -> 6, which --expect
# had already carried since 669459f4) after bisecting the move to c940609b, "Compute
# percent identity over the columns the mismatch tag counts".  That commit repaired an
# identity calculation which counted gaps in the numerator but not the denominator,
# leaving the value unbounded -- it measured -13900% -- and so discarded alignments at
# true identities of 99.0-99.9%.  This corpus is PacBio and this gate runs --HiFi, i.e.
# a high identity floor, which is exactly where those reads were being thrown away: one
# more read is now assigned (408 -> 409 distinct), carrying three more rows.  The +3 is
# this corpus's measured outcome rather than a general rule -- admitting a read changes
# EM and can change the survivor set, so row counts need not move monotonically.
#
# What did NOT move is the property this gate exists for: both arms shifted together,
# rows identical in order, every numeric column delta 0.0, severed symmetric difference
# 0.  Only the absolute control was stale, and it had been since c940609b, which is to
# say since before `make test` was last run against this target.
#
# The two arms select their cuts from different evidence -- one strand's coverage
# against both strands' union -- so their cut POSITIONS may differ and different
# positions sever different reads.  The gate therefore does not demand byte
# identity unconditionally: it demands that every differing row be attributable to
# a read in the symmetric difference of the two severed sets, and names any row
# that is not. On this corpus at 0.5 Mb the two arms' cut positions happen to
# coincide (symmetric difference 0) -- verified, not assumed; a future corpus or
# geometry change that makes them diverge is exactly what the attribution check
# above exists to keep interpretable.
test_gate:
	rm -rf __gate
	$(GATE) gate $(INPUTS) $(GEOMETRY) --contig minigenome --HiFi \
	    --cpu_budget 4 --output_dir __gate --report __gate.json \
	    --ordering_cost_region 'minigenome+:1014001-1502500' \
	    --expect expr_rows=555,tracking_rows=544,num_total_reads=2266,baseline_excluded_severed_reads=6


clean:
	rm -rf ./__*
