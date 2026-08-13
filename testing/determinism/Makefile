.PHONY: test smoke strong clean

# Smoke test on bundled SIRV data.  Fast, and its reach is measured rather than
# assumed: run against pristine 0301c30 it FAILS on 4 of 6 principal outputs
# with "bytes differ, content identical", and passes on the fixed tree.  So it
# is a verified guard against the merge-ordering regression.
#
# It is NOT a guard against the content regression: on this data the content
# was identical even on the unfixed tree.  A minigenome has no multipath
# component at or above min_mpgn_component_size_for_spawn (150), so it never
# engages the forked component path and its nondeterminism surface is strictly
# smaller than a real contig's.  Use the strong target for that.
test: smoke

smoke:
	python3 check_determinism.py \
		--genome ../sirvs/data/SIRVs1-7.genome.fa \
		--bam ../sirvs/data/sim.fasta.mm2.bam \
		--gtf ../sirvs/data/SIRVs1-7.annot.gtf \
		--replicates 3 \
		--threads 2 \
		-- --HiFi

# Strong guard.  Requires a region known to carry genuine content churn, not
# merely permutation churn: chr20:30-34Mb differs only by reordering between
# runs and would pass a content check while nondeterministic code sat
# underneath, whereas chr20:36,000,000-37,000,000 carries real content
# differences.  Budget for cold runs: each replicate uses its own output
# directory and therefore builds its own preprocessing cache, so the ~2 min
# figure quoted for that region elsewhere assumed a warm cache and does not
# apply here.  Expect roughly a full chr20 preprocessing pass per replicate.
#
#   make strong GENOME=/path/chr20.fa BAM=/path/chr20.bam GTF=/path/guide.gtf
strong:
	@test -n "$(GENOME)" -a -n "$(BAM)" || \
		{ echo "set GENOME= and BAM= (and optionally GTF=)"; exit 2; }
	python3 check_determinism.py \
		--genome $(GENOME) \
		--bam $(BAM) \
		$(if $(GTF),--gtf $(GTF),) \
		--region $(if $(REGION),$(REGION),chr20:36000000-37000000) \
		--replicates 3 \
		--threads 4 \
		-- --HiFi

clean:
	rm -rf det_work
