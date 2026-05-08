.PHONY: test build_input_bam run_quant run_reassign validate clean

test: clean build_input_bam run_quant run_reassign validate

build_input_bam:
	samtools faidx data/genome.fa
	samtools view -bS data/multigene_secondary.sam | samtools sort -o data/multigene_secondary.bam
	samtools index data/multigene_secondary.bam

run_quant:
	../../LRAA \
		--genome data/genome.fa \
		--bam data/multigene_secondary.bam \
		--gtf data/annot.gtf \
		--quant_only \
		--allow_secondary_alignments \
		--no_norm \
		--no_parallelize_contigs \
		--num_threads_per_worker 1 \
		--output_prefix LRAA \
		--verbose

run_reassign:
	python ../../util/reassign_multigene_tracking_reads.py \
		--quant_expr LRAA.quant.expr \
		--tracking LRAA.quant.tracking \
		--output_expr LRAA.cross_gene_em.quant.expr \
		--output_tracking LRAA.cross_gene_em.quant.tracking.gz

validate:
	python validate_outputs.py

clean:
	rm -f data/genome.fa.fai data/multigene_secondary.bam data/multigene_secondary.bam.bai
	rm -rf ./__*
	rm -f ./LRAA.quant.expr ./LRAA.quant.tracking ./LRAA.quant.tracking.gz
	rm -f ./LRAA.cross_gene_em.quant.expr ./LRAA.cross_gene_em.quant.tracking.gz
	rm -f ./LRAA.resources.summary.tsv ./LRAA.resources.tsv ./multigene_secondary.bam.count
