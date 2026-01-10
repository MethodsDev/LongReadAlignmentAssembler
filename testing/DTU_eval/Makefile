
test: basic cell_fraction_filter ignore_unspliced recip_delta_pi recip_top1ea recip_top1ea_ignore_unspliced pytest 


basic:
	../../util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu


cell_fraction_filter:
	../../util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.min_cell_frac --sc_cluster_fraction_matrix testdata/cluster_pseudobulk.cell_fractions.matrix  --min_cell_fraction 0.1

ignore_unspliced:
	../../util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.ignore_unspliced --ignore_unspliced


recip_delta_pi:
	../../util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.recip_delta_pi --reciprocal_delta_pi

recip_top1ea:
	../../util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.recip_delta_pi.top1ea --reciprocal_delta_pi --top_isoforms_each 1

recip_top1ea_ignore_unspliced:
	../../util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.recip_delta_pi.top1ea.ignore_unspliced --reciprocal_delta_pi --top_isoforms_each 1 --ignore_unspliced

pytest:
	python3 -m pytest test_diff_iso_tests.py


clean:
	rm -f ./*.diff_iso.tsv

