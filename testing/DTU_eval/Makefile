
teset: basic ignore_unspliced recip_delta_pi recip_top1ea recip_top1ea_ignore_unspliced 


basic:
	../../util/sc/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu


ignore_unspliced:
	../../util/sc/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.ignore_unspliced --ignore_unspliced


recip_delta_pi:
	../../util/sc/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.recip_delta_pi --reciprocal_delta_pi

recip_top1ea:
	../../util/sc/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.recip_delta_pi.top1ea --reciprocal_delta_pi --top_isoforms_each 1

recip_top1ea_ignore_unspliced:
	../../util/sc/sc_pseudobulk_test_isoform_DiffUsage.py  --sc_cluster_counts_matrix testdata/cluster_pseudobulk_read_counts.matrix_for_DTU.tsv --splice_hashcode_id_mappings testdata/gene_transcript_splicehashcode.tsv.wAnnotIDs --output_prefix test.dtu.recip_delta_pi.top1ea.ignore_unspliced --reciprocal_delta_pi --top_isoforms_each 1 --ignore_unspliced


clean:
	rm -f ./*.diff_iso.tsv

