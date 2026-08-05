#!/usr/bin/env python3

import pytest

from Transcript import GTF_contig_to_transcripts


@pytest.mark.parametrize(
    ("source_gene_id", "gene_name", "expected_output_gene_id"),
    [
        ("ENSG00000198888", "MT-ND1", "MT-ND1^ENSG00000198888"),
        ("PABPC4", "PABPC4", "PABPC4"),
        ("MT-ND1^ENSG00000198888", "MT-ND1", "MT-ND1^ENSG00000198888"),
        ("OLD^ENSG00000198888", "NEW", "OLD^ENSG00000198888"),
        ("ENSG00000198888", None, "ENSG00000198888"),
    ],
)
def test_gtf_gene_name_is_composed_only_for_output(
    tmp_path, source_gene_id, gene_name, expected_output_gene_id
):
    gtf_path = tmp_path / "reference.gtf"
    gene_name_attribute = f' gene_name "{gene_name}";' if gene_name else ""
    attributes = (
        f'gene_id "{source_gene_id}"; transcript_id "ENST00000361390";'
        f"{gene_name_attribute}"
    )
    gtf_path.write_text(
        f"chrM\tref\ttranscript\t1\t100\t.\t+\t.\t{attributes}\n"
        f"chrM\tref\texon\t1\t100\t.\t+\t.\t{attributes}\n"
    )

    transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(str(gtf_path))
    transcript = transcripts["chrM"][0]

    assert transcript.get_gene_id() == source_gene_id
    assert transcript.get_output_gene_id() == expected_output_gene_id

    output_gtf = transcript.to_GTF_format()
    assert output_gtf.count(f'gene_id "{expected_output_gene_id}";') == 2
