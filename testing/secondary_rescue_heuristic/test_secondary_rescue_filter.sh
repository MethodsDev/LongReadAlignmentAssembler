#!/usr/bin/env bash
set -euo pipefail

rm -rf work input.bam input.bam.bai filtered.bam filtered.bam.bai filtered.summary.tsv
samtools view -bS input.sam | samtools sort -o input.bam
samtools index input.bam

../../util/filter_bam_to_secondary_rescue.py \
  --input_bam input.bam \
  --output_bam filtered.bam \
  --summary_tsv filtered.summary.tsv \
  --threads 1 \
  --workdir work

samtools view filtered.bam | cut -f 1,2 > observed.records.tsv
cat > expected.records.tsv <<'EOF'
read_branch1	0
read_branch1	256
read_branch2	0
read_branch2	256
read_drop_nm	0
read_missing_ms	0
read_supp	0
EOF
diff -u expected.records.tsv observed.records.tsv

python - <<'PY'
import csv
with open("filtered.summary.tsv") as fh:
    row = next(csv.DictReader(fh, delimiter="\t"))
assert row["kept_non_secondary"] == "5", row
assert row["secondary_candidates"] == "4", row
assert row["kept_secondary_heuristic"] == "2", row
assert row["kept_secondary_branch1"] == "1", row
assert row["kept_secondary_branch2"] == "1", row
assert row["dropped_secondary_failed_delta_nm_rate"] == "1", row
assert row["dropped_secondary_missing_metrics"] == "1", row
assert row["dropped_supplementary"] == "1", row
PY
