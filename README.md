# Long Read Alignment Assembler (LRAA)

Isoform Discovery and/or Quantification from Long Read RNA-Seq

Visit the [LRAA wiki](https://github.com/MethodsDev/LongReadAlignmentAssembler/wiki) for user documentation


## Resource usage monitoring (optional)

Lightweight CPU and memory sampling is enabled by default to help size runs and tune parallelism:

- Disable monitoring: pass `--no_monitor_resources`
- Sampling interval: `--monitor_interval <seconds>` (default 60.0)
- Include child processes: `--monitor_children` (on by default)

Outputs:
- Main process: `<output_prefix>.resources.tsv`
- In contig-parallel runs (default), each worker also writes: `__<output_prefix>.contigtmp/<contig>/<strand>/<contig>.<strand>.resources.tsv`

Columns (TSV): `epoch_ts, elapsed_sec, rss_mb, cpu_percent, rss_mb_children, cpu_percent_children, note`

Notes:
- Monitoring uses `psutil`. In the provided Docker image it is installed. On bare-metal installs, if `psutil` is unavailable, monitoring will auto-disable with a warning.

## Progress updates during mapping and quantification

LRAA emits progress while:
1) Mapping read alignments to the splice graph ("map-reads")
2) Assigning reads to assembled isoforms ("quant-assign")

You can control this via config overrides:

- Enable/disable: `show_progress_quant_assign` (default: true)
- Enable/disable (mapping): `show_progress_mapping` (default: true)
- Update every N records: `progress_update_every_n` (default: 1000)
- Update at least every S seconds: `progress_update_interval_sec` (default: 5.0)
- Mapping stage update frequency: `mapping_update_every_n` (default: 10000), `mapping_update_interval_sec` (default: 2.0)
- Prefer tqdm progress bar if available: `use_tqdm_progress` (default: true)

Example using `--config_update`:

```
./LRAA \
	--bam sample.bam \
	--genome genome.fa \
	--gtf targets.gtf \
		--config_update '{"show_progress_mapping": true, "show_progress_quant_assign": true, "use_tqdm_progress": true, "mapping_update_every_n": 5000, "progress_update_every_n": 500, "progress_update_interval_sec": 2.0}'
```

Notes:
- If `tqdm` is installed, LRAA will show a dynamic progress bar by default. If `tqdm` is not available, LRAA falls back to a lightweight stderr progress line.
- The provided Docker image now includes `tqdm`. On bare-metal installs, you can add it via:
	- `python3 -m pip install tqdm`

## Quantification-only for single-cell clusters (shared splice graph)

Quantify multiple cell clusters separately while building a single shared splice graph from all cluster BAMs by supplying a BAM list file via `--bam_list` in quant-only mode:

```
./LRAA \
	--bam_list clusters.bams.txt \
	--genome genome.fa \
	--gtf targets.gtf \
	--quant_only \
	--output_prefix LRAA
```

The `--bam_list` file should contain one entry per line, either:

- `<cluster_id>\t/path/to/cluster.bam` (tab or whitespace separated), or
- `/path/to/cluster.bam` (cluster ID is the BAM basename without `.bam`)

Behavior:

- Builds one splice graph per contig/strand using the union of all listed BAMs together with the provided GTF.
- Quantifies each cluster BAM independently against that shared graph and writes per-cluster outputs:
	- `LRAA.<cluster_id>.quant.expr`
	- `LRAA.<cluster_id>.quant.tracking`
- With `--tag_bam`, the respective cluster BAM is annotated using its per-cluster tracking file.

Notes:

- `--bam_list` is only supported with `--quant_only`.
- `--CPU` is respected; contig-parallel execution is the default. Use `--no_parallelize_contigs` to disable contig-level parallelism and use inner component-level multithreading instead; resources are managed per cluster accordingly.
- `TPM` normalization uses the read count of each cluster BAM independently. `--num_total_reads` is not allowed with `--bam_list` and will raise an error.

