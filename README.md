# Long Read Alignment Assembler (LRAA)

Isoform Discovery and/or Quantification from Long Read RNA-Seq

Visit the [LRAA wiki](https://github.com/MethodsDev/LongReadAlignmentAssembler/wiki) for user documentation


## Resource usage monitoring (optional)

You can enable lightweight CPU and memory sampling to help size runs and tune parallelism:

- Enable monitoring: pass `--monitor_resources`
- Sampling interval: `--monitor_interval <seconds>` (default 2.0)
- Include child processes: `--monitor_children` (on by default)

Outputs:
- Main process: `<output_prefix>.resources.tsv`
- When `--parallelize_contigs` is used, each worker also writes: `<output_prefix>.contigtmp/<contig>/<strand>/<contig>.<strand>.resources.tsv`

Columns (TSV): `epoch_ts, elapsed_sec, rss_mb, cpu_percent, rss_mb_children, cpu_percent_children, note`

Notes:
- Monitoring uses `psutil`. In the provided Docker image it is installed. On bare-metal installs, if `psutil` is unavailable, monitoring will auto-disable with a warning.

