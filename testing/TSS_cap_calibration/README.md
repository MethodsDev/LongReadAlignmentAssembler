# TSS cap calibration: execution test

Exercises the decision step of `pylib/TSS_cap_calibration.py` -- whether a library carries
the untemplated-G cap signature, and with which baseline G run -- through the standalone
command a workflow runs before LRAA, `util/calibrate_TSS_cap_signature.py`. No bam is needed:
the command decides from a saved **scan table** (`--scan_table`), the per-peak and off-peak
5' G-run counts that every bam scan writes (`<prefix>.TSS_cap_calibration.scan.tsv.gz`, from
LRAA's own pre-pass or the standalone command).

Each case checks the decision, the config keys handed to LRAA (`.config.json`, for
`--config_update`) and the report files: `.TSS_cap_calibration.json` (decision and
statistics), `.hist.tsv` (per-peak capped-fraction histogram) and `.svg` (its plot).
`__test_runs/summary.md` tabulates the cases and shows each plot. `make test` takes a few
seconds.

| case | scan table | expected |
|---|---|---|
| `BT474_on` | bulk PacBio BT474 | on, baseline 1 G |
| `SIRV_undetermined` | SIRV spike-ins | undetermined (38 peaks), no config keys |
| `sim_off` | Arabidopsis simulation | off: one mode at 0 |
| `BT474_floor` | BT474, peak floor raised by `--config_update` | undetermined |
| `BT474_merge` | BT474 with `--merge` | caller's keys kept with the calibration's |

## Data

Scans of the LRAA paper's benchmark bams (2026-10-29), saved with
`write_scan_table(..., min_peak_reads=20)`: every contig strand's off-peak row, and the peaks
with at least 20 informative spliced reads, the default floor for judging (so the decisions
are those of the full scans):

- `BT474_PacBio.scan.tsv.gz`: BT474 E0 PacBio bulk, whole genome (11,922 peaks).
- `SIRV.scan.tsv.gz`: the SIRV spike-in reads of the same BT474 E0 library (38 peaks).
- `Arabidopsis_sim.Chr4_6M.scan.tsv.gz`: Arabidopsis isoseqsim (error 0.16) reads, Chr4
  1-6,000,000 (469 peaks).

The per-library study behind the decision rule, with ten libraries, is in the paper analysis
repo, `__TSS_5prime_Completeness/cap_calibration/`.
