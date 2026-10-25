# Long Read Alignment Assembler (LRAA)

Isoform Discovery and/or Quantification from Long Read RNA-Seq

Visit the [LRAA wiki](https://github.com/MethodsDev/LongReadAlignmentAssembler/wiki) for user documentation

## Docker images

```
us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest
```

`:latest` is the **current official public release** -- not the newest build. Use
it, or pin the release version it resolves to if you need a fixed reference. The
version it names moves only when a release is made.

Development builds are never published under `:latest`; they are tagged
`<version>-testing` or `<version>-<shortsha>`. See
[Docker/README.md](Docker/README.md) for the full tag policy and the build
scripts that enforce it.

## Fragmented and preassembled inputs

`--pasa_mode` chains fragmented transcript alignments (ESTs, partial cDNAs,
staggered reads) that share introns into longer models; `--preassembled_input`
accepts single-alignment support for preassembled transcripts (Trinity, StringTie).
Long-read RNA-seq needs neither. See [docs/pasa_mode.md](docs/pasa_mode.md).

## Developer utilities

`util/ascii_isoform_view.py` draws transcript structures as stacked ASCII rows —
reference GTF against `LRAA.gtf`, or a multipath path through the splice graph —
with per-row structural verdicts. See
[docs/ascii_isoform_illustrator.md](docs/ascii_isoform_illustrator.md).

