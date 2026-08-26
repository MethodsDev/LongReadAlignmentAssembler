# Validating a WDL change

Run **both** validators on every WDL you touch. They do not overlap, and each has
caught a real production failure the other missed.

```bash
cd WDL
miniwdl check LRAA.wdl
java -jar ~/BIN/womtool-91.jar validate LRAA.wdl
```

Then, from the repo root:

```bash
python3 -m pytest pylib/test_wdl_runtime_completeness.py -q
```

## Why two validators, and a test as well

### `miniwdl check` — what the local test harness accepts

`testing/**` runs on miniwdl, so this is the validator whose opinion the local
suite reflects. It is fast and it catches type errors, unbound identifiers and
bad call inputs.

What it does NOT do: **miniwdl does not require a container image.** A task with
no `docker` in its `runtime` block gets a default image and runs happily. So a
missing image is invisible to `miniwdl check`, to `make test_wdls`, and to every
smoke run in this repo.

### `womtool validate` — what Terra actually runs

Terra/Cromwell is a different engine, and constructs miniwdl accepts are not
always portable. Use womtool to prove a change is engine-portable before it
reaches a Terra submission: `File x = select_first([...])` inside a conditional,
cross-file task calls, nested-subworkflow input plumbing.

`womtool-91.jar` matches the Cromwell generation Terra runs. It lives at
`~/BIN/womtool-91.jar` on the dev workstation; download the matching release if
it is absent.

### The runtime-completeness test — what NEITHER validator catches

An absent runtime attribute is legal WDL, so neither validator objects, but GCP
Batch refuses the job:

```
No container image found in either 'container' or 'docker' runtime attributes.
```

This is not hypothetical. `require_annot_gtf` in `LRAA-cell_cluster_guided.wdl`
shipped with `cpu`, `memory` and `disks` and no `docker`. It passed
`miniwdl check`, passed a full `make test_wdls`, passed three complete `by_chunk`
single-cell smoke runs, and then failed a real Terra submission. It was reachable
only in `quant_only` mode, so nothing local that ran it was a Cromwell run.
Diagnosis cost several wrong hypotheses -- a stale image tag, an empty Terra
input, a mixed method snapshot -- before anyone asked the code the mechanical
question the error message already named.

`pylib/test_wdl_runtime_completeness.py` now asserts that every task in every WDL
under `WDL/` names `docker` or `container`, parametrized per file so a failure
names the file and the task. It carries a guard-the-guard case, because its
parsing could silently stop matching and every file would pass vacuously -- which
is how the original defect survived a green suite.

## Checklist before pushing a WDL change

1. `miniwdl check` on every file you touched **and** every file that imports it
   (imports are relative, so a change to `subwdls/X.wdl` affects its callers).
2. `womtool validate` on the same set.
3. `pytest pylib/test_wdl_runtime_completeness.py`.
4. If the change alters a task's runtime block, say in the commit message which
   backend you verified against. Local green does not imply Terra green.

## Note on images and Terra

`Docker/build_docker.testing.sh` refuses to build from an unpushed commit,
because images are labelled with the SHA. `make test_wdls` additionally refuses
to run when the image revision and the worktree disagree
(`check_test_image_revision`) -- that guard is doing its job when it fires, and
committing while a suite runs is what trips it.

For a Terra submission, bind the commit-pinned tag rather than `:testing`:

```
us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:0.28.0-<short-sha>
us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:0.28.0-<short-sha>
```

`:testing` moves on every rebuild, so a run bound to it is not reproducible.
