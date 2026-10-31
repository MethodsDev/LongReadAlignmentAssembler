# Releasing LRAA

How versions, development builds and official releases work, and the exact steps for an
official release. The Docker tag policy is in `Docker/README.md` ("Tags"); the short rules for
agents are in `CLAUDE.md`. This document is the full procedure.

## What a release is

An **official release** `vX.Y.Z` is three things that always name the same commit:

1. `main` fast-forwarded to that commit (`main` holds only the current official release);
2. an annotated tag `LRAA_vX.Y.Z` and a **published** GitHub release of that tag;
3. the Docker images `lraa`, `lraa-core`, `lraa-sc`, `lraa-orf` at the bare tag `X.Y.Z`, and
   `:latest` moved to them, written by `Docker/release_docker.OFFICIAL.sh`.

The GitHub release is the decision; the images follow it. Everything else is a development
build:

| build | tags | written by | meaning |
|---|---|---|---|
| official release | `X.Y.Z`, `latest` | `Docker/release_docker.OFFICIAL.sh` | a published release |
| development | `X.Y.Z-<shortsha>`, `X.Y.Z-testing`, `testing` | `Docker/build_docker.testing.sh` | a commit someone built; pin `X.Y.Z-<shortsha>` for runs you need to reattribute later |

**A version bump is not a release.** Bumping `VERSION` and building images marks a candidate;
nothing becomes official until the GitHub release is published.

## Version numbers

`vMAJOR.MINOR.PATCH` (pre-1.0, MAJOR stays 0):

- **MINOR** (0.47.0 -> 0.48.0): the change alters the outputs of a default run (models, sites,
  quantification, flags), or changes a workflow's inputs/outputs in a way callers notice.
- **PATCH** (0.47.0 -> 0.47.1): no change to the outputs of a default run. Runtime, cost,
  memory, logging, docs, new options that are off by default, and fixes to outputs that were
  plainly wrong in a narrow case (state it in the changelog).

The version is set in two places, changed together in one commit titled
`Version vX.Y.Z: <summary>`:

- `LRAA` (`VERSION = "vX.Y.Z"`)
- `Docker/VERSION.txt` (`X.Y.Z`)

## Changelog

`CHANGELOG.txt`, newest first, one section per version. Each section says up front whether run
outputs change ("Changes that alter run outputs (...)" / "No change to run outputs" / "Also in
this version, no change to default runs"), and gives MEASURED evidence for output changes: the
dataset, the region, and the numbers before and after. The release notes are a short,
user-facing summary of the changelog sections since the previous release.

## Before releasing

- The release commit is on `devel` and pushed. It need not be `devel`'s head: release the
  commit that was validated (e.g. the commit whose images produced an analysis), and let later
  devel commits go in the next release.
- `pytest pylib` passes (including the serial tests), and the execution tests under `testing/`
  pass on that commit.
- For a MINOR release: the changelog's MEASURED evidence exists for the output changes. For
  the LRAA paper this means the benchmark families and the single-cell evaluations were rerun
  on the candidate's development images (`X.Y.Z-<shortsha>`).
- For a PATCH release: outputs shown identical to the previous release on at least one sample
  per mode (bulk ref-free, ref-guided, quant-only, and the single-cell test pipeline).
- No Terra submission that tracks the WDL on Dockstore `devel` still needs a retry with the
  old images: Terra configs follow the branch's WDL, so a retry after the WDL changed would
  pair new WDL with old images.

## Steps (official release)

Only when the maintainer explicitly asks for an official / public release, and confirming
each outward step (1-4) with them first.

```bash
REL=<commit>                         # the validated release commit, already on origin/devel
V=X.Y.Z                              # must equal Docker/VERSION.txt at $REL

# 1. fast-forward main (a plain push refuses anything but a fast-forward)
git push origin ${REL}:refs/heads/main

# 2. annotated tag
git tag -a LRAA_v${V} ${REL} -m "LRAA v${V}"
git push origin LRAA_v${V}

# 3. GitHub release (notes: see "Release notes" below)
gh release create LRAA_v${V} --repo MethodsDev/LongReadAlignmentAssembler \
    --title "LRAA Release v${V}" --notes-file release_notes.md --latest

# 4. official images, from a checkout whose HEAD is the release commit
git worktree add --detach /tmp/LRAA_release_v${V} ${REL}
cd /tmp/LRAA_release_v${V}/Docker && bash release_docker.OFFICIAL.sh
#    then remove the worktree: git worktree remove /tmp/LRAA_release_v${V}
```

`release_docker.OFFICIAL.sh` refuses unless HEAD equals `origin/main` and `LRAA_v<VERSION>` is
a published GitHub release whose tag resolves to HEAD. **If it refuses, stop and ask.** Never
move `main`, create a tag or publish a release to get past it, and never hand-retag `:latest`
or a bare version (a `docker tag` + `docker push` bypasses the guard).

After the build, check the registry: `0.X.Y` and `latest` of `lraa-core`, `lraa-sc`, `lraa-orf`
and `lraa` all name the new digest.

## Release notes

Same layout every time (see `LRAA_v0.44.3`, `LRAA_v0.47.0`):

```
LRAA vX.Y.Z, <Month D, YYYY>

`CHANGELOG.txt` has the full detail (vA.B.C through vX.Y.Z). Docker:
`us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:X.Y.Z` (also `:latest`; `lraa-core`,
`lraa-sc`, `lraa-orf` at the same tag).

**Output changes**      (what a rerun of a default run will show differently, with numbers)
**Workflows**           (WDL inputs/outputs/tasks)
**New options and tools** (no change to default runs)
**Fixes and performance**
```

Write for users of the tool: what changed in their results and why, not the implementation.
No internal paths or commit hashes in the notes.

## After releasing

- `devel` continues from where it was; the next version bump starts the next candidate.
- Dockstore: the `main` branch's WDL is now the release's; Terra configs that should stay on
  this release pin the image `X.Y.Z` (a bare version asserts a published release).
- Analyses that used development images of the release commit (`X.Y.Z-<shortsha>`) can cite
  the release version: the official image is built from the same commit.

## History

- `LRAA_v0.47.0` (2026-10-31): released from `09025bb2`, the commit whose development images
  (`0.47.0-09025bb`) produced the LRAA paper's v0.47.0 analyses; the two later devel commits
  (docs and comments) go into the next release. This document was written at that release.
