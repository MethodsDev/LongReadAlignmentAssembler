# Notes for agents working in this repo

## Docker images and releases

- **"Rebuild the dockers" means testing images**: `cd Docker && bash build_docker.testing.sh`.
  It writes `testing`, `<version>-testing` and `<version>-<shortsha>`. This holds after a
  version bump too: a version bump is not a release.
- Only when the user explicitly says this is an **official / public release**:
  fast-forward `main`, push an `LRAA_v<version>` tag, publish the GitHub release, then run
  `Docker/release_docker.OFFICIAL.sh` (bare `<version>` and `:latest`). Confirm each of
  those outward steps with the user first.
- If a release script refuses (HEAD is not `origin/main`, or there is no published GitHub
  release for HEAD), **stop and ask**. Never move `main`, create a tag or publish a GitHub
  release to get past the check; that check is the release decision.
- Never hand-retag `:latest` or bare version tags. See `Docker/README.md`, "Tags".

## Branches

- Work happens on `devel`. `main` holds only the current official release, and moves
  only by fast-forward at release time.
- Several agents may work in this checkout at once. Commit only the files you changed.
