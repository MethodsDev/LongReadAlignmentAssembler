#!/bin/bash
#
# Retired.  The bare-<version> and :latest builds are now ONE script that also
# requires a published GitHub release of HEAD: release_docker.OFFICIAL.sh.
# Kept as a stub so muscle memory lands on the explanation, not on a build.

echo "" >&2
echo "$(basename $0) is retired." >&2
echo "" >&2
echo "  Not an official public release (including a version bump):" >&2
echo "      bash build_docker.testing.sh       # testing, <version>-testing, <version>-<shortsha>" >&2
echo "" >&2
echo "  Official public release, after its GitHub release is published:" >&2
echo "      bash release_docker.OFFICIAL.sh    # <version> and :latest" >&2
echo "" >&2
exit 1
