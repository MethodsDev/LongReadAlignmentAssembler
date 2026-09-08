# Streaming quant-only canonical-path collision (LRAA v0.31.0).  See README.md.
#
# test        the pytest leg, against the LRAA in THIS checkout.  LRAA_HOME is
#             pinned on the recipe rather than left to the parent walk in
#             _lraa_cmd(), because _lraa_cmd() honours an inherited LRAA_HOME
#             FIRST: without the pin, a caller with LRAA_HOME set to another
#             install would have `make test` silently exercise that install
#             instead of the working tree.  ~3 s.
# test_image  the container leg, which needs an image:
#                 make test_image IMG=/path/to/lraa-core_<tag>.sif
#
# python3 -m pytest rather than bare pytest, for the sys.path reason spelled out
# in testing/Makefile.

.PHONY: test test_image clean

test:
	LRAA_HOME=$(abspath $(CURDIR)/../../..) python3 -m pytest -q test_streaming_quant_canonical_path_collision.py

test_image:
	@[ -n "$(IMG)" ] || { echo "usage: make test_image IMG=<lraa-core sif or docker uri>" >&2; exit 1; }
	./run_test.sh $(IMG)

# Both legs run entirely inside throwaway temp directories (pytest's
# tmp_path_factory and run_test.sh's mktemp -d), so nothing lands here.
clean:
	@:
