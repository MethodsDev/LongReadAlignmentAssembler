# Execution regression test for the whole-genome alignment-mismapping filter
# (v0.40.0). Runs the checked-out ../../LRAA on de-identified minigenomes, with
# the filter off vs on, and asserts the artifacts are removed and real models
# kept. See README.md. SLOW target (runs LRAA 8x); not part of the fast pytest.

.PHONY: test clean

test:
	python3 test_mismapping_filter.py

clean:
	rm -rf __test_runs
