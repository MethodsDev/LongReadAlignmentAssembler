# Execution test for the library-level TSS cap calibration (v0.47.0): the standalone
# step decides from saved scan tables (no bam) and writes its report and plot; runs are
# kept under __test_runs/ (summary.md shows every case's histogram). See README.md.

.PHONY: test clean

test:
	python3 test_TSS_cap_calibration.py

clean:
	rm -rf __test_runs
