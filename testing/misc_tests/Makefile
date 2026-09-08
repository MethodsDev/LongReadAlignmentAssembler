# regression_test_canonical_path_collision first: it is a 3 s fixture-scale test
# of the working tree, so a regression there is reported before the multi-minute
# runs below rather than after them.
DIRS = regression_test_canonical_path_collision internal_priming_removal LRAA_to_ORFs merge_LRAA_GTFs EMP3_alt_termini

test:
	@for i in $(DIRS); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

clean:
	@for i in $(DIRS) $(QC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done
