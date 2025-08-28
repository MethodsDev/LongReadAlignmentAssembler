DIRS = internal_priming_removal

test:
	@for i in $(DIRS); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

clean:
	@for i in $(DIRS) $(QC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done
