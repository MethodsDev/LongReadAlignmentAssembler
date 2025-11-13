DIRS = sc_most_basic  sc_add_gene_symbols sc_cluster_guided_pipe  sc_full_pipe


test:
	@for i in $(DIRS); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

test_wdls:
	@for i in $(DIRS); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test_wdls) || exit $$?; done

clean:
	@for i in $(DIRS); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done
