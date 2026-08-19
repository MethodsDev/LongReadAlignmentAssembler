DIRS = sc_most_basic sc_basic_quant_only  sc_add_gene_symbols sc_cluster_guided_pipe  sc_full_pipe sc_full_pipe_scattered sc_monoexonic_cell_filter_deferral


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
