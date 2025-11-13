DIRS = sirvs quant_eval single_contig sep_contigs ont_sep_contigs read_to_isoform_class_assignments single_cells

QC = read_to_isoform_class_assignments

MISC = misc_tests

test:
	@for i in $(DIRS) $(MISC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

test_QC:
	@for i in $(QC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

test_docker:
	cd sirvs && ${MAKE} test_docker


test_wdls:
	cd sep_contigs && $(MAKE) test_wdls
	cd single_contig && $(MAKE) test_wdls
	cd single_cells && $(MAKE) test_wdls


clean:
	@for i in $(DIRS) $(QC) $(MISC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done
