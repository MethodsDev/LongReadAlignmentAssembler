DIRS = sirvs quant_eval single_contig sep_contigs ont_sep_contigs read_to_isoform_class_assignments

QC = read_to_isoform_class_assignments

test:
	@for i in $(DIRS); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

test_QC:
	@for i in $(QC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done


test_docker:
	cd sirvs && ${MAKE} test_docker

clean:
	@for i in $(DIRS) $(QC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done
