.PHONY: test clean $(DIRS)

DIRS = Mpc1_ENSMUST00000155364.8 Phb2_ENSMUST00000004375.16 Psap_ENSMUST00000177779.8 Septin7_ENSMUST00000165594.4

test:
	@for dir in $(DIRS); do \
		$(MAKE) -C $$dir test || exit $$?; \
	done

clean:
	@for dir in $(DIRS); do \
		$(MAKE) -C $$dir clean || exit $$?; \
	done
