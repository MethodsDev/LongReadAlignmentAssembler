.PHONY: test clean $(DIRS)

DIRS = microexon repetitive

test:
	@for dir in $(DIRS); do \
		$(MAKE) -C $$dir test || exit $$?; \
	done

clean:
	@for dir in $(DIRS); do \
		$(MAKE) -C $$dir clean || exit $$?; \
	done
