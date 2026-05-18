.PHONY: test clean $(DIRS)

DIRS = AT3G28550 AT4G05320

test:
	@for dir in $(DIRS); do \
		$(MAKE) -C $$dir test || exit $$?; \
	done

clean:
	@for dir in $(DIRS); do \
		$(MAKE) -C $$dir clean || exit $$?; \
	done
