BLOSUM-METRIC = $(shell pwd)

all:
	@make -C src BLOSUM-METRIC=$(BLOSUM-METRIC) OPT=1

install:
	@make -C src BLOSUM-METRIC=$(BLOSUM-METRIC) OPT=1 install

test:
	@make -C src BLOSUM-METRIC=$(BLOSUM-METRIC) test
.PHONY: test

clean:
	@make -C src BLOSUM-METRIC=$(BLOSUM-METRIC) clean
.PHONY: clean

distclean: clean
	@rm -rf $(BLOSUM-METRIC)/bin
	@rm -rf $(BLOSUM-METRIC)/include
.PHONY: distclean
