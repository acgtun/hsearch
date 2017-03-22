ifndef ACGTUN_UTIL
ACGTUN_UTIL=$(BLOSUM-METRIC)/src/acgtun_util
endif

app_subdirs=BLOSUM-metric

all:
	@for i in $(app_subdirs); do \
		make -C $${i} ACGTUN_UTIL=$(ACGTUN_UTIL) SRC_ROOT=$(BLOSUM-METRIC) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} ACGTUN_UTIL=$(ACGTUN_UTIL) SRC_ROOT=$(BLOSUM-METRIC) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} ACGTUN_UTIL=$(ACGTUN_UTIL) SRC_ROOT=$(BLOSUM-METRIC) test; \
	done;

clean:
	@make -C acgtun_util OPT=1 clean
	@for i in $(app_subdirs); do \
		make -C $${i} ACGTUN_UTIL=$(ACGTUN_UTIL) SRC_ROOT=$(BLOSUM-METRIC) clean; \
	done;
