ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef ACGTUN_UTIL
ACGTUN_UTIL=./../acgtun_util
endif

PROGS = evaluate

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50

TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

SRC = k_nearest_neighbor_similarity.cpp k_nearest_neighbor_distance.cpp
OBJ = $(patsubst %.cpp,%.o,$(SRC))

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<  -I$(ACGTUN_UTIL)

evaluate : $(OBJ) \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

%: %.cpp
		$(CXX) $(CFLAGS) -o $@ $^  -I$(ACGTUN_UTIL)
	
install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
