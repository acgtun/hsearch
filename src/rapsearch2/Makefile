# Makefile
#

ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=./../smithlab_cpp
endif

objects := HashSearch.o BlastStat.o Seg.o mergeUnit.o read_proteins.o lsh.o

CXX := g++
CXXFLAGs :=

ifeq ($(dbg), 1)
	CXXFLAGS += -g -Wall -std=c++11 
else
	CXXFLAGS += -O3 -w -std=c++11
endif

ifeq ($(pg), 1)
	CXXFLAGS += -pg
else
	CXXFLAGS +=
endif


INC := -I ./
LIBS := -L ./

all: rapsearch prerapsearch pcluster

rapsearch: main.cpp $(objects)
	$(CXX) $(CXXFLAGS) -o rapsearch main.cpp $(objects) $(INC) $(LIBS) -lboost_serialization -lpthread -lboost_thread -lboost_system -lboost_chrono

prerapsearch: BuildHash.cpp $(objects)
	$(CXX) $(CXXFLAGS) -o prerapsearch BuildHash.cpp $(objects) $(INC) $(LIBS) -lboost_serialization -lpthread -lboost_thread -lboost_system -lboost_chrono

pcluster: pcluster.cpp $(objects)
	$(CXX) $(CXXFLAGS) -o pcluster pcluster.cpp $(objects) $(INC) $(LIBS) $(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o OptionParser.o) \
	-lboost_serialization -lpthread -lboost_thread -lboost_system -lboost_chrono

$(objects) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)

install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

clean:
	$(RM) *.o rapsearch prerapsearch pcluster
