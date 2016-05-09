#
#    Copyright (C) 2015 University of Southern California
#
#    Authors: Haifeng Chen and Ting Chen
#
#    This file is part of the PCLUSTER.
#
#    PCLUSTER is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PCLUSTER is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PCLUSTER.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(PCLUSTER)/src/smithlab_cpp
endif

app_subdirs=pcluster rapsearch2

all:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(PCLUSTER) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(PCLUSTER) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(PCLUSTER) test; \
	done;

clean:
	@make -C smithlab_cpp OPT=1 clean
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(PCLUSTER) clean; \
	done;
