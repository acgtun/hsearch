#
#    Copyright (C) 2016 University of Southern California
#
#    Authors: Haifeng Chen and Ting Chen
#
#    This file is part of the HCLUST.
#
#    HCLUST is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HCLUST is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HCLUST.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(HCLUST)/src/smithlab_cpp
endif

app_subdirs=hclust

all:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(HCLUST) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(HCLUST) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(HCLUST) test; \
	done;

clean:
	@make -C smithlab_cpp OPT=1 clean
	@for i in $(app_subdirs); do \
		make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) SRC_ROOT=$(HCLUST) clean; \
	done;
