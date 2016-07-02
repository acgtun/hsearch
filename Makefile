#
#    Copyright (C) 2016 University of Southern California
#
#    Authors: Haifeng Chen and Ting Chen
#
#    This file is part of the PMF.
#
#    PMF is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PMF is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PMF.  If not, see <http://www.gnu.org/licenses/>.
#

PMF = $(shell pwd)

all:
	@make -C src PMF=$(PMF) OPT=1

install:
	@make -C src PMF=$(PMF) OPT=1 install

test:
	@make -C src PMF=$(PMF) test
.PHONY: test

clean:
	@make -C src PMF=$(PMF) clean
.PHONY: clean

distclean: clean
	@rm -rf $(PMF)/bin
	@rm -rf $(PMF)/include
.PHONY: distclean
