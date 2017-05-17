## This file is part of JSeqArray.
##
## JSeqArray is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License Version 3 as
## published by the Free Software Foundation.
##
## JSeqArray is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public
## License along with JSeqArray.
## If not, see <http://www.gnu.org/licenses/>.


## The command which runs the C compiler.
CC ?= gcc
## The command which runs the C++ compiler.
CXX ?= g++

## Debugging and optimization options for the C compiler.
CFLAGS = -O2
## Options for the C++ compiler.
CXXFLAGS = -O2

## Options for linker
LDFLAGS =
LDLIBS = -lz -llzma


## Remove directory entries
RM = rm -rf


## Determine the filename extension
platform=$(shell uname)
ifeq ($(platform),Linux)
	LIBEXT=so
	LINKLIBS=-lpthread
else ifeq ($(platform),Darwin)
	LIBEXT=dylib
	LINKLIBS=-lpthread
else
	LIBEXT=dll
	LINKLIBS=
endif


## julia program
JL_SHARE = $(shell julia -e 'print(joinpath(JULIA_HOME,Base.DATAROOTDIR,"julia"))')
CFLAGS   += $(shell $(JL_SHARE)/julia-config.jl --cflags)
CXXFLAGS += $(shell $(JL_SHARE)/julia-config.jl --cflags)
LDFLAGS  += $(shell $(JL_SHARE)/julia-config.jl --ldflags)
LDLIBS   += $(shell $(JL_SHARE)/julia-config.jl --ldlibs)


## Define lib macro and include
GDS_SHARE = $(shell julia -e 'using jugds; print(gds_get_include())')
LIBMACRO =
LIBINCLUDE = -I$(GDS_SHARE)


CFLAGS   += $(LIBMACRO) $(LIBINCLUDE)
CXXFLAGS += $(LIBMACRO) $(LIBINCLUDE)


## JSeqArray library object files
LIB_OBJS = GetData.o Index.o JSeqArray.o LinkGDS.o ReadByVariant.o vectorization.o


## all jobs
all: $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) -shared \
		$(LIB_OBJS) $(LDFLAGS) $(LINKLIBS) $(LDLIBS) -o libJSeqArray.$(LIBEXT)


## Clean and remove files
clean:
	$(RM) *.o
	$(RM) libJSeqArray.$(LIBEXT)




##########################################################################

GetData.o: GetData.cpp
	$(CXX) $(CXXFLAGS) GetData.cpp -c -o $@

Index.o: Index.cpp Index.h
	$(CXX) $(CXXFLAGS) Index.cpp -c -o $@

JSeqArray.o: JSeqArray.cpp
	$(CXX) $(CXXFLAGS) JSeqArray.cpp -c -o $@

LinkGDS.o: LinkGDS.c
	$(CC) $(CFLAGS) LinkGDS.cpp -c -o $@

ReadByVariant.o: ReadByVariant.cpp ReadByVariant.h
	$(CXX) $(CXXFLAGS) ReadByVariant.cpp -c -o $@

vectorization.o: vectorization.c vectorization.h
	$(CC) $(CFLAGS) vectorization.c -c -o $@

