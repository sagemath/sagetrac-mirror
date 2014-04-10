# CILK_ROOT must contains the GCC/Cilk root directory
OS	 		= $(shell uname) 
#CXX         = $(CILK_ROOT)/bin/g++
CXX         = g++
CXX         = clang++
CPPFLAGS    = -DNDEBUG -DMAX_GENUS=$(MAX_GENUS)
#CXXFLAGS    = -std=c++11 -fcilkplus -g -Wall -O3 # -Winline
CXXFLAGS    = -std=c++11 -g -Wall -O3 # -Winline

ifeq ($(shell uname), Darwin)
TARGET_ARCH = -march=corei7 -march=corei7
LDFLAGS 	= -lcilkrts
else
TARGET_ARCH = -march=native -mtune=native
LDFLAGS     = -lcilkrts -Wl,-rpath=$(CILK_ROOT)/lib64
endif

TARGET = canon.so

SAGE = sage
PYTHON = $(SAGE) -python


all: $(TARGET)

canon.so: borie.pxd borie.hpp canon.pxd canon.pyx
	$(PYTHON) setup.py build_ext --inplace

# borie.o: borie.cpp

clean:
	rm -rf $(TARGET) *.o build canon.cpp canon.html sage cmp test.py 

test: all
	$(SAGE) test.sage
