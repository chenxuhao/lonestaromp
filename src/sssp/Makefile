VARIANTS=sssp-base sssp-atomic sssp-wln sssp-wlc
CXX=g++
CFLAGS=-O3
EXTRA=-fopenmp -DENABLE_OPENMP
SRC=main.cc

all:
	make -f makefile.serial
	make -f makefile.omp

sssp-wln:
	$(CXX) $(CFLAGS) $(EXTRA) -DVARIANT=1 $(SRC) -o sssp-wln

sssp-wlc:
	$(CXX) $(CFLAGS) $(EXTRA) -DVARIANT=2 $(SRC) -o sssp-wlc

sssp-atomic:
	$(CXX) $(CFLAGS) $(EXTRA) -DVARIANT=3 $(SRC) -o sssp-atomic

clean:
	rm sssp-*
