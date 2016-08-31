VARIANTS=bfs bfs-wlw bfs-wlc
CXX=g++
CFLAGS=-O3
EXTRA=-fopenmp -DENABLE_OPENMP
SRC=main.cc

all:
	make -f makefile.serial
	make -f makefile.omp


bfs-wlw:
	$(CXX) $(CFLAGS) $(EXTRA) -DVARIANT=3 $(SRC) -o bfs-wlw

bfs-wlc:
	$(CXX) $(CFLAGS) $(EXTRA) -DVARIANT=6 $(SRC) -o bfs-wlc

clean:
	rm bfs-*
