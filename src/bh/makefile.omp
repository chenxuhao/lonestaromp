CXX=g++ 
CXX_FLAGS=-O3 -fopenmp -DENABLE_OPENMP
SRC=main.cc
EXE=bh-omp
INC=-I../../include

all:
	$(CXX) $(CXX_FLAGS) $(INC) $(SRC) -o $(EXE)
	mv $(EXE) ../../bin
