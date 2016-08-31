#ifndef SCAN_H
#define SCAN_H
#include "common.h"
#define MAX_NUM_THREADS 1024
class BlockScan {
private:
	int num_threads;
public:
	static int in[MAX_NUM_THREADS];
	static int out[MAX_NUM_THREADS];
	BlockScan() {}
	BlockScan(int nthreads):
		num_threads(nthreads) {
	}
	void ExclusiveSum(int input, int &output, int &block_aggregate, unsigned id) {
		in[id] = input;
		__syncthreads();
/*
		if(id == 0) {
			for(int i=0; i<num_threads; i++)
				printf("t%d: input=%d ", i, in[i]);
			printf("\n");
		}
*/
		if(id == 0) {
			out[0] = 0;
			for(int i = 1; i < num_threads; i++)
				out[i] = out[i-1] + in[i-1];
		}
		__syncthreads();
		output = out[id];
		//printf("output=%d, id =%d\n", output, id);
		block_aggregate = out[num_threads-1] + in[num_threads-1]; 
	}
};

int BlockScan::in[MAX_NUM_THREADS] = {0};
int BlockScan::out[MAX_NUM_THREADS] = {0};
#endif
