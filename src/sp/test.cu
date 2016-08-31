#include <stdio.h>
__global__ void test(float *g_summag){
	float maxmag = 0.560288846;
	printf("maxmag=%.9f\n", maxmag);
	printf("g_summag=%.9f\n", *g_summag);
	*g_summag += maxmag;
	printf("g_summag=%.9f\n", *g_summag);
}

int main() {
	float summag = 0.555002689;
	printf("summag=%.9f\n", summag);
	float *g_summag;
	cudaMalloc(&g_summag, sizeof(float));
	cudaMemcpy(g_summag, &summag, sizeof(summag), cudaMemcpyHostToDevice);
	test<<<1, 1>>>(g_summag);
	cudaMemcpy(&summag, g_summag, sizeof(summag), cudaMemcpyDeviceToHost);
	printf("summag=%.9f\n", summag);
	return 0;
}
