#include <stdio.h>
int main() {
	float summag = 0.555002689;
	printf("summag=%.9f\n", summag);
	float maxmag = 0.560288846;
	printf("maxmag=%.9f\n", maxmag);
	summag += maxmag;
	printf("summag=%.9f\n", summag);
	return 0;
}
