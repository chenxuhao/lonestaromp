#pragma once
#include "sharedptr.h"
typedef unsigned int uint;
typedef struct {
	uint x;
	uint y;
	uint z;
}uint3;
#define MINANGLE	30
#define PI		3.14159265358979323846	// from C99 standard.
#define FORD		double
#define DIMSTYPE	unsigned

#define INVALIDID	1234567890
#define MAXID		INVALIDID

// "usual" ratio of final nodes to final elements, determined empirically
// used to adjust maxfactor for nodes
// use 1 to be conservative
#define MAX_NNODES_TO_NELEMENTS 2
#define INIT_OWNER -1
struct Mesh {
	uint maxnelements;
	uint maxnnodes;
	uint ntriangles;
	uint nnodes;
	uint nsegments;
	uint nelements;
	FORD *nodex; // could be combined
	FORD *nodey;
	uint3 *elements;
	volatile bool *isdel;  
	bool *isbad;
	uint3 *neighbours;
	int *owners;

	Mesh() {}
	void alloc() {
		printf("Allocating memory for mesh elements...\n");
		elements = (uint3 *) calloc(maxnelements, sizeof(uint3));
		neighbours = (uint3 *) calloc(maxnelements, sizeof(uint3));
		isdel = (bool *) calloc(maxnelements, sizeof(bool));
		isbad = (bool *) calloc(maxnelements, sizeof(bool));
		owners = (int *) calloc(maxnelements, sizeof(int));
	}

	void init() {
		for(int i=0; i<maxnelements; i++) {
			owners[i] = INIT_OWNER;
		}
	}
	void free() {
	}
};

#define IS_SEGMENT(element) (((element).z == INVALIDID))
