/** Single source shortest paths -*- C++ -*-
 * @author Xuhao Chen <cxh.nudt@gmail.com>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cassert>
#include <inttypes.h>
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "common.h"
typedef unsigned foru;
#include "graph.h"
#include "variants.h"
int num_omp_threads;

void verifysolution(foru *dist, Graph graph, unsigned *nerr) {
	for (int nn = 0; nn < graph.nnodes; nn ++) {
		unsigned int nsrcedges = graph.getOutDegree(nn);
		for (unsigned ii = 0; ii < nsrcedges; ++ii) {
			unsigned int u = nn;
			unsigned int v = graph.getDestination(u, ii);
			foru wt = graph.getWeight(u, ii);
			if (wt > 0 && dist[u] + wt < dist[v]) {
				++*nerr;
			}
		}
	}	
}

void print_output(const char *filename, foru *dist, Graph graph) {
	printf("Writing output to %s\n", filename);
	FILE *o = fopen(filename, "w");
	for(int i = 0; i < graph.nnodes; i++) {
		fprintf(o, "%d: %d\n", i, dist[i]);
	}
	fclose(o);
}

int main(int argc, char *argv[]) {
	foru *dist;
	Graph graph;
	unsigned nerr = 0;
	if (argc != 3) {
		printf("Usage: %s <nThreads> <graph>\n", argv[0]);
		exit(1);
	}
	num_omp_threads = atoi(argv[1]);
#ifdef ENABLE_OPENMP
	omp_set_num_threads(num_omp_threads);
#endif
	graph.read(argv[2]);
	long unsigned totalcommu;
	dist = (foru *)malloc(graph.nnodes * sizeof(foru));
	printf("initializing.\n");
	initialize(dist, graph.nnodes);
	sssp(dist, graph, totalcommu);
	printf("verifying.\n");
	verifysolution(dist, graph, &nerr);
	printf("\tno of errors = %d.\n", nerr);
	print_output("sssp-output.txt", dist, graph);
	free(dist);
	return 0;
}
