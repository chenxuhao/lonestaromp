/** Breadth-first search -*- C++ -*-
 * @author Xuhao Chen  <cxh.nudt@gmail.com>
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
//#include "bfs_ls.h"

int num_omp_threads;

void verifysolution(foru *dist, Graph graph, unsigned *nerr) {
	for (int nn = 0; nn < graph.nnodes; nn ++) {
		unsigned int nsrcedges = graph.getOutDegree(nn);
		for (unsigned ii = 0; ii < nsrcedges; ++ii) {
			unsigned int u = nn;
			unsigned int v = graph.getDestination(u, ii);
			foru wt = 1;
			if (wt > 0 && dist[u] + wt < dist[v]) {
				//printf("%d %d %d %d\n", u, v, dist[u], dist[v]);
				++*nerr;
			}
		}
	}	
}

void write_solution(const char *fname, Graph &graph, foru *dist) {
	printf("Writing solution to %s\n", fname);
	FILE *f = fopen(fname, "w");
	// formatted like Merrill's code for comparison
	fprintf(f, "Computed solution (source dist): [");

	for(int node = 0; node < graph.nnodes; node++)
		fprintf(f, "%d:%d\n ", node, dist[node]);
	fprintf(f, "]");
}

int main(int argc, char *argv[]) {
	unsigned intzero = 0;
	Graph graph;
	unsigned nerr = 0;
	foru *dist;

	if (argc != 3) {
		printf("Usage: %s <nThreads> <graph>\n", argv[0]);
		exit(1);
	}
	num_omp_threads = atoi(argv[1]);
#ifdef ENABLE_OPENMP
	omp_set_num_threads(num_omp_threads);
#endif
	graph.read(argv[2]);
	dist = (foru *) malloc(graph.nnodes * sizeof(foru));
	memset(dist, 0, graph.nnodes * sizeof(foru));
	bfs(graph, dist);
	printf("verifying.\n");
	verifysolution(dist, graph, &nerr);
	printf("\tno of errors = %d.\n", nerr);
	write_solution("bfs-output.txt", graph, dist);
	free(dist);
	return 0;
}
