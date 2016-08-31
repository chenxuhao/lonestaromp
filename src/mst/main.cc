/** Minimum spanning tree -*- C++ -*-
 * Computes minimum spanning tree of a graph using Boruvka's algorithm.
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
#include "component.h"
#include "devel.h"

int num_omp_threads;

void init(Graph graph, ComponentSpace cs, foru *eleminwts, 
		foru *minwtcomponent, unsigned *partners, 
		bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
#ifdef TIMING
	double starttime = rtclock();
#endif
	#pragma omp parallel for schedule(static)
	for (int id=0; id < graph.nnodes; id++) {
		eleminwts[id] = MYINFINITY;
		minwtcomponent[id] = MYINFINITY;
		goaheadnodeofcomponent[id] = graph.nnodes;
		partners[id] = id;
		processinnextiteration[id] = false;
	}
#ifdef TIMING
	double endtime = rtclock();
	printf("\truntime [init] = %f ms.\n", 1000 * (endtime - starttime));
#endif
}

//find the minimum edge of each node whose other endpoint is in another component
void findelemin(Graph graph, ComponentSpace cs, foru *eleminwts, 
		foru *minwtcomponent, unsigned *partners, 
		bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
#ifdef TIMING
	double starttime = rtclock();
#endif
	#pragma omp parallel for schedule(static)
	for(int id = 0; id < graph.nnodes; id ++) {
		// if I have a cross-component edge,
		// find my minimum wt cross-component edge,
		// inform my boss about this edge e (atomicMin).
		unsigned src = id;
		unsigned srcboss = cs.find(src);
		unsigned dstboss = graph.nnodes;
		foru minwt = MYINFINITY;
		unsigned degree = graph.getOutDegree(src);
		for (unsigned ii = 0; ii < degree; ++ii) {
			foru wt = graph.getWeight(src, ii);
			if (wt < minwt) {
				unsigned dst = graph.getDestination(src, ii);
				unsigned tempdstboss = cs.find(dst);
				if (srcboss != tempdstboss) {
					minwt = wt;
					dstboss = tempdstboss;
				}
			}
		}
		//printf("\tminwt[%d] = %d\n", id, minwt);
		eleminwts[id] = minwt;
		partners[id] = dstboss;

		if (minwt < minwtcomponent[srcboss] && srcboss != dstboss) {
			// inform boss.
			#pragma omp critical
			{
//			minwtcomponent[srcboss] = minwt < minwtcomponent[srcboss] ? minwt : minwtcomponent[srcboss];
			if (minwt < minwtcomponent[srcboss])
				minwtcomponent[srcboss] = minwt;
			}
		}
	}
//	for(int id = 0; id < 8; id ++)
//		printf("%d[%d]->%d (%d)\n", id, graph.getOutDegree(id), partners[id], eleminwts[id]);
#ifdef TIMING
	double endtime = rtclock();
	printf("\truntime [findelemin] = %f ms.\n", 1000 * (endtime - starttime));
#endif
}

// put the minimum weight node in goaheadnodeofcomponent[]
void findelemin2(Graph graph, ComponentSpace cs, foru *eleminwts, 
		foru *minwtcomponent, unsigned *partners, 
		bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
#ifdef TIMING
	double starttime = rtclock();
#endif
	#pragma omp parallel for schedule(static)
	for (int id=0; id < graph.nnodes; id++) {
		unsigned src = id;
		unsigned srcboss = cs.find(src);
		if(eleminwts[id] == minwtcomponent[srcboss] && 
				srcboss != partners[id] && partners[id] != graph.nnodes) {
			unsigned degree = graph.getOutDegree(src);
			for (unsigned ii = 0; ii < degree; ++ii) {
				foru wt = graph.getWeight(src, ii);
				if (wt == eleminwts[id]) {
					unsigned dst = graph.getDestination(src, ii);
					unsigned tempdstboss = cs.find(dst);
					if (tempdstboss == partners[id]) {
						my_compare_swap<unsigned>(&goaheadnodeofcomponent[srcboss], graph.nnodes, id);
					}
				}
			}
		}
	}
//	for(int id = 0; id < 8; id ++)
//		printf("go_ahead:%d:%d\n", id, goaheadnodeofcomponent[id]);
#ifdef TIMING
	double endtime = rtclock();
	printf("\truntime [findelemin2] = %f ms.\n", 1000 * (endtime - starttime));
#endif
}

void verify_min_elem(Graph graph, ComponentSpace cs, foru *eleminwts, 
		foru *minwtcomponent, unsigned *partners, 
		bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
#ifdef TIMING
	double starttime = rtclock();
#endif
	#pragma omp parallel for schedule(static)
	for (int id=0; id < graph.nnodes; id++) {
		if(cs.isBoss(id)) {
			if(goaheadnodeofcomponent[id] == graph.nnodes) {
				continue;
			}
			unsigned minwt_node = goaheadnodeofcomponent[id];
			unsigned degree = graph.getOutDegree(minwt_node);
			foru minwt = minwtcomponent[id];
			if(minwt == MYINFINITY)
				continue;
			bool minwt_found = false;
			//printf("%d: looking at %d def %d minwt %d\n", id, minwt_node, degree, minwt);
			for (unsigned ii = 0; ii < degree; ++ii) {
				foru wt = graph.getWeight(minwt_node, ii);
				//printf("%d: looking at %d edge %d wt %d (%d)\n", id, minwt_node, ii, wt, minwt);

				if (wt == minwt) {
					minwt_found = true;
					unsigned dst = graph.getDestination(minwt_node, ii);
					unsigned tempdstboss = cs.find(dst);
					if(tempdstboss == partners[minwt_node] && tempdstboss != id) {
						processinnextiteration[minwt_node] = true;
						break;
					}
				}
			}
			assert(minwt_found);
			//printf("component %d is wrong %d\n", id, minwt_found);
		}
	}
//	for(int id = 0; id < 8; id ++)
//		printf("process:%d:%s\n", id, processinnextiteration[id]?"true":"false");
#ifdef TIMING
	double endtime = rtclock();
	printf("\truntime [verify] = %f ms.\n", 1000 * (endtime - starttime));
#endif
}
/*
void elim_dups(Graph graph, ComponentSpace cs, foru *eleminwts, 
		foru *minwtcomponent, unsigned *partners, 
		bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
	for (int id=0; id < graph.nnodes; id++) {
		if(processinnextiteration[id]) {
			unsigned srcc = cs.find(id);
			unsigned dstc = partners[id];

			if(minwtcomponent[dstc] == eleminwts[id]) {
				if(id < goaheadnodeofcomponent[dstc]) {
					processinnextiteration[id] = false;
					//printf("duplicate!\n");
				}
			}
		}
	}
}

void findcompmin(Graph graph, ComponentSpace cs, foru *eleminwts, 
		foru *minwtcomponent, unsigned *partners, 
		bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
	for (int id=0; id < graph.nnodes; id++) {
		if(partners[id] == graph.nnodes)
			return;
		unsigned srcboss = cs.find(id);
		unsigned dstboss = cs.find(partners[id]);
		if (id != partners[id] && srcboss != dstboss && 
			eleminwts[id] != MYINFINITY && 
			minwtcomponent[srcboss] == eleminwts[id] && dstboss != id && 
			goaheadnodeofcomponent[srcboss] == id) {// my edge is min outgoing-component edge.
			if(!processinnextiteration[id]);
			//printf("whoa!\n");
			//= true;
		}
		else {
			if(processinnextiteration[id]);
			//printf("whoa2!\n");
		}
	}
}
*/
void findcompmintwo(unsigned *mstwt, Graph graph, ComponentSpace &csw, 
		foru *eleminwts, foru *minwtcomponent, unsigned *partners, 
		bool *processinnextiteration, unsigned *goaheadnodeofcomponent, 
		bool *repeat, unsigned *count) {
#ifdef TIMING
	double starttime = rtclock();
#endif
	unsigned up = graph.nnodes;
#ifdef ENABLE_OPENMP
	up = (up+num_omp_threads-1)/num_omp_threads*num_omp_threads;
#endif
//	printf("up = %d\n", up);
//	#pragma omp parallel for
	for(int id = 0; id < up; id ++) {
		unsigned srcboss, dstboss;
		if(id < graph.nnodes && processinnextiteration[id]) {
			srcboss = csw.find(id);
			dstboss = csw.find(partners[id]);
		}
		__syncthreads();
		if(id < graph.nnodes && processinnextiteration[id] && srcboss != dstboss) {
//			printf("trying unify id=%d (%d -> %d)\n", id, srcboss, dstboss);
			if (csw.unify(srcboss, dstboss)) {
#ifdef ENABLE_OPENMP
				my_fetch_add<unsigned>(mstwt, eleminwts[id]);
#else
				*mstwt += eleminwts[id];
#endif
#ifdef ENABLE_OPENMP
				my_fetch_add<unsigned>(count, 1);
#else
				(*count) ++;
#endif
//				printf("u %d -> %d (%d)\n", srcboss, dstboss, eleminwts[id]);
				processinnextiteration[id] = false;
				eleminwts[id] = MYINFINITY;	// mark end of processing to avoid getting repeated.
			}
			else {
				*repeat = true;
			}
//			printf("\tcomp[%d] = %d.\n", srcboss, csw.find(srcboss));
		}
		__syncthreads();
	}
#ifdef TIMING
	double endtime = rtclock();
	printf("\truntime [findcompmin] = %f ms.\n", 1000 * (endtime - starttime));
#endif
}

int main(int argc, char *argv[]) {
	unsigned mstwt = 0;
	int iteration = 0;
	Graph graph;
	unsigned *partners;
	foru *eleminwts, *minwtcomponent;
	bool *processinnextiteration;
	unsigned *goaheadnodeofcomponent;
	double starttime, endtime;
	if (argc != 3) {
		printf("Usage: %s <nThreads> <graph>\n", argv[0]);
		exit(1);
	}
	num_omp_threads = atoi(argv[1]);
	graph.read(argv[2]);
	//graph.print();
	ComponentSpace cs(graph.nnodes);

	eleminwts = (foru *)malloc(graph.nnodes * sizeof(foru));
	minwtcomponent = (foru *)malloc(graph.nnodes * sizeof(foru));
	partners = (unsigned *)malloc(graph.nnodes * sizeof(unsigned));
	processinnextiteration = (bool *)malloc(graph.nnodes * sizeof(bool));
	goaheadnodeofcomponent = (unsigned *)malloc(graph.nnodes * sizeof(unsigned));

	unsigned prevncomponents, currncomponents = graph.nnodes;

	bool repeat = false;
	unsigned edgecount = 0;
#ifdef ENABLE_OPENMP
	omp_set_num_threads(num_omp_threads);
#endif
	printf("finding mst.\n");
	starttime = rtclock();
	do {
		++iteration;
		prevncomponents = currncomponents;
		init(graph, cs, eleminwts, minwtcomponent, partners, processinnextiteration, goaheadnodeofcomponent);
		findelemin(graph, cs, eleminwts, minwtcomponent, partners, processinnextiteration, goaheadnodeofcomponent);
		findelemin2(graph, cs, eleminwts, minwtcomponent, partners, processinnextiteration, goaheadnodeofcomponent);
		verify_min_elem(graph, cs, eleminwts, minwtcomponent, partners, processinnextiteration, goaheadnodeofcomponent);
		if(0) print_comp_mins(cs, graph, minwtcomponent, goaheadnodeofcomponent, partners, processinnextiteration);
		do {
			repeat = false;
			findcompmintwo(&mstwt, graph, cs, eleminwts, minwtcomponent, partners, processinnextiteration, goaheadnodeofcomponent, &repeat, &edgecount);
		} while (repeat); // only required for quicker convergence?
		currncomponents = cs.numberOfComponents();
		printf("\titeration %d, number of components = %d (%d), mstwt = %u mstedges = %u\n", iteration, currncomponents, prevncomponents, mstwt, edgecount);
	} while (currncomponents != prevncomponents);

	endtime = rtclock();
	printf("\tmstwt = %u, iterations = %d.\n", mstwt, iteration);
	printf("\t%s result: weight: %u, components: %u, edges: %u\n", argv[2], mstwt, currncomponents, edgecount);
	printf("\truntime [mst] = %f ms.\n", 1000 * (endtime - starttime));

	cs.deallocate();
	free(eleminwts);
	free(minwtcomponent);
	free(partners);
	free(processinnextiteration);
	free(goaheadnodeofcomponent);
	return 0;
}
