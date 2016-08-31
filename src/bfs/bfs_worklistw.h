#define BFS_VARIANT "worklistw"
#include "worklist.h"

void initialize(foru *dist, unsigned int nv) {
	for (int i = 0; i < nv; i++)
		dist[i] = MYINFINITY;
}

foru processedge(foru *dist, Graph &graph, unsigned src, unsigned ii, unsigned &dst) {
	dst = graph.getDestination(src, ii);
	if (dst >= graph.nnodes) return 0;
	foru wt = 1;
	if (wt >= MYINFINITY) return 0;
	foru altdist = dist[src] + wt;
	if(dist[dst] == MYINFINITY) {
		dist[dst] = altdist;
		return MYINFINITY;
	}
	return 0;
}

unsigned processnode(foru *dist, Graph &graph, Worklist *outwl, unsigned work) {
	unsigned nn = work;
	if (nn >= graph.nnodes) return 0;
	unsigned neighborsize = graph.getOutDegree(nn);
	for (unsigned i = 0; i < neighborsize; ++ i) {
		unsigned dst = graph.nnodes;
		foru olddist = processedge(dist, graph, nn, i, dst);
		if (olddist) {
			if (outwl->push(dst)) {	// buffer oveflow.
				dist[dst] = olddist;	// no atomicMax required.
				return 1;
			}
		}
	}
	return 0;
}

void drelax(foru *dist, Graph graph, Worklist *inwl, Worklist *outwl, unsigned *err) {
	unsigned end = inwl->count();
	#pragma omp parallel for
	for (unsigned i = 0; i < end; ++i) {
		unsigned work = inwl->getItem(i);
		if (processnode(dist, graph, outwl, work)) {
			*err = 1;
#ifndef ENABLE_OPENMP
			return;
#endif
		}
	}
}
#define SWAP(a, b)	{ tmp = a; a = b; b = tmp; }

void bfs(Graph &graph, foru *dist) {
	int iteration = 0;
	Worklist inwl, outwl, *inwlptr, *outwlptr, *tmp;
	unsigned nerr;
	double starttime, endtime;
	double runtime;
	initialize(dist, graph.nnodes);
	dist[0] = 0;
	unsigned wlsz = 0;
	inwl.ensureSpace(13421760);
	outwl.ensureSpace(13421760);
	inwl.push(0);	// source.
	inwlptr = &inwl;
	outwlptr = &outwl;
	printf("solving.\n");
	starttime = rtclock();
	do {
		++iteration;
//		printf("Iteration %d: %d nodes in the input worklist\n", iteration, inwlptr->count());
		nerr = 0;
		drelax(dist, graph, inwlptr, outwlptr, &nerr);
		wlsz = outwlptr->getSize();
//		printf("%d nodes in the output worklist\n", wlsz);
		if (nerr == 0) {
			SWAP(inwlptr, outwlptr);
			outwlptr->noverflows = inwlptr->noverflows;
		} else {
			printf("Error: currently only buffer oveflow\n");
			if (++outwlptr->noverflows == MAXOVERFLOWS) {
				unsigned cap = inwlptr->getCapacity();
				inwlptr->ensureSpace(2 * cap);	// double the capacity.
				outwlptr->ensureSpace(2 * cap);
				inwlptr->append(*outwlptr);
				outwlptr->noverflows = 0;
			} else {
				// defer increasing worklist capacity.
				printf("\tdeferred increasing worklist capacity.\n");
			}
		}
		outwlptr->clear();	// clear it whether overflow or not.
	} while (wlsz);
	endtime = rtclock();
	printf("\titerations = %d.\n", iteration);
	runtime = (1000.0f * (endtime - starttime));
	printf("\truntime [%s] = %f ms.\n", BFS_VARIANT, runtime);
	return;
}
