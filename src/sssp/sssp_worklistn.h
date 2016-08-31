#define SSSP_VARIANT "worklistN"
#include "worklist.h"

void initialize(foru *dist, unsigned int nv) {
	for (int i = 0; i < nv; i++)
		dist[i] = MYINFINITY;
}

foru processedge(foru *dist, Graph &graph, unsigned src, unsigned ii, unsigned &dst) {
	dst = graph.getDestination(src, ii);
	if (dst >= graph.nnodes) return 0;
	foru wt = graph.getWeight(src, ii);
	if (wt >= MYINFINITY) return 0;
	 foru altdist = dist[src] + wt;
	 if (altdist < dist[dst]) {
	 	foru olddist = atomicMin(&dist[dst], altdist);
		if (altdist < olddist) {
			return olddist;
		} 
		// someone else updated distance to a lower value.
	 }
	  return 0;
}

unsigned processnode(foru *dist, Graph &graph, Worklist *outwl, unsigned work) {
	unsigned nn = work;	// wl.getItem(wlii);
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
	for (unsigned i = 0; i < end; ++ i) {
		unsigned work = inwl->getItem(i);
		if (processnode(dist, graph, outwl, work)) {// buffer oveflow.
			*err = 1;
#ifndef ENABLE_OPENMP
			return;
#endif
		}
	}
}
#define SWAP(a, b)	{ tmp = a; a = b; b = tmp; }

void sssp(foru *dist, Graph &graph, long unsigned totalcommu) {
	int iteration = 0;
	Worklist inwl, outwl, *inwlptr, *outwlptr, *tmp;
	unsigned nerr;
	double starttime, endtime;
	double runtime;

	dist[0] = 0;
	unsigned wlsz = 0;
	inwl.ensureSpace(graph.nnodes / 5);
	outwl.ensureSpace(graph.nnodes / 5);
	inwl.push(0);
	inwlptr = &inwl;
	outwlptr = &outwl;
	printf("solving.\n");
	starttime = rtclock();
	do {
		++iteration;
		nerr = 0;
		drelax(dist, graph, inwlptr, outwlptr, &nerr);
		wlsz = outwlptr->getSize();
		if (nerr == 0) {
			if (iteration % 500 == 0) printf("iteration=%d, outwl.size=%d.\n", iteration, wlsz);
			SWAP(inwlptr, outwlptr);
			outwlptr->noverflows = inwlptr->noverflows;
		} else {	// error: currently only buffer oveflow.
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
	printf("\truntime [%s] = %f ms.\n", SSSP_VARIANT, runtime);
}
