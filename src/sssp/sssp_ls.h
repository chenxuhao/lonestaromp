#pragma once
#define SSSP_VARIANT "topology_base"

void initialize(foru *dist, unsigned int nv) {
	for (int i = 0; i < nv; i++)
		dist[i] = MYINFINITY;
}

bool processedge(foru *dist, Graph &graph, unsigned src, unsigned ii, unsigned &dst) {
	dst = graph.getDestination(src, ii);
	if (dst >= graph.nnodes) return 0;
	foru wt = graph.getWeight(src, ii);
	if (wt >= MYINFINITY) return 0;
	foru altdist = dist[src] + wt;
	if (altdist < dist[dst]) {
#ifndef ENABLE_OPENMP
		dist[dst] = altdist;
		return true;
#else
	 	foru olddist = atomicMin(&dist[dst], altdist);
		if (altdist < olddist) {
			return true;
		} 
		// someone else updated distance to a lower value.
#endif
	}
	return false;
}

bool processnode(foru *dist, Graph &graph, unsigned work) {
	unsigned nn = work;
	if (nn >= graph.nnodes) return 0;
	bool changed = false;
	unsigned neighborsize = graph.getOutDegree(nn);
	for (unsigned i = 0; i < neighborsize; ++ i) {
		unsigned dst = graph.nnodes;
		if(processedge(dist, graph, nn, i, dst))
			changed = true;
	}
	return changed;
}

void drelax(foru *dist, Graph graph, bool *changed) {
	#pragma omp parallel for
	for (unsigned i = 0; i < graph.nnodes; ++ i) {
		if (processnode(dist, graph, i)) {
			*changed = true;
		}
	}
}

void sssp(foru *dist, Graph &graph, long unsigned totalcommu) {
	double starttime, endtime;
	bool changed;
	int iteration = 0;
	dist[0] = 0;
	printf("solving.\n");
	starttime = rtclock();
	do {
		++iteration;
		changed = false;
		drelax(dist, graph, &changed);
	} while (changed);
	endtime = rtclock();
	totalcommu += graph.nnodes * sizeof(foru);
	printf("\titerations = %d communication = %.3lf MB.\n", iteration, totalcommu * 1.0 / 1000000);
	printf("\truntime [%s] = %f ms\n", SSSP_VARIANT, 1000 * (endtime - starttime));
}
