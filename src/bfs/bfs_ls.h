#define BFS_VARIANT "topology-base"

void initialize(foru *dist, unsigned int nv) {
	for (int i = 0; i < nv; i++)
		dist[i] = MYINFINITY;
}

bool processedge(foru *dist, Graph &graph, unsigned src, unsigned ii, unsigned &dst) {
	dst = graph.getDestination(src, ii);
	if (dst >= graph.nnodes) return 0;
	foru wt = 1;
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
	bool changed = false;
	unsigned neighborsize = graph.getOutDegree(work);
	for (unsigned i = 0; i < neighborsize; ++ i) {
		unsigned dst = graph.nnodes;
		if (processedge(dist, graph, work, i, dst))
			changed = true;
	}
	return changed;
}

void drelax(foru *dist, Graph graph, bool *changed) {
	#pragma omp parallel for
	for (unsigned i = 0; i < graph.nnodes; ++i) {
		if (processnode(dist, graph, i)) {
			*changed = true;
		}
	}
}

void bfs(Graph &graph, foru *dist) {
	double starttime, endtime;
	bool changed;
	int iteration = 0;
	printf("initializing.\n");
	initialize(dist, graph.nnodes);
	dist[0] = 0;
	printf("solving.\n");
	starttime = rtclock();
	do {
		++iteration;
		//printf("Iteration %d\n", iteration);
		changed = false;
		drelax(dist, graph, &changed);
	} while (changed);
	endtime = rtclock();
	printf("\titerations = %d.\n", iteration);
	printf("\truntime [%s] = %f ms.\n", BFS_VARIANT, 1000 * (endtime - starttime));
}
