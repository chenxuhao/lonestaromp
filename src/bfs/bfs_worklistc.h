#define BFS_VARIANT "worklistc"
#include "worklistc.h"
#include "scan.h"
const int IN_CORE = 0;// set this to zero to disable global-barrier version

void initialize(foru *dist, unsigned int nv) {
	for (int i = 0; i < nv; i++)
		dist[i] = MYINFINITY;
}

foru processedge2(foru *dist, Graph &graph, unsigned iteration, unsigned edge, unsigned &dst) {
	dst = graph.edgessrcdst[edge];
	if (dst >= graph.nnodes) return 0;
	foru wt = 1;
	if (wt >= MYINFINITY) return 0;
	if(dist[dst] == MYINFINITY) {
		dist[dst] = iteration;
		return MYINFINITY;
	}
	return 0;
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
/*
void expandByCTA(foru *dist, Graph &graph, Worklist2 &inwl, Worklist2 &outwl, unsigned iteration) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	int nn;
	typedef cub::BlockScan<int, BLKSIZE> BlockScan;
	__shared__ int owner;
	__shared__ int shnn;
	int total_inputs = (*inwl.dindex + gridDim.x * blockDim.x - 1)/(gridDim.x * blockDim.x);
	owner = -1;
	while(total_inputs-- > 0)
	{      
		int neighborsize = 0;
		int neighboroffset = 0;
		int nnsize = 0;

		if(inwl.pop_id(id, nn))
		{	  
			neighborsize = nnsize = graph.getOutDegree(nn);
			neighboroffset = tex1Dfetch(row_offsets, graph.srcsrc[nn]);	  
		}
		while(true)
		{
			if(nnsize > BLKSIZE)
				owner = threadIdx.x;
			__syncthreads();
			if(owner == -1)
				break;
			if(owner == threadIdx.x)
			{
				shnn = nn;
				cub::ThreadStore<cub::STORE_CG>(inwl.dwl + id, -1);
				owner = -1;
				nnsize = 0;
			}
			__syncthreads();
			neighborsize = graph.getOutDegree(shnn);
			neighboroffset = tex1Dfetch(row_offsets, graph.srcsrc[shnn]);
			int xy = ((neighborsize + blockDim.x - 1) / blockDim.x) * blockDim.x;
			for(int i = threadIdx.x; i < xy; i+= blockDim.x)
			{
				int ncnt = 0;
				unsigned to_push = 0;

				if(i < neighborsize)
					if(processedge2(dist, graph, iteration, neighboroffset + i, to_push))
					{
						ncnt = 1;
					}
				outwl.push_1item<BlockScan>(ncnt, (int) to_push, BLKSIZE);
			}
		}
		id += gridDim.x * blockDim.x;
	}
}
*/
unsigned processnode2(foru *dist, Graph &graph, Worklist2 &inwl, Worklist2 &outwl, unsigned iteration) {
/*
	//expandByCTA(dist, graph, inwl, outwl, iteration);
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	int nn;
	typedef cub::BlockScan<int, BLKSIZE> BlockScan;
	const int SCRATCHSIZE = BLKSIZE;
	BlockScan::TempStorage temp_storage;
	int gather_offsets[SCRATCHSIZE];
	gather_offsets[threadIdx.x] = 0;
	int total_inputs = (*inwl.dindex + gridDim.x * blockDim.x - 1)/(gridDim.x * blockDim.x);
	while(total_inputs-- > 0) {
		int neighborsize = 0;
		int neighboroffset = 0;
		int scratch_offset = 0;
		int total_edges = 0;
		if(inwl.pop_id(id, nn)) {
			if(nn != -1) {
				neighborsize = graph.getOutDegree(nn);
				neighboroffset = tex1Dfetch(row_offsets, graph.srcsrc[nn]);
			}
		}
		BlockScan(temp_storage).ExclusiveSum(neighborsize, scratch_offset, total_edges);
		int done = 0;
		int neighborsdone = 0;
		while(total_edges > 0) {
			__syncthreads();
			int i;
			for(i = 0; neighborsdone + i < neighborsize && (scratch_offset + i - done) < SCRATCHSIZE; i++)
			{
				gather_offsets[scratch_offset + i - done] = neighboroffset + neighborsdone + i;
			}
			neighborsdone += i;
			scratch_offset += i;
			__syncthreads();
			int ncnt = 0;
			unsigned to_push = 0;
			if(threadIdx.x < total_edges) {
				if(processedge2(dist, graph, iteration, gather_offsets[threadIdx.x], to_push)) {
					ncnt = 1;
				}
			}
			outwl.push_1item<BlockScan>(ncnt, (int) to_push, BLKSIZE);
			total_edges -= BLKSIZE;
			done += BLKSIZE;
		}
		id += blockDim.x * gridDim.x;
	}
	return 0;
*/
}

bool processnode(foru *dist, Graph &graph, Worklist2 *inwl, Worklist2 *outwl, unsigned bid, unsigned tid) {
	int neighbours[256];
	int ncnt = 0;
	int nn;
	int id = bid * num_omp_threads + tid;
	//printf("bid=%d, tid=%d, id=%d\n", bid, tid, id);
	if(inwl->pop_id(id, nn)) {
		unsigned neighborsize = graph.getOutDegree(nn);
		if(neighborsize > 256)
			printf("whoa! out of local space");
		for(unsigned i = 0; i < neighborsize; ++ i) {
			unsigned dst = graph.nnodes;
			foru olddist = processedge(dist, graph, nn, i, dst);
			if(olddist) {
				neighbours[ncnt] = dst;
				ncnt++;
			}
		}
	}
	printf("t%d: ncnt=%d\n", id, ncnt);
	for(int i = 0; i < ncnt; i ++)
		printf("t%d: neighbours[%d] = %d\n", id, i, neighbours[i]);
	return outwl->push_nitems<BlockScan>(ncnt, neighbours, tid) == 0 && ncnt > 0;
}

void relax(foru *dist, Graph& graph, unsigned *errno, Worklist2 *inwl, Worklist2 *outwl, int iteration) {
	if(iteration == 0) {
		inwl->push(0);
		return;
	}
	unsigned nitems = inwl->nitems();
	unsigned limit = (nitems + num_omp_threads - 1) / num_omp_threads * num_omp_threads;
	unsigned count = 0;
	#pragma omp parallel for
	for(unsigned i = 0; i < limit; ++i) {
		int tid = omp_get_thread_num();
		//if(processnode2(dist, graph, inwl, outwl, iteration))
		if(processnode(dist, graph, inwl, outwl, count/num_omp_threads, tid)) {
			//printf("out worklist size: %d\n", outwl->nitems());
			*errno = 1;
		}
		my_fetch_add<unsigned>(&count, 1);
		__syncthreads();
	}
}

void relax3(foru *dist, Graph graph, unsigned *errno, Worklist2 *inwl, Worklist2 *outwl, int iteration) {
	relax(dist, graph, errno, inwl, outwl, iteration);
}

void relax2(foru *dist, Graph graph, unsigned *errno, Worklist2 *inwl, Worklist2 *outwl, int iteration) {
	if(iteration == 0)
		relax(dist, graph, errno, inwl, outwl, iteration);
	else {
		Worklist2 *tmp;
		while(inwl->index > 0) {
			relax(dist, graph, errno, inwl, outwl, iteration);
			__syncthreads();
			tmp = inwl;
			inwl = outwl;
			outwl = tmp;
			outwl->index = 0;
			iteration++;
		}
	}
}

void bfs(Graph &graph, foru *dist) {
	int iteration = 0;
	unsigned nerr;
	double starttime, endtime;
	double runtime;
	initialize(dist, graph.nnodes);
	dist[0] = 0;
	Worklist2 wl1(graph.nedges * 2), wl2(graph.nedges * 2);
	Worklist2 *inwl = &wl1, *outwl = &wl2;
	int nitems = 1;
	printf("solving.\n");
	printf("starting...\n");
	starttime = rtclock();
	if(IN_CORE) {
		relax2(dist, graph, &nerr, inwl, outwl, 0);
		relax2(dist, graph, &nerr, inwl, outwl, 1);
	}
	else {
		relax3(dist, graph, &nerr, inwl, outwl, 0);
		nitems = inwl->nitems();
		printf("in worklist size: %d\n", nitems);
		while(nitems > 0) {
			++iteration;
			printf("ITERATION: %d\n", iteration);
			inwl->display_items();
			relax3(dist, graph, &nerr, inwl, outwl, iteration);
			nitems = outwl->nitems();
			printf("out worklist size: %d\n", nitems);
			Worklist2 *tmp = inwl;
			inwl = outwl;
			outwl = tmp;
			outwl->reset();
			if(iteration>4) return;
		};
	}
	endtime = rtclock();
	printf("\titerations = %d.\n", iteration);
	runtime = (1000.0f * (endtime - starttime));
	printf("\truntime [%s] = %f ms.\n", BFS_VARIANT, runtime);
	return;
}
