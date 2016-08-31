/* Survey propagation -*- C++ -*-
 * Sequential implementation of the Survey Propagation Algorithm
 * @author Xuhao Chen <cxh.nudt@gmail.com>
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>
#include <math.h>
#include <cassert>
#include "sp-s.h"
#include "db.h"
#include "sort.h"
#include "common.h"
#define EPSILON 0.01 // in the source, in the paper 10^-3
#define MAXITERATION 1000 
#define PARTIAL "partial.cnf"
#define PARAMAGNET 0.01

int num_omp_threads;

void init_from_file(const char *F, 
		int max_lit_per_clause,
		struct CSRGraph &clauses,
		struct CSRGraph &vars, 
		struct Edge &ed) {
	FILE *f = fopen(F, "r");
	int nclauses, nvars, ret;
	char line[255];
	if(!f) {
		fprintf(stderr, "unable to read file %s.\n", F);
		exit(1);
	}
	while(true) {
		if(fgets(line, 255, f)) {
			if(line[0] != 'c')
				break;
			printf("%s", line);
		}
		else {
			fprintf(stderr, "unable to read %s\n", F);
			exit(1);
		}
	}
	ret = sscanf(line, "p cnf %d %d", &nvars, &nclauses);
	assert(ret == 2);
	clauses.nnodes = nclauses;
	vars.nnodes = nvars;
	ed.nedges = clauses.nedges = vars.nedges = nclauses * max_lit_per_clause; // over-estimate
	assert(clauses.alloc());
	assert(vars.alloc());
	assert(ed.alloc());
	int newlit, lit;
	int clndx = 0, litndx = 0, edndx = 0;

	// read lines of literals terminated by 0
	// assumes literals numbered from 1
	do {
		ret = fscanf(f, "%d", &newlit);
		if(ret == EOF) break;
		if(newlit == 0) {
			assert(clndx < nclauses);
			clndx++;
			litndx = 0;
			clauses.row_offsets[clndx] = edndx;
			continue;
		}
		assert(litndx < max_lit_per_clause);
		// convert to zero-based
		lit = ((newlit < 0) ? -newlit : newlit) - 1;
		assert(lit >= 0);
		ed.src[edndx] = clndx;
		ed.dst[edndx] = lit;
		ed.bar[edndx] = newlit < 0;
		ed.eta[edndx] = (float)(rand()) / (float)RAND_MAX;

		// essentially clause -> edge
		clauses.columns[clauses.row_offsets[clndx] + litndx] = edndx;

		// record size of every var node
		vars.row_offsets[lit]++;
		litndx++;
		edndx++;
	} while(true);

	clauses.nedges = vars.nedges = ed.nedges = edndx;

	clauses.set_last_offset();
	vars.set_last_offset();

	// populate vars
	// exclusive-sum
	for(int i = 0, sum = 0; i < vars.nnodes; i++) {
		int size = vars.row_offsets[i];
		vars.row_offsets[i] = sum;
		sum += size;
	}
	int *varndx = (int *) calloc(vars.nedges, sizeof(int));
	for(int i = 0; i < ed.nedges; i++) {
		unsigned var = ed.dst[i];
		vars.columns[vars.row_offsets[var] + varndx[var]++] = i;
	}
	printf("read %d clauses, %d variables, %d literals\n", clauses.nnodes, vars.nnodes, ed.nedges);
}

void print_solution(const char *sol, const CSRGraph &vars) {
	FILE *f = fopen(sol, "w");
	int i;
	for(i = 0; i < vars.nnodes; i++) {
		if(vars.sat[i])
			fprintf(f, "%d\n", vars.value[i] ? (i + 1) : -(i + 1));
	}
	fclose(f);
}

void dump_formula(const char *output, const CSRGraph &clauses, const CSRGraph &vars, const Edge &ed) {
	FILE *of = fopen(output, "w");
	fprintf(of, "p cnf %d %d\n", vars.nnodes, clauses.nnodes);
	for(int cl = 0; cl < clauses.nnodes; cl++) {
		unsigned offset = clauses.row_offsets[cl];
		for(int i = 0; i < clauses.degree(cl); i++) {
			unsigned edndx = clauses.columns[offset + i];
			fprintf(of, "%d ", ed.bar[edndx] ? -(ed.dst[edndx]+1) : (ed.dst[edndx]+1));
		}
		fprintf(of, "0\n");
	}
}

void dump_partial(const char *output, const CSRGraph &clauses, const CSRGraph &vars, const Edge &ed) {
	FILE *of = fopen(output, "w");
	int sat = 0;
	for(int cl = 0; cl < clauses.nnodes; cl++)
		if(clauses.sat[cl]) sat++;
	fprintf(of, "p cnf %d %d\n", vars.nnodes, clauses.nnodes - sat);
	for(int cl = 0; cl < clauses.nnodes; cl++) {
		if(clauses.sat[cl])
			continue;
		unsigned offset = clauses.row_offsets[cl];
		for(int i = 0; i < clauses.degree(cl); i++) {
			unsigned edndx = clauses.columns[offset + i];
			if(vars.sat[ed.dst[edndx]])
				continue;
			fprintf(of, "%d ", ed.bar[edndx] ? -(ed.dst[edndx]+1) : (ed.dst[edndx]+1));
		}
		fprintf(of, "0\n");
	}
}

void calc_pi_values(CSRGraph clauses, CSRGraph vars, Edge ed) {
    double starttime, endtime;
	starttime = rtclock();
	// over all a -> j
	// for each edge
#ifdef ENABLE_OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for(int edndx = 0; edndx < ed.nedges; edndx ++) {
		//get the src and dst nodes id
		int j = ed.dst[edndx];
		int a = ed.src[edndx];

		//if any of them is satisfied, continue
		if(clauses.sat[a] || vars.sat[j])
			continue;

		int V_j = vars.row_offsets[j];//get the starting index of edge id of node j
		int V_j_len = vars.degree(j);//get the number of outgoing edges of node j

		float pi_0 = 1.0;
		float V_s_a = 1.0;
		float V_u_a = 1.0;

		// over all b E V(j)
		//for each outgoing edge of node j
		for(int i = 0; i < V_j_len; i++) {
			int ed_btoj = vars.columns[V_j + i];//get the edge id
			int b = ed.src[ed_btoj];//get the src node of this edge
			if(clauses.sat[b])
				continue;
			if(b != a) {
				pi_0 *= (1 - ed.eta[ed_btoj]);
				if(ed.bar[ed_btoj] == ed.bar[edndx])
					V_s_a *= (1 - ed.eta[ed_btoj]);
				else
					V_u_a *= (1 - ed.eta[ed_btoj]); 
			}
		}
		ed.pi_0[edndx] = pi_0;
		ed.pi_U[edndx] = (1 - V_u_a) * (V_s_a);
		ed.pi_S[edndx] = (1 - V_s_a) * (V_u_a);
		//printf("[%d] %f %f %f\n", edndx, ed.pi_0[edndx], ed.pi_U[edndx], ed.pi_S[edndx]);
	}
    endtime = rtclock();
#ifdef TIMING
    printf("\truntime[calc_pi] = %f ms.\n", 1000*(endtime-starttime));
#endif
}

//iterates over the clauses and the literals of the formula updating 'surveys'
void update_eta(CSRGraph clauses, CSRGraph vars, Edge ed, float *max_eps) {
#ifdef ENABLE_OPENMP
	float *lmaxeps;//local maximum eps
//	printf("numThreads=%d\n", num_omp_threads);
	lmaxeps = (float *)malloc(num_omp_threads*sizeof(float));
	for(int i=0; i<num_omp_threads; i++) {
		lmaxeps[i] = 0;
	}
#else
	float lmaxeps = 0;
#endif
	double starttime, endtime;
	starttime = rtclock();
	//for each edge
#ifdef ENABLE_OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for(int edndx = 0; edndx < ed.nedges; edndx++) {
		float eps;
		int a = ed.src[edndx];//get the src node id
		int i = ed.dst[edndx];//get the dst node id

		// as these are "removed"
		if(clauses.sat[a] || vars.sat[i])
			continue;

		int clndx = clauses.row_offsets[a];//get the starting index of edge id of node a
		int nlit = clauses.degree(a);//get the number of outgoing edges of node a

		float new_eta = 1.0;

		//for each outgoing edge of node a
		for(int aedndx = 0; aedndx < nlit; aedndx++) {
			int jedndx = clauses.columns[clndx + aedndx];//get the edge id 
			int j = ed.dst[jedndx];//get the dst node of this edge
			if(j == i)
				continue;
			if(vars.sat[j])
				continue;
			float sum = ed.pi_0[jedndx] + ed.pi_S[jedndx] + ed.pi_U[jedndx];
			if(sum == 0.0) { // TODO: non-standard ...
				new_eta = 0;
				break;
			}
			new_eta *= ed.pi_U[jedndx] / sum;//calculate the new eta
		}
		eps = fabs(new_eta - ed.eta[edndx]);
#ifdef ENABLE_OPENMP
		int tid = omp_get_thread_num();
		if(eps > lmaxeps[tid])
			lmaxeps[tid] = eps;
#else
		if(eps > lmaxeps)
			lmaxeps = eps;
#endif
		ed.eta[edndx] = new_eta;//update eta
	}

#ifdef ENABLE_OPENMP
	for(int i=0; i<num_omp_threads; i++) {
		if(lmaxeps[i]>*max_eps)
			*max_eps = lmaxeps[i];
	}
	free(lmaxeps);
#else
	if(lmaxeps>*max_eps)
		*max_eps = lmaxeps;
#endif
	endtime = rtclock();
#ifdef TIMING
	printf("\truntime[update_eta] = %f ms.\n", 1000*(endtime-starttime));
#endif
}

void update_bias(CSRGraph clauses, CSRGraph vars, Edge ed, float *bias_list, int *bias_list_vars, int *bias_len, float &summag) {
    double starttime, endtime;
#ifdef ENABLE_OPENMP
	float *lsummag = (float *) malloc(num_omp_threads*sizeof(float));
	for(int i=0; i<num_omp_threads; i++)
		lsummag[i] = 0.0f;
#endif
	starttime = rtclock();
	//for each node
#ifdef ENABLE_OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for(int v = 0; v < vars.nnodes; v++) {
		float maxmag = 0.0f;
		if(vars.sat[v])
			continue;
		float pi_0_hat = 1.0, V_minus = 1.0, V_plus = 1.0;
		float pi_P_hat, pi_M_hat;
		int edoff = vars.row_offsets[v];//get the starting index of edge id of node v
		int ncl = vars.degree(v);//get the number of outgoing edges of node v

		//for each outgoing edge of node v 
		// a E v(i)
		for(int edndx = 0; edndx < ncl; edndx++) {
			int edge = vars.columns[edoff + edndx];//get the edge id
			int cl = ed.src[edge];//get the src node of this edge
			if(clauses.sat[cl])
				continue;
			pi_0_hat *= (1 - ed.eta[edge]);
			if(ed.bar[edge])
				V_minus *= (1 - ed.eta[edge]);
			else
				V_plus *= (1 - ed.eta[edge]);
		}
		pi_P_hat = (1 - V_plus) * V_minus;
		pi_M_hat = (1 - V_minus) * V_plus;  
		float W_plus, W_minus; //, W_zero;
		if (((pi_0_hat + pi_P_hat + pi_M_hat)) == 0.0) {
			W_plus = 0.0;
			W_minus = 0.0;
		}
		else {
			W_plus = pi_P_hat / (pi_0_hat + pi_P_hat + pi_M_hat);
			W_minus = pi_M_hat / (pi_0_hat + pi_P_hat + pi_M_hat);
		}

		//W_zero = 1 - W_plus - W_minus;

		//update bias
		vars.bias[v] = fabs(W_plus - W_minus);
		vars.value[v] = (W_plus > W_minus);

		maxmag = W_plus > W_minus ? W_plus : W_minus;
#ifdef ENABLE_OPENMP
		lsummag[omp_get_thread_num()] += maxmag;
#else
		summag += maxmag;
#endif
		int ndx = my_fetch_add<int>(bias_len, 1);
		bias_list[ndx] = vars.bias[v];//update bias list
		bias_list_vars[ndx] = v;
//		printf("%d\t%.8f\n", v, maxmag);
//		if(v%1000==0||v<10) printf("summag=%.8f\n", summag);
//		if(v<10)printf("bias_list[%d]=%f\n", v, bias_list[v]);
//		if(v>3989)printf("bias_list[%d]=%f\n", v, bias_list[v]);
	}
#ifdef ENABLE_OPENMP
	for(int i=0; i<num_omp_threads; i++)
		summag += lsummag[i];
#endif
    endtime = rtclock();
#ifdef TIMING 
    printf("\truntime[update_bias] = %f ms.\n", 1000*(endtime-starttime));
#endif
}

//removes any literals whose bias is close to true or false
void decimate_2 (CSRGraph clauses, CSRGraph vars, Edge ed, int *bias_list_vars, 
		const int * bias_list_len, int fixperstep) {
	// NOTE: this is slower because the lower-work computation does not use all SMs.
	//for each vars node in the bias list
#ifdef ENABLE_OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for(int i = *bias_list_len - fixperstep; i < *bias_list_len; i++)  {
		int v = bias_list_vars[i];//get the node id
		vars.sat[v] = true;
		int edoff = vars.row_offsets[v];//get the starting index of edge id of node v
		int cllen = vars.degree(v);//get the number of outgoing edges of node v

		// for all a E V(i)
		//for each outgoing edge of node v 
		for(int edndx = 0; edndx < cllen; edndx++) {
			int edge = vars.columns[edoff + edndx];//get the edge id
			int cl = ed.src[edge];//get the src node of this edge 
			if(!clauses.sat[cl])
				if(ed.bar[edge] != vars.value[v])
					clauses.sat[cl] = true;
		}
	}
}

int compare_float(const void *x, const void *y) {
	float xx = *(float *)x, yy = *(float *)y;
	// reverse order
	if(xx > yy)
		return -1;
	else if (xx < yy)
		return 1;
	return 0;
}

float sort_bias_list(struct DoubleBuffer<float> &db_bias_list, 
		struct DoubleBuffer<int> &db_bias_list_vars, 
		int bias_list_len, float summag, int& fixperstep) {
	float r = 0;
	if(bias_list_len) {
		r = (summag / bias_list_len);
		printf("<bias>:%f\n", r);
		if(r < PARAMAGNET)
			return r;

		// Run sorting operation
		sortPairs(db_bias_list, db_bias_list_vars, bias_list_len);

		if(fixperstep > bias_list_len)
			fixperstep = 1;
	}
	return r;
}

void usage(char *argv[]) {
	fprintf(stderr, "usage: %s formula.cnf MAXLITERALS\n", argv[0]);
	exit(1);
}

//This function iterates over the clauses and the literals of the formula
//updating 'surveys' until all updates are blow EPSILON
int converge(CSRGraph &cl, CSRGraph &vars, Edge &ed) {
	float max_eps;
	int i = 0;
	do {
		calc_pi_values(cl, vars, ed);
		max_eps = 0;
		update_eta(cl, vars, ed, &max_eps);
//		printf("max_eps[%d]=%.8f\n", i, max_eps);
	} while(max_eps > EPSILON && i++ < MAXITERATION);
	if(max_eps <= EPSILON) {
		printf("converged in %d iterations max eps %f\n", i, max_eps);
		return 1;
	} else {
		printf("SP UN-CONVERGED, max eps %f\n", max_eps);
		//TODO write out formula?
		exit(1);
	}
	return 0;  
}

int build_list(CSRGraph &cl, CSRGraph &vars, Edge &ed, struct DoubleBuffer<float> &db_bias_list,
        struct DoubleBuffer<int> &db_bias_list_vars, int *bias_list_len, int &fixperstep) {
    double starttime, endtime;
	starttime = rtclock();
	float summag = 0.0f;
	*bias_list_len = 0;
//	printf("Updating bias...\n");
	update_bias(cl, vars, ed, db_bias_list.Current(), db_bias_list_vars.Current(), bias_list_len, summag);
	printf("Bias updated, bias_list_len=%d, summag=%f\n", *bias_list_len, summag);
    endtime = rtclock();
#ifdef TIMING
	printf("\truntime[update_bias] = %f ms.\n", 1000*(endtime-starttime));
#endif

//	printf("Sorting bias...\n");
	starttime = rtclock();
	float limitbias = sort_bias_list(db_bias_list, db_bias_list_vars, *bias_list_len, (float)summag, fixperstep);

//	for(int i=3990; i<4000; i++)
//		printf("bias_list[%d]=%f\n", i, db_bias_list.Current()[i]);
    endtime = rtclock();
#ifdef TIMING
	printf("\truntime[sort] = %f ms.\n", 1000*(endtime-starttime));
#endif
//	printf("Bias sorted, fixperstep=%d\n", fixperstep);
	if(limitbias < PARAMAGNET) {
		printf("paramagnetic state\n");
		return 1;
	}
	return 0;
}

int main(int argc, char *argv[]) {
	if(argc < 4)
		usage(argv);
	//srand(7);
	num_omp_threads = atoi(argv[1]);
	int max_literals = atoi(argv[3]);
	CSRGraph cl, vars;
	Edge ed;
	float *bias_list, *bias_list_2;
	int *bias_list_vars, *bias_list_vars_2;
	int bias_list_len;
	double starttime, endtime, runtime;

	init_from_file(argv[2], max_literals, cl, vars, ed);

	bias_list = (float *)malloc(vars.nnodes * sizeof(float));
	bias_list_2 = (float *)malloc(vars.nnodes * sizeof(float));
	bias_list_vars = (int *)malloc(vars.nnodes * sizeof(int));
	bias_list_vars_2 = (int *)malloc(vars.nnodes * sizeof(int));

	struct DoubleBuffer<float> db_bias_list(bias_list, bias_list_2);
	struct DoubleBuffer<int> db_bias_list_vars(bias_list_vars, bias_list_vars_2);

	int canfix = 0.01 * vars.nnodes;
	if(canfix < 1) canfix = 1;
#ifdef ENABLE_OPENMP
	omp_set_num_threads(num_omp_threads);
#endif
	starttime = rtclock();
	int round = 0;
	while(converge(cl, vars, ed)) {
		printf("round = %d\n", round++);
		if(build_list(cl, vars, ed, db_bias_list, db_bias_list_vars, &bias_list_len, canfix))
			break;
//		return 0;
		decimate_2(cl, vars, ed, db_bias_list_vars.Current(), &bias_list_len, canfix);

	};
	endtime = rtclock();
	runtime = (1000.0 * (endtime - starttime));
	printf("\truntime [nsp] = %f ms.\n", runtime);
	print_solution("sp_sol.dat", vars);
	dump_partial(PARTIAL, cl, vars, ed);
	return 0;
}

