/** Delaunay refinement -*- C++ -*-
 * Copyright (C) 2013, Xuhao Chen. All rights reserved.
 * Refinement of an initial, unrefined Delaunay mesh to eliminate triangles
 * with angles < 30 degrees
 * @author: Xuhao Chen <cxh.nudt@gmail.com>
 */

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include "common.h"
#include "meshfiles.h"
#include "dmr.h"
#include "geomprim.h"
#include "worklistc.h"
#include "devel.h"
#include <map>

#define CAVLEN 256
#define BCLEN 1024
#define debug 0

int num_omp_threads;
void check_is_bad(Mesh &mesh, int ele) {
	uint3 *el = &mesh.elements[ele];
	mesh.isbad[ele] = (angleLT(mesh, el->x, el->y, el->z)
			|| angleLT(mesh, el->z, el->x, el->y) 
			|| angleLT(mesh, el->y, el->z, el->x));
}

bool shares_edge(uint nodes1[3], uint nodes2[3]) {
	int i;
	int match = 0;
	uint help;

	for (i = 0; i < 3; i++) {
		if ((help = nodes1[i]) != INVALIDID) {
			if (help == nodes2[0]) match++;
			else if (help == nodes2[1]) match++;
			else if (help == nodes2[2]) match++;
		}
	} 
	return match == 2;
}
/*
void find_neighbours(Mesh mesh, int start, int end) {
	//int id = threadIdx.x + blockDim.x * blockIdx.x;
	//int threads = blockDim.x * gridDim.x;
	int ele;
	int oele;
	int nc = 0;
	uint nodes1[3], nodes2[3], neigh[3] = {INVALIDID, INVALIDID, INVALIDID};

	for(int x = 0; x < mesh.nelements; x ++) {
		// currently a n^2 algorithm -- launch times out for 250k.ele!
		for(ele = start; ele < end; ele ++) {
			if(x == 0) {
				neigh[0] = INVALIDID;
				neigh[1] = INVALIDID;
				neigh[2] = INVALIDID;
			}
			else {
				neigh[0] = mesh.neighbours[ele].x;
				neigh[1] = mesh.neighbours[ele].y;
				neigh[2] = mesh.neighbours[ele].z;
			}
			if(neigh[2] != INVALIDID) continue;
			//TODO: possibly remove uint3 from Mesh/ShMesh
			nodes1[0] = mesh.elements[ele].x;
			nodes1[1] = mesh.elements[ele].y;
			nodes1[2] = mesh.elements[ele].z;
			nc = (neigh[0] == INVALIDID) ? 0 : ((neigh[1] == INVALIDID) ? 1 : 2);

			//TODO: block this
			for(oele=0; oele < mesh.nelements; oele++) {
				nodes2[0] = mesh.elements[oele].x; 
				nodes2[1] = mesh.elements[oele].y; 
				nodes2[2] = mesh.elements[oele].z;
				if(shares_edge(nodes1, nodes2)) {
					assert(nc < 3);
					neigh[nc++] = oele;
				}
				if((IS_SEGMENT(mesh.elements[ele]) && nc == 2) || nc == 3)
					break;
			}

			mesh.neighbours[ele].x = neigh[0];
			mesh.neighbours[ele].y = neigh[1];
			mesh.neighbours[ele].z = neigh[2];
		}
	}
}
*/
void dump_mesh_element(Mesh &mesh, uint3 &ele, int element) {
	if(IS_SEGMENT(ele))
		printf("[ %.17f %.17f %.17f %.17f %d]\n", 
				mesh.nodex[ele.x], mesh.nodey[ele.x],
				mesh.nodex[ele.y], mesh.nodey[ele.y], element);
	else
		printf("[ %.17f %.17f %.17f %.17f %.17f %.17f %d]\n", 
				mesh.nodex[ele.x], mesh.nodey[ele.x],
				mesh.nodex[ele.y], mesh.nodey[ele.y],
				mesh.nodex[ele.z], mesh.nodey[ele.z], element);
}

void dump_mesh_element(Mesh &mesh, int element) {
	dump_mesh_element(mesh, mesh.elements[element], element);
}

bool encroached(Mesh &mesh, int element, uint3 &celement, FORD centerx, FORD centery, bool &is_seg) {
	if(element == INVALIDID)
		return false;
	assert(!mesh.isdel[element]);
	uint3 ele = mesh.elements[element];
	if(IS_SEGMENT(ele)) {
		FORD cx, cy, radsqr;
		uint nsp;
		is_seg = true;
		nsp = (celement.x == ele.x) ? ((celement.y == ele.y) ? celement.z : celement.y) : celement.x;

		// check if center and triangle are on opposite sides of segment
		// one of the ccws does not return zero
		if(counterclockwise(mesh.nodex[ele.x], 
			mesh.nodey[ele.x], mesh.nodex[ele.y], 
			mesh.nodey[ele.y], mesh.nodex[nsp], 
			mesh.nodey[nsp]) > 0 != 
			counterclockwise(mesh.nodex[ele.x], 
			mesh.nodey[ele.x], mesh.nodex[ele.y], 
			mesh.nodey[ele.y], centerx, centery) > 0)
			return true; 

		// nope, do a distance check
		cx = (mesh.nodex[ele.x] + mesh.nodex[ele.y]) / 2;
		cy = (mesh.nodey[ele.x] + mesh.nodey[ele.y]) / 2;
		radsqr = distanceSquare(cx, cy, mesh.nodex[ele.x], mesh.nodey[ele.x]);

		return distanceSquare(centerx, centery, cx, cy) < radsqr;
	} else
		return gincircle(mesh.nodex[ele.x], mesh.nodey[ele.x],
				mesh.nodex[ele.y], mesh.nodey[ele.y],
				mesh.nodex[ele.z], mesh.nodey[ele.z],
				centerx, centery) > 0.0;
}

void add_to_cavity(uint cavity[], uint &cavlen, int element) {
	int i;
	for(i = 0; i < cavlen; i++)
		if(cavity[i] == element)
			return;
	cavity[cavlen++] = element;
}

void add_to_boundary(uint boundary[], uint &boundarylen, uint sn1, uint sn2, uint src, uint dst) {
	int i;
	for(i = 0; i < boundarylen; i+=4)
		if((sn1 == boundary[i] && sn2 == boundary[i+1]) ||
				(sn1 == boundary[i+1] && sn2 == boundary[i]))
			return;
	boundary[boundarylen++] = sn1;  
	boundary[boundarylen++] = sn2;
	boundary[boundarylen++] = src;
	boundary[boundarylen++] = dst;
}

unsigned add_node(Mesh &mesh, FORD x, FORD y, uint ndx) {
	assert(ndx < mesh.maxnnodes);
//	printf("ndx=%d, maxnnodes=%d\n", ndx, mesh.maxnnodes);
	mesh.nodex[ndx] = x;
	mesh.nodey[ndx] = y;  
	return ndx;
}

uint add_segment(Mesh &mesh, uint n1, uint n2, uint ndx) {
	//TODO: parallelize
	uint3 ele;
	ele.x = n1; ele.y = n2; ele.z = INVALIDID;
	assert(ndx < mesh.maxnelements);
	mesh.isbad[ndx] = false;
	mesh.isdel[ndx] = false;
	mesh.elements[ndx] = ele;
	mesh.neighbours[ndx].x = mesh.neighbours[ndx].y = mesh.neighbours[ndx].z = INVALIDID;
	mesh.owners[ndx] = INIT_OWNER;//cxh
	return ndx;
}

uint add_triangle(Mesh &mesh, uint n1, uint n2, uint n3, uint nb1, uint oldt, uint ndx) {
	uint3 ele;
	if(counterclockwise(mesh.nodex[n1], mesh.nodey[n1], 
				mesh.nodex[n2], mesh.nodey[n2],
				mesh.nodex[n3], mesh.nodey[n3]) > 0) {
		ele.x = n1; ele.y = n2; ele.z = n3;
	}
	else {
		ele.x = n3; ele.y = n2; ele.z = n1;
	}
	assert(ndx < mesh.maxnelements);
	mesh.isbad[ndx] = false;
	mesh.isdel[ndx] = false;
	mesh.elements[ndx] = ele;
	mesh.neighbours[ndx].x = nb1;
	mesh.neighbours[ndx].y = mesh.neighbours[ndx].z = INVALIDID;
	mesh.owners[ndx] = INIT_OWNER;//cxh
	uint3 *nb = &mesh.neighbours[nb1];
	if(mesh.neighbours[nb1].x == oldt)
		nb->x = ndx;
	else {
		if(mesh.neighbours[nb1].y == oldt)
			nb->y = ndx;
		else {
			if(mesh.neighbours[nb1].z != oldt)
				printf("%u %u %u %u %u %u\n", ndx, oldt, nb1, mesh.neighbours[nb1].x, 
						mesh.neighbours[nb1].y, mesh.neighbours[nb1].z);
			assert(mesh.neighbours[nb1].z == oldt);
			nb->z = ndx;
		}
	}
	return ndx;
}

bool adjacent(uint3 &elem1, uint3 &elem2) {
	int sc = 0;
	if(elem1.x == elem2.x || elem1.x == elem2.y || elem1.x == elem2.z)
		sc++;
	if(elem1.y == elem2.x || elem1.y == elem2.y || elem1.y == elem2.z)
		sc++;
	if(!IS_SEGMENT(elem1) && (elem1.z == elem2.x || elem1.z == elem2.y || elem1.z == elem2.z))
		sc++;
	return sc == 2;
}

void find_shared_edge(uint3 &elem1, uint3 &elem2, uint se[2]) {
	int sc = 0;
	if(elem1.x == elem2.x || elem1.x == elem2.y || elem1.x == elem2.z)
		se[sc++] = elem1.x;

	if(elem1.y == elem2.x || elem1.y == elem2.y || elem1.y == elem2.z)
		se[sc++] = elem1.y;

	if(!IS_SEGMENT(elem1) && (elem1.z == elem2.x || elem1.z == elem2.y || elem1.z == elem2.z))
		se[sc++] = elem1.z;
	assert(sc == 2);
	assert(se[0] != INVALIDID);
	assert(se[1] != INVALIDID);
}

bool build_cavity(Mesh &mesh, uint cavity[], uint &cavlen, int max_cavity, uint boundary[], uint &boundarylen, FORD &cx, FORD &cy) {
	int ce = 0;
	uint3 ele = mesh.elements[cavity[0]];
	bool is_seg = false;
	if(IS_SEGMENT(ele)) {
		cx = (mesh.nodex[ele.x] + mesh.nodex[ele.y]) / 2;
		cy = (mesh.nodey[ele.x] + mesh.nodey[ele.y]) / 2;
	}
	else {
		circumcenter(mesh.nodex[ele.x], mesh.nodey[ele.x],
				mesh.nodex[ele.y], mesh.nodey[ele.y],
				mesh.nodex[ele.z], mesh.nodey[ele.z],
				cx, cy);
	}
	//printf("highlight %d %d [%f %f]\n", cavity[0], IS_SEGMENT(ele), cx, cy);
	while (ce < cavlen) {
		if(mesh.isdel[cavity[ce]])
			printf("deleted: %d\n", cavity[ce]);
		assert(cavlen < max_cavity);
		assert(!mesh.isdel[cavity[ce]]);
		uint3 neighbours = mesh.neighbours[cavity[ce]];
		uint neighb[3] = {neighbours.x, neighbours.y, neighbours.z};
		for(int i = 0; i < 3; i++) {
			if(neighb[i] == cavity[0])
				continue;
			if(neighb[i] == INVALIDID)
				continue;
			//printf("neigbour %d\n", neighb[i]);
			is_seg  = false;
			if(!(IS_SEGMENT(ele) && IS_SEGMENT(mesh.elements[neighb[i]])) && 
					encroached(mesh, neighb[i], ele, cx, cy, is_seg)) {
				if(!is_seg)
					add_to_cavity(cavity, cavlen, neighb[i]);
				else {
					assert(!IS_SEGMENT(ele));
					cavity[0] = neighb[i];
					cavlen = 1;
					boundarylen = 0;
					return false;
				}
			} else {
				uint se[2];
				if(!adjacent(mesh.elements[cavity[ce]], mesh.elements[neighb[i]])) {
					dump_mesh_element(mesh, cavity[ce]);
					dump_mesh_element(mesh, neighb[i]);
					printf("%d %d\n", cavity[ce], neighb[i]);
				}
				assert(boundarylen < BCLEN);
				find_shared_edge(mesh.elements[cavity[ce]], mesh.elements[neighb[i]], se);
				add_to_boundary(boundary, boundarylen, se[0], se[1], neighb[i], cavity[ce]);
			}
		}
		ce++;
	}
	return true;
}

void addneighbour(Mesh &mesh, uint3 &neigh, uint elem) {
	// TODO
	if(neigh.x == elem || neigh.y == elem || neigh.z == elem) return;
	assert(neigh.x == INVALIDID || neigh.y == INVALIDID || neigh.z == INVALIDID);
	if(neigh.x == INVALIDID) { neigh.x = elem; return; }
	if(neigh.y == INVALIDID) { neigh.y = elem; return; }
	if(neigh.z == INVALIDID) { neigh.z = elem; return; }
}

void setup_neighbours(Mesh &mesh, uint start, uint end) {
	// relies on all neighbours being in start--end
	for(uint i = start; i < end; i++) {
		uint3 &neigh = mesh.neighbours[i];
		for(uint j = i+1; j < end; j++) {
			if(adjacent(mesh.elements[i], mesh.elements[j])) {
				addneighbour(mesh, neigh, j);
				addneighbour(mesh, mesh.neighbours[j], i);
			}
		}    
	}
}

uint opposite(Mesh &mesh, uint element) {
	bool obtuse = false;
	int obNode = INVALIDID;
	uint3 el = mesh.elements[element];
	if(IS_SEGMENT(el))
		return element;
	// figure out obtuse node
	if(angleOB(mesh, el.x, el.y, el.z)) {
		obtuse = true;
		obNode = el.z;
	} else {
		if(angleOB(mesh, el.z, el.x, el.y)) {
			obtuse = true;
			obNode = el.y;
		} else {
			if(angleOB(mesh, el.y, el.z, el.x)) {
				obtuse = true;
				obNode = el.x;
			}
		}
	}
	if(obtuse) {
		// find the neighbour that shares an edge whose points do not include obNode
		uint se_nodes[2];
		uint nobneigh;
		uint3 neigh = mesh.neighbours[element];
		//printf("obtuse node [%f %f]\n", mesh.nodex[obNode], mesh.nodey[obNode]);
		assert(neigh.x != INVALIDID && neigh.y != INVALIDID && neigh.z != INVALIDID);
		nobneigh = neigh.x;
		find_shared_edge(el, mesh.elements[neigh.x], se_nodes);
		if(se_nodes[0] == obNode || se_nodes[1] == obNode) {
			nobneigh = neigh.y;
			find_shared_edge(el, mesh.elements[neigh.y], se_nodes);
			if(se_nodes[0] == obNode || se_nodes[1] == obNode) {
				nobneigh = neigh.z;
			}
		}
		return nobneigh;
	}
	return element;
}

void refine(Mesh &mesh, uint *nnodes, uint *nelements,
		Worklist2 wl, Worklist2 &owl) {
#ifdef ENABLE_OPENMP
	uint ulimit = ((wl.nitems()+num_omp_threads-1)/num_omp_threads)*num_omp_threads;
#else
	uint ulimit = wl.nitems();
#endif
	//for each triangle 'ele' in my worklist
#ifdef ENABLE_OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for(int eleit = 0; eleit < ulimit; eleit++) {
		uint cavity[CAVLEN];
		uint boundary[BCLEN];
		int ele = 0;
		uint nc = 0, bc = 0;
		bool repush = false;
		int id = eleit;
		FORD cx, cy;
		int elems_added = 0;
#ifdef ENABLE_OPENMP
		int stage = 0;
		id = omp_get_thread_num();
#endif
//		printf("Iteration: %d, thread_id: %d\n", eleit, id);
		int haselem = wl.pop_id(eleit, ele);
		if(debug) printf("%d:%d:%d\n", eleit, ele, haselem);
		// if 'ele' is bad and it is not deleted
		if(haselem && ele < mesh.nelements && mesh.isbad[ele] && !mesh.isdel[ele]) {
			// create and expand cavity around 'ele'
			cavity[nc++] = ele;
			if(debug) {
				printf("original center element ");
				dump_mesh_element(mesh, cavity[0]);
			}
			uint oldcav;
			do {
				oldcav = cavity[0];
				cavity[0] = opposite(mesh, ele);
			} while(cavity[0] != oldcav);
			if(!build_cavity(mesh, cavity, nc, CAVLEN, boundary, bc, cx, cy))
				build_cavity(mesh, cavity, nc, CAVLEN, boundary, bc, cx, cy);
			if(debug) {
				printf("center element [%f %f] ", cx, cy);
				dump_mesh_element(mesh, cavity[0]);
				printf("pre-graph %d\n", nc);
				for(int i = 1; i < nc; i++) {
					dump_mesh_element(mesh, cavity[i]);
				}
				printf("boundary-edges %d\n", bc);
				for(int i = 0; i < bc; i+=4) {
					printf("[%f %f %f %f]\n", 
						mesh.nodex[boundary[i]], mesh.nodey[boundary[i]],
						mesh.nodex[boundary[i+1]], mesh.nodey[boundary[i+1]]);
					dump_mesh_element(mesh, boundary[i+2]);
				}
			}
///*
#ifdef ENABLE_OPENMP
			// try to claim ownership
			// mark all triangles in the cavity with my (thread) id
			for(int i = 0; i < nc; i++)
				mesh.owners[cavity[i]] = id;

			for(int i = 0; i < bc; i+=4)
				mesh.owners[boundary[i + 2]] = id;
			stage = 1;
		}
		__syncthreads();
		if(stage == 1) {
			// check and re-mark with priority
			// check for conflicts (TODO: atomic)
			for(int i = 0; i < nc; i++) {
				if(mesh.owners[cavity[i]] > id)
					#pragma omp critical
					{
					if(mesh.owners[cavity[i]] > id)
						mesh.owners[cavity[i]] = id;
					}
//					mesh.owners[cavity[i]] = id;
			}
			for(int i = 0; i < bc; i+=4) {
				if(mesh.owners[boundary[i + 2]] > id)
					#pragma omp critical
					{
					if(mesh.owners[boundary[i + 2]] > id)
						mesh.owners[boundary[i + 2]] = id;
					}
//					mesh.owners[boundary[i + 2]] = id;
			}
			stage = 2;
		}
//		printf("id=%d: end stage 1\n", id);
		__syncthreads();
//		printf("id=%d: start stage 2\n", id);
		//In stage 2, if all cavity-triangles are marked with my thread id, goto stage 3, otherwise back-off
		if(stage == 2) {
			for(int i = 0; i < nc; i++) {
				if(mesh.owners[cavity[i]] != id) {
					repush = true;
					//printf("%d conflict\n", ele);
					//printf("%d: %d owned by %d\n", 
					//	id, cavity[i], mesh.owners[cavity[i]]);
					break;
				}
			}
			if(!repush)
				for(int i = 0; i < bc; i+=4) {
					if(mesh.owners[boundary[i + 2]] != id) {
						repush = true;
						//printf("%d conflict\n", ele);
						//printf("%d: %d owned by %d\n",
						//	id, boundary[i + 2], mesh.owners[boundary[i + 2]]);
						break;
					}
				}
			if(!repush) {
				stage = 3;
#endif
				elems_added = (bc>>2) + (IS_SEGMENT(mesh.elements[cavity[0]]) ? 2:0);
		
#ifdef ENABLE_OPENMP
			}
		}
		if(stage == 3) {
#endif
//*/
/*
#ifdef ENABLE_OPENMP
			#pragma omp critical
			{
			for(int i = 0; i < nc; i++) {
				if(mesh.owners[cavity[i]] != INIT_OWNER) {
					repush = true;
					break;
				}
			}
			if(!repush) {
				for(int i = 0; i < bc; i+=4) {
					if(mesh.owners[boundary[i+2]] != INIT_OWNER) {
						repush = true;
						break;
					}
				}
			}
			if(!repush) {
				for(int i = 0; i < nc; i++)
					mesh.owners[cavity[i]] = id;
				for(int i = 0; i < bc; i+=4)
					mesh.owners[boundary[i+2]] = id;
				stage = 3;
#endif
				elems_added = (bc>>2) + (IS_SEGMENT(mesh.elements[cavity[0]]) ? 2:0);
#ifdef ENABLE_OPENMP
			}
			}
		}
		if(stage == 3) {
#endif
//*/
			//create new cavity by retriangulating
			uint old_nnodes = my_fetch_add<uint>(nnodes,1);
			uint cnode = add_node(mesh, cx, cy, old_nnodes);
			uint cseg1 = 0, cseg2 = 0;
			uint nelements_added = elems_added;
			//printf("start: %d %d %d %d %d\n", 
			//	id, elems_added, start, offset, start+offset);
			uint oldelements = my_fetch_add<uint>(nelements, nelements_added);
			uint newelemndx = oldelements;
			if(debug) printf("post-graph\n");
			if(IS_SEGMENT(mesh.elements[cavity[0]])) {
				cseg1 = add_segment(mesh, mesh.elements[cavity[0]].x, cnode, newelemndx++);
				cseg2 = add_segment(mesh, cnode, mesh.elements[cavity[0]].y, newelemndx++);
				if(debug) {
					dump_mesh_element(mesh, cseg1);
					dump_mesh_element(mesh, cseg2);
				}
			}
			// add triangles from new cavity to mesh
			for(int i = 0; i < bc; i+=4) {
				uint ntri  = add_triangle(mesh, boundary[i], boundary[i+1], 
						cnode, boundary[i+2], boundary[i+3], newelemndx++);
				if(debug) dump_mesh_element(mesh, ntri);
			}
			assert(oldelements + nelements_added == newelemndx);
			setup_neighbours(mesh, oldelements, newelemndx);
			repush = true;
			//delete reiangles in old cavity from mesh
			for(int i = 0; i < nc; i++) {
				mesh.isdel[cavity[i]] = true;
				if(cavity[i] == ele) repush = false;  
				//printf("%d: deleting %d\n", id, cavity[i]);
			}
			if(debug) printf("update over\n");
		}
		// push it to the worklist
		if(repush) owl.push(ele);
#ifdef ENABLE_OPENMP
		__syncthreads();
#endif
	}
}

void check_triangles(Mesh mesh, uint &bad_triangles, Worklist2 *wl, int start) {
#ifdef ENABLE_OPENMP
	int *count = (int *) malloc(num_omp_threads*sizeof(int));
	for(int i=0; i<num_omp_threads; i++)
		count[i] = 0;
#else
	int count = 0;
#endif
//	printf("start %d, nelements %d\n", start, mesh.nelements);
	#pragma omp parallel for schedule(static)
	for(int ele = start; ele < mesh.nelements; ele ++) {
		bool push = false;
		if(mesh.isdel[ele])
			continue;
		if(IS_SEGMENT(mesh.elements[ele]))
			continue;
		if(!mesh.isbad[ele]) {
			uint3 *el = &mesh.elements[ele];
			mesh.isbad[ele] = (angleLT(mesh, el->x, el->y, el->z) 
					|| angleLT(mesh, el->z, el->x, el->y) 
					|| angleLT(mesh, el->y, el->z, el->x));
		}
		if(mesh.isbad[ele]) {
			push = true;
#ifdef ENABLE_OPENMP
			count[omp_get_thread_num()] += 1;
#else
			count++;
#endif
		}
		if(push) wl->push(ele); // not really as bad as it looks
	}
#ifdef ENABLE_OPENMP
	for(int i=0; i<num_omp_threads; i++)
		bad_triangles += count[i];
#else
	bad_triangles = count;
#endif
}

void verify_mesh(Mesh &mesh) {
	// TODO: check for delaunay property
}

void addneighbour_cpu(uint3 &neigh, uint elem) {
	// TODO
	if(neigh.x == elem || neigh.y == elem || neigh.z == elem) return;
	assert(neigh.x == INVALIDID || neigh.y == INVALIDID || neigh.z == INVALIDID);
	if(neigh.x == INVALIDID) { neigh.x = elem; return; }
	if(neigh.y == INVALIDID) { neigh.y = elem; return; }
	if(neigh.z == INVALIDID) { neigh.z = elem; return; }
}

void find_neighbours_cpu(Mesh &mesh) {
	std::map<std::pair<int, int>, int> edge_map;
	uint nodes1[3];
	uint3 *elements = mesh.elements;
	uint3 *neighbours = mesh.neighbours;
	int ele;
	for(ele = 0; ele < mesh.nelements; ele++) {
		uint3 *neigh = &neighbours[ele];
		neigh->x = INVALIDID;
		neigh->y = INVALIDID;
		neigh->z = INVALIDID;
		nodes1[0] = elements[ele].x;
		nodes1[1] = elements[ele].y;
		nodes1[2] = elements[ele].z;
		if(nodes1[0] > nodes1[1]) std::swap(nodes1[0], nodes1[1]);
		if(nodes1[1] > nodes1[2]) std::swap(nodes1[1], nodes1[2]);
		if(nodes1[0] > nodes1[1]) std::swap(nodes1[0], nodes1[1]);
		assert(nodes1[0] <= nodes1[1] && nodes1[1] <= nodes1[2]);
		std::pair<int, int> edges[3];
		edges[0] = std::make_pair<int, int>(nodes1[0], nodes1[1]);
		edges[1] = std::make_pair<int, int>(nodes1[1], nodes1[2]);
		edges[2] = std::make_pair<int, int>(nodes1[0], nodes1[2]);
		int maxn = IS_SEGMENT(elements[ele]) ? 1 : 3;
		for(int i = 0; i < maxn; i++) {
			if(edge_map.find(edges[i]) == edge_map.end())
				edge_map[edges[i]] = ele;
			else {	  
				int node = edge_map[edges[i]];
				addneighbour_cpu(neighbours[node], ele);
				addneighbour_cpu(neighbours[ele], node);
				edge_map.erase(edges[i]);
			}
		}
	}
}

void refine_mesh(Mesh &mesh) {
	uint nbad, nelements, nnodes;
	Worklist2 wl1(mesh.nelements), wl2(mesh.nelements);
	find_neighbours_cpu(mesh);
	//dump_neighbours(mesh);
	nelements = mesh.nelements;
	nnodes = mesh.nnodes;

	double starttime, endtime;
	int lastnelements = 0;
	Worklist2 *inwl = &wl1, *outwl = &wl2;
	starttime = rtclock();
	nbad = 0;
//	printf("checking triangles...\n");
	check_triangles(mesh, nbad, inwl, 0); 
	printf("%d initial bad triangles\n", nbad);
	while(nbad) {   
//		if(debug) inwl->display_items();
		lastnelements = mesh.nelements;
		printf("in elements: %d\n", inwl->nitems());
		refine(mesh, &nnodes, &nelements, *inwl, *outwl);
		printf("refine over\n");
		mesh.nnodes = nnodes;
		mesh.nelements = nelements;
		printf("out elements: %d\n", outwl->nitems());
		std::swap(inwl, outwl);
		outwl->reset();
		nbad = 0;
//		printf("checking triangles...\n");
		// need to check only new triangles
		check_triangles(mesh, nbad, inwl, lastnelements); 
		nbad = inwl->nitems();
		printf("%u bad triangles.\n", nbad);
		if(0)
			debug_isbad(*inwl, mesh);
		if(nbad == 0) 
			break;
//		return;
	}
	endtime = rtclock();
	nbad = 0;
	check_triangles(mesh, nbad, inwl, 0); 
	printf("%d final bad triangles\n", nbad);
	printf("time: %f ms\n", (endtime-starttime)*1000);
}

void read_mesh(const char *basefile, Mesh &mesh, int maxfactor) {
	readNodes(basefile, mesh, maxfactor);
	readTriangles(basefile, mesh, maxfactor);
	assert(mesh.maxnelements > 0);
	printf("Number of elements: %d\n", mesh.nelements);
	//printf("memory for owners: %d MB\n", mesh.maxnelements * sizeof(int) / 1048576);
	//mesh.owners.alloc(mesh.maxnelements);
	mesh.init();
	printf("memory for worklists: %d MB\n", 2 * mesh.nelements * sizeof(int) / 1048576);

	printf("%s: %d nodes, %d triangles, %d segments read\n", basefile, mesh.nnodes, mesh.ntriangles, mesh.nsegments);
	assert(mesh.nnodes > 0);
	assert(mesh.ntriangles > 0);
	assert(mesh.nsegments > 0);
	assert(mesh.nelements > 0);
}

int main(int argc, char *argv[]) {
	Mesh mesh;
	int maxfactor = 2;
	int mesh_nodes, mesh_elements;
	if(argc < 3) {
		printf("Usage: %s numThreads basefile [maxfactor]\n", argv[0]);
		exit(0);
	}
	if(argc == 4) {
		maxfactor = atoi(argv[3]);
	}
	num_omp_threads = atoi(argv[1]);
	printf("numThreads: %d, InputFile: %s, maxfactor: %d\n", num_omp_threads, argv[2], maxfactor);
	read_mesh(argv[2], mesh, maxfactor);
	mesh_nodes = mesh.nnodes;
	mesh_elements = mesh.ntriangles + mesh.nsegments;
#ifdef ENABLE_OPENMP
	omp_set_num_threads(num_omp_threads);
#endif
	refine_mesh(mesh);
	printf("%f increase in number of elements (maxfactor hint)\n", 1.0 * mesh.nelements / mesh_elements);
	printf("%f increase in number of nodes (maxfactor hint)\n", 1.0 * mesh.nnodes / mesh_nodes);
	verify_mesh(mesh);
	write_mesh(argv[2], mesh);
	return 0;
}
