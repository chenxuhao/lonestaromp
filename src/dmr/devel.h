#pragma once

void dump_element(FORD *nodex, FORD *nodey, uint3 &element, int ele) {
	printf("%f %f %f %f", nodex[element.x], nodey[element.x],
			nodex[element.y], nodey[element.y]);
	if(!IS_SEGMENT(element))
		printf(" %f %f", nodex[element.z], nodey[element.z]);
	printf(" %d", ele);
}

void dump_neighbours(Mesh &mesh) {
	FORD *nodex = mesh.nodex;
	FORD *nodey = mesh.nodey;
	uint3 *elements = mesh.elements;
	uint3 *neighbours = mesh.neighbours;

	for(int i = 0; i < mesh.nelements; i++) {
		printf("center element [");
		dump_element(nodex, nodey, elements[i], i);
		printf("]\npre-graph\n");
		if(neighbours[i].x != INVALIDID) {
			printf("["); dump_element(nodex, nodey, elements[neighbours[i].x], i); printf("]\n");
		}
		if(neighbours[i].y != INVALIDID) {
			printf("["); dump_element(nodex, nodey, elements[neighbours[i].y], i); printf("]\n");
		}
		if(neighbours[i].z != INVALIDID) {
			printf("["); dump_element(nodex, nodey, elements[neighbours[i].z], i); printf("]\n");
		}	
		printf("post-graph\n");
		printf("update over\n");      
	}
}

void debug_isbad(Worklist2 &wl, Mesh &mesh) {
	bool *isbad = mesh.isbad;
	volatile bool *isdel = mesh.isdel;
	int i;
	int badcount = 0;

	int wlitems = wl.length;

	printf("checking %d elements\n", mesh.nelements);
	for(i = 0; i < mesh.nelements; i++) {
		if(isdel[i])
			continue;
		if(isbad[i]) {
			badcount++;
			int j = 0;
			for(j = 0; j < wlitems; j++) {
				if(wl.wl[j] == i)
					break;
			}
			if(j == wlitems)
				printf("\tnot found: %d\n", i);	
		}
	}
	printf("bad count: %d\n", badcount);  
}
