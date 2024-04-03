/*
---------------------------------------------------------------------
 This file is part of BADIOS framework
 Copyright (c) 2012,
 By:    Ahmet Erdem Sariyuce,
        Erik Saule,
        Kamer Kaya,
        Umit V. Catalyurek
---------------------------------------------------------------------
 For license info, please see the README.txt and LICENSE.txt files in
 the main directory.
---------------------------------------------------------------------
*/

#include "bc-seq-brandes.h"



void extract_cc_info (int nVtx, int* xadj, int* adj, int* numof_CCs, int* maxCC_nvtx,
		int* maxCC_nedge, int* total_nvtx, int* total_nedge) {
	int* mark = (int*) calloc (nVtx, sizeof(int));
	int* component = (int*) malloc (sizeof(int) * nVtx);

	int comp_id = 0;
	for (vertex i = 0; i < nVtx; i++) {
		if(mark[i] == 0) {
			memset(component, -1, sizeof(vertex) * nVtx);
			// detecting connected component
			int compcounter = 0;
			component[compcounter++] = i;
			mark[i] = 1;
			int cur = compcounter - 1;
			while (cur != compcounter) {
				vertex v = component[cur];
				for (myindex j = xadj[v]; j < xadj[v+1]; j++) {
					vertex w = adj[j];
					if ((w != -1) && (mark[w] == 0)) {
						mark[w] = 1;
						component[compcounter] = w;
						compcounter++;
					}
				}
				cur++;
			}
			// count nedges of comp
			int nedges = 0;
			for (int i = 0; i < compcounter; i++) {
				int u = component[i];
				for (myindex j = xadj[u]; j < xadj[u+1]; j++) {
					if (adj[j] != -1)
						nedges++;
				}
			}
			nedges /= 2;

			if (compcounter > 1) {
				*total_nvtx += compcounter;
				if (compcounter > *maxCC_nvtx)
					*maxCC_nvtx = compcounter;

				*total_nedge += nedges;
				if (nedges > *maxCC_nedge)
					*maxCC_nedge = nedges;

				printf("COMP: %d, NVTX: %d and NEDGE: %d\n", comp_id++, compcounter, nedges);
			}
		}
	}
	*numof_CCs = comp_id;
}


void* myre_alloc (void* ptr, size_t size) {

	void* retval = realloc(ptr, size);
	return retval;
}

/**
 * 從 v 探出去的所有 node 跟 v 一樣都是同個 component id (cid)
 * @brief
 * 這裡也有改成，在探訪AP本尊所連接的每一個不同的component時，也有把每個component的ff累加
*/
void assign_component_ids (int cid, vertex v, vertex u, int* componentid, vertex* adj, int* xadj, int nVtx, double* comp_dist_from_u, int* dist_arr, double* weight, double* ff) {
	
	int* bfsorder = new int[nVtx];
	int* mark = new int[nVtx];

	#pragma region CC
	
	// for(int i = 121421 ; i < nVtx ; i ++){
	// 	printf("ff[%d] = %f\n", i, ff[i]);
	// }

	memset(dist_arr, -1, sizeof(int) * nVtx);
	for(myindex j = xadj[u] ; j < xadj[u + 1] ; j++){
		vertex y = adj[j];
		if(y >= nVtx){
			// printf("AP_ID : %d, neighborID = %d, nVtx = %d\n", u, y, nVtx);
			exit(1);
		}
		if(y != -1){
			dist_arr[y] = 1;
			// printf("dist[%d] = %f\n", y, dist_arr[y]);
		}
	}

	#pragma endregion //CC	

	memset (bfsorder, 0, sizeof(int) * nVtx);
	memset (mark, 0, sizeof(int) * nVtx);
	int endofbfsorder = 1;
	bfsorder[0] = v;
	int cur = 0;
	mark[v] = 1;
	
	// printf("\tcid = %d\n", cid);
	while (cur != endofbfsorder) {
		vertex x = bfsorder[cur];
		componentid[x] = cid;
		
		#pragma region CC
		comp_dist_from_u[cid + 1] += ((double)dist_arr[x]) * weight[x] + ff[x];
		// if(std::isnan(comp_dist_from_u[cid + 1])){
		// 	printf("currentNodeID = %d, dist = %d, weight = %f, ff = %f, comp_dist_from_u[%d] = NAN\n", x, dist_arr[x], weight[x], ff[x], cid + 1);
		// 	exit(1);
		// }
		// printf("comp_dist_from_u[%d] = %f\n", cid + 1, comp_dist_from_u[cid + 1]);
		#pragma endregion //CC

		for (myindex j = xadj[x]; j < xadj[x + 1]; j++) {
			vertex y = adj[j];
			
			if ((y != -1) && (mark[y] == 0) && (y != u)) {
				//如果探訪到的 node y 是 u 的話就跳過，這樣就不會探訪到跟 u 有關的所有 node 了
				bfsorder[endofbfsorder++] = y;
				mark[y] = 1;
				
				#pragma region CC
				if(dist_arr[y] == -1){
					dist_arr[y] = dist_arr[x] + 1;
				}
				#pragma endregion //CC
			}
		}
		cur++;
	}
	if((cid + 1) == 121446){
		printf("comp_dist_from_u[%d] = %f\n", cid + 1, comp_dist_from_u[cid + 1]);
	}
	// exit(1);
	delete[] bfsorder;
	delete[] mark;
	
	#pragma region CC
	// free(dist_arr);
	// printf("[assign_component_ids] : free(dist)\n");
	// exit(1);
	#pragma endregion //CC
}


int compare(const void *p1, const void *p2) {
	/* The actual arguments to this function are "pointers to
pointers to char", but strcmp(3) arguments are "pointers
to char", hence the following cast plus dereference */
	int a = *((int *)p1);
	int b = *((int *)p2);
	return (a > b)?1:0;
}



/**
 * @brief 以 u 當起點對 u 所在的 component 進行 BFS 取得該 component 的 weight
 * 
 * @return the weight of the component that u belongs to
*/
double totalw_idv(vertex u, int nVtx, vertex* adj, int* xadj, double* weight, int* idv_track, idv_info** identical_sets,
		int* mark, int* bfsorder) {
	double total_weight = 0;
	memset (mark, 0, sizeof(int) * nVtx);
	memset (bfsorder, 0, sizeof(int) * nVtx);
	int endofbfsorder = 1;
	bfsorder[0] = u;
	int cur = 0;
	mark[u] = 1;
	while (cur != endofbfsorder) {
		vertex v = bfsorder[cur];
		if (idv_track[v] == -1)
			total_weight += weight[v];
		else {
			int idx_of_v = idv_track[v];
			total_weight += identical_sets[idx_of_v][0].weight;
		}
		for (myindex j = xadj[v]; j < xadj[v + 1]; j++) {
			vertex w = adj[j];
			if ((w != -1) && (mark[w] == 0)) {
				bfsorder[endofbfsorder++] = w;
				mark[w] = 1;
			}
		}
		cur++;
	}

	return total_weight;
}


double totalw(vertex u, int nVtx, vertex* adj, int* xadj, double* weight) {
	double total_weight = 0;
	int* mark = new int[nVtx];
	int* bfsorder = new int[nVtx];
	memset (mark, 0, sizeof(int) * nVtx);
	int endofbfsorder = 1;
	bfsorder[0] = u;
	int cur = 0;
	mark[u] = 1;
	while (cur != endofbfsorder) {
		vertex v = bfsorder[cur];
		total_weight += weight[v];
		for (myindex j = xadj[v]; j < xadj[v + 1]; j++) {
			vertex w = adj[j];
			if ((w != -1) && (mark[w] == 0)) {
				bfsorder[endofbfsorder++] = w;
				mark[w] = 1;
			}
		}
		cur++;
	}

	delete[] mark;
	delete[] bfsorder;
	return total_weight;
}


int totalsize(vertex u, int nVtx, vertex* adj, int* xadj) {
	int total_size = 0;
	int* mark = new int[nVtx];
	int* bfsorder = new int[nVtx];
	memset (mark, 0, sizeof(int) * nVtx);
	int endofbfsorder = 1;
	bfsorder[0] = u;
	int cur = 0;
	mark[u] = 1;
	while (cur != endofbfsorder) {
		vertex v = bfsorder[cur];
		total_size++;
		for (myindex j = xadj[v]; j < xadj[v + 1]; j++) {
			vertex w = adj[j];
			if ((w != -1) && (mark[w] == 0)) {
				bfsorder[endofbfsorder++] = w;
				mark[w] = 1;
			}
		}
		cur++;
		if (cur > 1)
			break;
	}

	delete[] mark;
	delete[] bfsorder;
	return total_size;
}


void graph_check (int* xadj, vertex* adj, int nVtx) {
	for (int i = 0; i < nVtx; i++) {
		for (myindex j =(xadj)[i]; j <(xadj)[i+1]; j++) {
			int v = (adj)[j];
			if (v != -1) {
				int flag = 0;
				for (myindex k = (xadj)[v]; k < (xadj)[v+1]; k++) {
					if ((adj)[k] == i) {
						flag = 1;
						break;
					}
				}
				if (flag == 0) {
					printf("%d and %d do not match each other\n",i+1,v+1);
					break;
				}
			}
		}
	}
}


int base_bc(int nVtx, int *xadj, int *adj, Betweenness *bc,
		int maxvertex, int nTry, //algo parameter
		util::timestamp& totaltime, util::timestamp& phase1time, util::timestamp& phase2time
#ifdef PARTIAL_BC
, double partial_bc_factor
#endif
) {

	int  cnt=0;
	char *p = (char*) &totaltime;
	char *p1 = (char*) &phase1time;
	char *p2 = (char*) &phase2time;

	util::timestamp *pout = (util::timestamp*) p;
	util::timestamp *p1time = (util::timestamp*) p1;
	util::timestamp *p2time = (util::timestamp*) p2;

	for(int Try = 0; Try < THROW_AWAY+nTry; Try++) {
		for (vertex i=0; i<nVtx; i++)
			bc[i] = 0.;

		util::timestamp t1;

		vertex* bfsorder = new vertex[nVtx];
		vertex* Pred = new vertex[xadj[nVtx]];
		int* endpred = new int[nVtx];
		int* level = new int[nVtx];
		pathnumber* sigma = new pathnumber[nVtx];
		Betweenness* delta = new Betweenness[nVtx];

#ifdef PARTIAL_BC
		maxvertex = (int)(partial_bc_factor * (double)(maxvertex));
#endif

		for (vertex source=0; source<maxvertex; source++) {
			util::timestamp t_det_1;
			int endofbfsorder = 1;
			bfsorder[0] = source;

			for (int i=0; i<nVtx; i++)
				endpred[i] = xadj[i];

			for (int i=0; i< nVtx; i++)
				level[i] = -2;
			level[source] = 0;

			for (int i=0; i< nVtx; i++)
				sigma[i] = 0;
			sigma[source] = 1;

			//step 1: build shortest path graph
			int cur = 0;
			while (cur != endofbfsorder) {
				vertex v = bfsorder[cur];
				assert (level[v] >= 0);
				for (myindex j=xadj[v]; j!=xadj[v+1]; j++) {
					vertex w = adj[j];
					if (level[w] < 0) {
						level[w] = level[v]+1;
						bfsorder[endofbfsorder++] = w;
						//                         assert (endofbfsorder <= nVtx);
					}
					if (level[w] == level[v]+1) {
						sigma[w] += sigma[v];
					}
					else if (level[w] == level[v] - 1) {
						Pred[endpred[v]++] = w;
					}
				}
				cur++;
			}

			for (int i=0; i< nVtx; i++) {
				delta[i] = 0.;
			}

			util::timestamp t_det_2;

			//step 2: compute betweenness
			for (int i = endofbfsorder-1; i>0; i--) {
				vertex w = bfsorder[i];
				for (int j=xadj[w]; j != endpred[w]; j++) {
					vertex v = Pred[j];
					delta[v] +=( sigma[v]*(1+delta[w]))/sigma[w];
				}
				bc[w] += delta[w];
#ifdef DEBUG
				printf("source %d adds bc[%d]: %lf\n",source+1,w+1,delta[w]);
#endif
			}
#ifdef DEBUG
			printf("source:%d, delta[4]:%lf\n",source+1,delta[3]);
#endif

			util::timestamp t_det_3;
			*p2time += (t_det_3 - t_det_2);
			*p1time += (t_det_2 - t_det_1);
		}

		util::timestamp t2;

		for (vertex source=0; source<nVtx; source++) {
			printf("bc[%d] : %lf\n",source+1,bc[source]);
		}

		delete[] bfsorder;
		delete[] Pred;
		delete[] level;
		delete[] sigma;
		delete[] delta;

		if (Try >= THROW_AWAY) {
			*pout += ((t2-t1));
		}
	}

	return cnt;
}


