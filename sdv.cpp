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


void remove_covs (int len, vertex* component, int cov_i, vertex* clique_only_v, vertex* reversecomp, Bucket* bs, int nVtx,
		double* weight, int* xadj,
		vertex* adj, vertex* bfsorder, int* endpred, int* level, pathnumber* sigma,
		vertex* Pred, Betweenness* delta, Betweenness* bc, int* numof_removed_edges, int* idv_track,
		int* identical_sets_c, idv_info** identical_sets, int next_idvset_id, double* totalw_of_covs,
		double* total_weights_of_each_comp, int* comp_ids_of_each_v, double* CCs, double* ff) {

	int* tmark = (int *)malloc(sizeof(int) * 2 * nVtx);
	vertex* tbfsorder = (vertex *)malloc(sizeof(vertex) * 2 * nVtx);

	/**
	 * @todo 
	 * 為甚麼要 temp giving back : 因為要先算 clique_vertex 的 BC
	 * 
	 * [Note]
	 * next_idvset	: 總共有多少個 idvset
	*/
	// Temp giving back
	for (int i = 0; i < next_idvset_id; i++) {
		if (identical_sets_c[i] > 0) {
			Betweenness bc_of_idv_repr = bc[identical_sets[i][1].id];
			for (int j = 2; j < identical_sets_c[i]; j++) {
				bc[identical_sets[i][j].id] += bc_of_idv_repr;
#ifdef BCCOMP_DBG
				printf("%lf is added back to bc[%d] from bc[%d], at temp give back\n",
						bc_of_idv_repr, identical_sets[i][j].id+1, identical_sets[i][1].id+1);
#endif
			}
		}
	}


	/**
	 * 
	*/
	for (int i = 0; i < cov_i; i++) {
		int u = clique_only_v[i];
#ifdef DEBUG
		printf("cov: %d\n",u+1);
#endif
		int rev_u = reversecomp[u];
		int degree = bs->values[rev_u]; //current degree of u
		if ((degree == 2) || (degree == 3)) {

			double weight_of_u;
			if (idv_track[u] == -1)
				weight_of_u = weight[u];
			else
				weight_of_u = identical_sets[idv_track[u]][0].weight; //u所在的 identical set 的全部 weight


			//這邊把 u 的 bc先算完
			bc_comp_for_cov (u, nVtx, weight, xadj, adj, bfsorder, endpred, level, sigma, Pred, delta, bc, idv_track, identical_sets_c, identical_sets, CCs, ff);

			/**
			 * 斷開 clique_vertex u 跟他的所有鄰居的連結
			*/
			for (myindex j = xadj[u]; j < xadj[u+1]; j++) {
				int w = adj[j];
				if (w != -1) {
					adj[j] = -1;
					(*numof_removed_edges)++;
					Zoltan_Bucket_DecVal(bs, rev_u);
					int rev_w = reversecomp[w];
					for (myindex k = xadj[w]; k < xadj[w+1]; k++) {
						if (adj[k] == u) {
							adj[k] = -1;
							Zoltan_Bucket_DecVal(bs, rev_w);
							break;
						}
					}
				}
			}

			//把 component 的 weight 扣掉 u的weight
			total_weights_of_each_comp[comp_ids_of_each_v[u]] -= weight_of_u;
		}
	}

	free(tmark);
	free(tbfsorder);


	/**
	 * @todo
	 * 這邊又把剛剛 + 的 BC 扣回去
	*/
	// Taking back
	for (int i = 0; i < next_idvset_id; i++) {
		if (identical_sets_c[i] > 0) {
			Betweenness bc_of_idv_repr = bc[identical_sets[i][1].id];
			for (int j = 2; j < identical_sets_c[i]; j++) {
				bc[identical_sets[i][j].id] -= bc_of_idv_repr;
#ifdef BCCOMP_DBG
				printf("%lf is removed from to bc[%d] (as amount of bc[%d]), at take back\n",
						bc_of_idv_repr, identical_sets[i][j].id+1, identical_sets[i][1].id+1);
#endif
			}
		}
	}

	for (int i = 0; i < len; i++) {
		int u = component[i];
		if (u != -1) {
			int flag = 0;
			for (myindex j = xadj[u]; j < xadj[u+1]; j++) {
				if (adj[j] != -1) {
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				component[i] = -1;
		}
	}
}

/**
 * @brief
 * 輸入一個source，算出該點的BC
 * @param
 * source 要當做起點的 node
*/
void bc_comp_for_cov (int source, int nvtx, double* weight, int* xadj,
		vertex* adj, vertex* bfsorder, int* endpred, int* level, pathnumber* sigma,
		vertex* Pred, Betweenness* delta, Betweenness* bc, int* idv_track,
		int* identical_sets_c, idv_info** identical_sets, double* CCs, double* ff) {

#ifdef BCCOMP_DBG
	printf("%d is cov!!!\n",source+1);
#endif

	if (idv_track[source] == -1) {
		int endofbfsorder = 1;
		bfsorder[0] = source;

		// for (int i = 0; i < nvtx; i++)
		// 	endpred[i] = xadj[i];

		for (int i = 0; i < nvtx; i++)
			level[i] = -2;
		level[source] = 0;

		// for (int i = 0; i < nvtx; i++)
		// 	sigma[i] = 0;
		// sigma[source] = 1;

		// step 1: build shortest path graph
		int cur = 0;
		while (cur != endofbfsorder) {
			vertex v = bfsorder[cur];
			assert (level[v] >= 0);
			int idx_of_v = idv_track[v];
			for (myindex j = xadj[v]; j < xadj[v + 1]; j++) {
				vertex w = adj[j];
				if (w != -1)
				{
					if (level[w] < 0) {
						level[w] = level[v]+1;
						bfsorder[endofbfsorder++] = w;
						
						#pragma region CC

						if(idv_track[w] == -1){
							CCs[source] += ff[w] + level[w] * weight[w];
						}
						else{
							int idx_of_w = idv_track[w];
							CCs[source] += identical_sets[idx_of_w][0].idv_ff + level[w] * identical_sets[idx_of_w][0].weight;
						}

						CCs[w] += ff[source] + level[w] * weight[source];

						#pragma endregion //CC
					}

					// if (level[w] == level[v]+1) {
					// 	if (idv_track[v] == -1){
					// 		CCs[w] += ff[v] + level[w] * weight[v];
					// 		// sigma[w] += sigma[v];
					// 	}
					// 	else {
					// 		CCs[w] += identical_sets[idx_of_v][0].idv_ff + level[w] * identical_sets[idx_of_v][0].weight;
					// 		// sigma[w] += sigma[v] * (identical_sets_c[idx_of_v] - 1);
					// 	}
					// }
					// else if (level[w] == level[v] - 1) {
					// 	Pred[endpred[v]++] = w;
					// }
				}
			}
			cur++;
		}

		// for (int i = 0; i < nvtx; i++) {
		// 	delta[i] = 0.;
		// }

		// for (int i = 0; i < nvtx; i++) {
		// 	delta[i] += (weight[i]) - 1;
		// }

		// step 2: compute betweenness
// 		for (int i = endofbfsorder - 1; i > 0; i--) {
// 			vertex w = bfsorder[i];
// 			if (idv_track[w] == -1) {
// 				for (myindex j = xadj[w]; j < endpred[w]; j++) {
// 					vertex v = Pred[j];
// 					delta[v] += sigma[v] * (1 + delta[w]) / sigma[w];
// 				}
// 			}
// 			else { // phase-2 hack for handling idv vertices
// 				//                 printf("%d is idv\n",w+1);
// 				int idx_of_w = idv_track[w];
// 				for (myindex j = xadj[w]; j < endpred[w]; j++) {
// 					vertex v = Pred[j];
// 					delta[v] += sigma[v] * (((identical_sets_c[idx_of_w] - 1) *
// 							(delta[w] + 1 - identical_sets[idx_of_w][1].weight)) + identical_sets[idx_of_w][0].weight) / sigma[w];
// 				}
// 			}

// 			bc[w] += weight[source] * delta[w];
// 			bc[w] += weight[source] * (delta[w] - (weight[w] - 1));

// #ifdef BCCOMP_DBG
// 			printf("weight[%d]: %lf, delta[%d]: %lf\n",source+1,weight[source],source+1,delta[source]);
// 			printf("weight[%d]: %lf, delta[%d]: %lf\n",w+1,weight[w],w+1,delta[w]);
// 			printf("source %d adds bc[%d]: %lf\n",source+1,w+1,((weight[source])) * (delta[w]));
// 			printf("source %d adds bc[%d]: %lf\n",source+1,w+1,((weight[source])) * (delta[w] - ((weight[w]) - 1)));
// #endif


// 			if (idv_track[w] != -1) {
// 				int idv_of_w = idv_track[w];
// 				for (int i = 2; i < identical_sets_c[idv_of_w]; i++) {
// 					bc[identical_sets[idv_of_w][i].id] += weight[source] * (delta[w] - weight[w] + identical_sets[idv_of_w][i].weight) ;
// 					bc[identical_sets[idv_of_w][i].id] += weight[source] * (delta[w] - (weight[w] - 1));
// #ifdef BCCOMP_DBG
// 					printf("%lf is added to bc[%d]\n", weight[source] * (delta[w] - weight[w] + identical_sets[idv_of_w][i].weight),
// 							identical_sets[idv_of_w][i].id+1);
// 					printf("%lf is added to bc[%d]\n", weight[source] * (delta[w] - (weight[w] - 1)), identical_sets[idv_of_w][i].id+1);
// #endif
// 				}
// 			}
// 		}
// 		bc[source] += ((weight[source] - 1)) * (delta[source] - ((weight[source]) - 1));

#ifdef BCCOMP_DBG
		printf("source %d adds bc[%d]: %lf\n",source+1,source+1,((weight[source]) - 1) * (delta[source] - ((weight[source]) - 1)));
#endif

	}
	else {

		int idx_of_source = idv_track[source];

		int endofbfsorder = 1;
		bfsorder[0] = source;

		// for (int i = 0; i < nvtx; i++)
		// 	endpred[i] = xadj[i];

		for (int i = 0; i < nvtx; i++)
			level[i] = -2;
		level[source] = 0;

		// for (int i = 0; i < nvtx; i++)
		// 	sigma[i] = 0;
		// sigma[source] = 1;

		// step 1: build shortest path graph
		int cur = 0;
		while (cur != endofbfsorder) {
			vertex v = bfsorder[cur];
			assert (level[v] >= 0);
			// int idx_of_v = idv_track[v];
			for (myindex j = xadj[v]; j < xadj[v + 1]; j++) {
				vertex w = adj[j];
				if (w != -1)
				{
					if (level[w] < 0) {
						level[w] = level[v]+1;
						bfsorder[endofbfsorder++] = w;

						#pragma region CC

						if(idv_track[w] == -1){
							CCs[source] += ff[w] + level[w] * weight[w];
						}
						else{
							int idx_of_w = idv_track[w];
							CCs[source] += identical_sets[idx_of_w][0].idv_ff + level[w] * identical_sets[idx_of_w][0].weight;
						}
						CCs[w] += identical_sets[idx_of_source][0].idv_ff + level[w] * identical_sets[idx_of_source][0].weight;

						#pragma endregion //CC

					}
					// if (level[w] == level[v]+1) {
					// 	if ((idv_track[v] != -1) && (v != source)) {
					// 		sigma[w] += sigma[v] * (identical_sets_c[idx_of_v] - 1);
					// 	}
					// 	else
					// 		sigma[w] += sigma[v];
					// }
					// else if (level[w] == level[v] - 1) {
					// 	Pred[endpred[v]++] = w;
					// }
				}
			}
			cur++;
		}

		// for (int i = 0; i < nvtx; i++) {
		// 	delta[i] = 0.;
		// }

		// for (int i = 0; i < nvtx; i++) {
		// 	delta[i] += (weight[i]) - 1;
		// }

		// step 2: compute betweenness
// 		for (int i = endofbfsorder - 1; i > 0; i--) {
// 			vertex w = bfsorder[i];
// 			if ((idv_track[w] != -1) && (w != source)) {
// 				int idx_of_w = idv_track[w];
// 				for (myindex j = xadj[w]; j < endpred[w]; j++) {
// 					vertex v = Pred[j];
// 					delta[v] += sigma[v] * (((identical_sets_c[idx_of_w] - 1) *
// 							(delta[w] + 1 - identical_sets[idx_of_w][1].weight)) + identical_sets[idx_of_w][0].weight) / sigma[w];
// 				}
// 			}
// 			else { // phase-2 hack for handling idv vertices
// 				for (myindex j = xadj[w]; j < endpred[w]; j++) {
// 					vertex v = Pred[j];
// 					delta[v] += sigma[v] * (1+delta[w])/sigma[w];
// 				}
// 			}

// 			bc[w] += identical_sets[idx_of_source][0].weight * delta[w];
// 			bc[w] += identical_sets[idx_of_source][0].weight * (delta[w] - (weight[w] - 1));
// #ifdef BCCOMP_DBG
// 			printf("source %d adds bc[%d]: %lf\n",source+1, w+1, identical_sets[idx_of_source][0].weight * delta[w]);
// 			printf("source %d adds bc[%d]: %lf\n",source+1, w+1, identical_sets[idx_of_source][0].weight * (delta[w] - (weight[w] - 1)));
// #endif

// 			if (idv_track[w] != -1) {
// 				int idv_of_w = idv_track[w];
// 				for (int i = 2; i < identical_sets_c[idv_of_w]; i++) {
// 					bc[identical_sets[idv_of_w][i].id] += identical_sets[idx_of_source][0].weight * (delta[w] - weight[w] + identical_sets[idv_of_w][i].weight);
// 					bc[identical_sets[idv_of_w][i].id] += identical_sets[idx_of_source][0].weight * (delta[w] - (weight[w] - 1));
// #ifdef BCCOMP_DBG
// 					printf("%lf is added to bc[%d]\n", identical_sets[idx_of_source][0].weight * (delta[w] - weight[w] + identical_sets[idv_of_w][i].weight), identical_sets[idv_of_w][i].id+1);
// 					printf("%lf is added to bc[%d]\n", identical_sets[idx_of_source][0].weight * (delta[w] -
// 							weight[w] - 1), identical_sets[idv_of_w][i].id+1);
// #endif
// 				}
// 			}
// 		}

// 		bc[source] += ((weight[source] - 1)) * (delta[source] - ((weight[source]) - 1));
// #ifdef BCCOMP_DBG
// 		printf("%lf is added to bc[%d]\n", ((weight[source] - 1)) * (delta[source] - ((weight[source]) - 1)), source+1);
// #endif

// 		for (int i = 2; i < identical_sets_c[idx_of_source]; i++) {
// 			bc[identical_sets[idx_of_source][i].id] += ((identical_sets[idx_of_source][i].weight - 1)) * (delta[source] - (weight[source] - 1));
// #ifdef BCCOMP_DBG
// 			printf("%lf , come ion, is added to bc[%d]\n", ((identical_sets[idx_of_source][i].weight - 1)) *
// 					(delta[source] - (weight[source] - 1)), identical_sets[idx_of_source][i].id+1);
// #endif
// 		}

		idv_track[source] = -1;
		identical_sets_c[idx_of_source] = 0;
		identical_sets[idx_of_source][0].id = -1;

	}
}

/**
 * @brief
 * 用三個for迴圈 判斷一個點是否是clique的一部分
*/
int is_part_of_clique (int u, int* xadj, vertex* adj, int clique_neig_num) {

	for (myindex i = xadj[u]; i < xadj[u+1]; i++) {
		int v = adj[i];
		if (v != -1) {
			if (xadj[v+1] - xadj[v] < clique_neig_num)
				return 0;
			for (myindex k = xadj[u]; k < xadj[u+1]; k++) {
				if ((adj[k] != v) && (adj[k] != -1)) {
					int flag = 0;
					for (myindex j = xadj[v]; j < xadj[v+1]; j++) {
						if (adj[j] == adj[k]) {
							flag = 1;
							break;
						}
					}
					if (flag == 0)
						return 0;
				}
			}
		}
	}
	return 1;
}

/**
 * @brief
 * 找出這個component中的所有nodes並存入clique_only_v中，
 * 這些nodes 一定可以組成某些clique
 * 
 * cov_i : clique_only_vertex 的個數
*/
void find_clique_only_vertices (vertex* clique_only_v, int* cov_i, vertex* component, int len, int* xadj, vertex* adj) {
	// 4 and 3-clique vertices are detected
	for (int i = 0; i < len; i++) {
		int u = component[i];
		int t = 0;
		if (u != -1) {
			//計數 u 有幾個 neighbor
			for (myindex j = xadj[u]; j < xadj[u+1]; j++) {
				if (adj[j] != -1)
					t++;
			}
			if ((t == 2) || (t == 3)) {
				if (is_part_of_clique (u, xadj, adj, t) == 1) {
					clique_only_v[(*cov_i)++] = u;
				}
			}
		}
	}
	return;
}


