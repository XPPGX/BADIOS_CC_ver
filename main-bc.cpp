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

/**
 * @brief 
 * [Func OK]		:	D1, bridge
 * [Func done yet]	:	AP, idv, side, order
*/

#include "bc-seq-brandes.h"


int main_bc(int nVtx, int **pxadj, int **padj, int **ptadj, Betweenness **bc,
		int nTry,
		util::timestamp& totaltime, util::timestamp& preproctime, util::timestamp& phase1time, util::timestamp& phase2time,
		util::timestamp& deg1remtime, util::timestamp& bridgedettime, util::timestamp& bridgeremtime,
		util::timestamp& cliquedettime,
		util::timestamp& cliqueremtime,
		util::timestamp& artpdettime, util::timestamp& artpremtime,
		util::timestamp& idvdettime, util::timestamp& idvremtime, util::timestamp& bfsordertime,
		int* numof_removed_edges, int* numof_art_points, int* numof_newly_created_vertices,
		int* numof_identical_vertices,
		int* biggest_cc_before, int* biggest_cc_after, int* num_comp, bool dc, bool da, bool db, bool dd, bool di
#ifdef PARTIAL_BC
, double partial_bc_factor
#endif
) {

#ifdef CC_PROFILE
	int numof_CCs_atb= 0;
	int maxCC_nvtx_atb = 0;
	int maxCC_nedge_atb = 0;
	int total_nvtx_atb = 0;
	int total_nedge_atb = 0;
	printf("---------CC INFORMATION AT THE BEGINNING---------\n");
	extract_cc_info (nVtx, (*pxadj), (*padj), &numof_CCs_atb, &maxCC_nvtx_atb,
					&maxCC_nedge_atb, &total_nvtx_atb, &total_nedge_atb);
	printf("-------------------------------------------------\n\n\n\n\n");
#endif

	int  cnt=0;
	int updatedlen;
	char *p = (char*) &totaltime;
	char *p1 = (char*) &preproctime;
	char *p2 = (char*) &cliquedettime;
	char *p3 = (char*) &cliqueremtime;

	char *p4 = (char*) &deg1remtime;

	char *p5 = (char*) &bridgedettime;
	char *p6 = (char*) &bridgeremtime;

	char *p7 = (char*) &artpdettime;
	char *p8 = (char*) &artpremtime;

	char *p9 = (char*) &bfsordertime;

	char *p10 = (char*) &idvdettime;
	char *p11 = (char*) &idvremtime;

	util::timestamp *pout = (util::timestamp*) p;
	util::timestamp *preproc = (util::timestamp*) p1;
	util::timestamp *clique_det = (util::timestamp*) p2;
	util::timestamp *clique_rem = (util::timestamp*) p3;

	util::timestamp *deg1_rem = (util::timestamp*) p4;

	util::timestamp *bridge_det = (util::timestamp*) p5;
	util::timestamp *bridge_rem = (util::timestamp*) p6;

	util::timestamp *artp_det = (util::timestamp*) p7;
	util::timestamp *artp_rem = (util::timestamp*) p8;

	util::timestamp *bfsorder_t = (util::timestamp*) p9;

	util::timestamp *idv_det = (util::timestamp*) p10;
	util::timestamp *idv_rem = (util::timestamp*) p11;


	for(int Try = 0; Try < THROW_AWAY+nTry; Try++) {
		for (vertex i=0; i<nVtx; i++)
			(*bc)[i] = 0.;

		util::timestamp t1;

		int num_edges =(*pxadj)[nVtx];

		int initnVtx = nVtx;
		int initnedges =(*pxadj)[nVtx];

		// may myre_allocated if art_p is enabled
		int sizeofnewxadj = 2 * (nVtx + 1);
		int sizeofnewadjorpred = 2 * (*pxadj)[nVtx];
		vertex* newxadj = (vertex *)malloc(sizeof(vertex) * 2 * (nVtx + 1));
		vertex* newadj = (vertex *)malloc(sizeof(vertex) * 2 * (*pxadj)[nVtx]);
		memset(newxadj, 0, sizeof(int) * (2*(nVtx+1)));
		memset(newadj, 0, sizeof(int) * (2*(*pxadj)[nVtx]));
		vertex* bfsorder = (vertex *)malloc(sizeof(vertex) * 2 * nVtx);
		vertex* Pred = (vertex *)malloc(sizeof(vertex) * 2 *(*pxadj)[nVtx]);
		int* endpred = (int *)malloc(sizeof(int) * 2 * nVtx);
		int* level = (int *)malloc(sizeof(int) * 2 * nVtx);
		pathnumber* sigma = (pathnumber *)malloc(sizeof(pathnumber) * 2 * nVtx);
		Betweenness* delta = (Betweenness *)malloc(sizeof(Betweenness) * 2 * nVtx);
		int* mark = (int *)malloc(sizeof(int) * 2 * nVtx);

		#pragma region CC_var
		double* CCs = (double*)calloc(sizeof(double), 2 * nVtx);
		double* ff = (double*)calloc(sizeof(double), 2 * nVtx);
		//weight 就是 reach 
		#pragma endregion //CC_var

		int* tmark = (int *)malloc(sizeof(int) * 2 * nVtx);
		vertex* tbfsorder = (vertex *)malloc(sizeof(vertex) * 2 * nVtx);

		memset (mark, 0, sizeof(int) * 2 * nVtx);
		double* weight = (double *)malloc(sizeof(double) * 2 * nVtx);

		// degree-0 vertex elimination
		int zerodegrees = 0;
		for (vertex i = 0; i < nVtx; i++) {
			if ((*pxadj)[i+1] ==(*pxadj)[i]) {
				(*bc)[i] = 0.;
				mark[i] = 1;
				zerodegrees++;
			}
		}

		vertex* ordered_comp = (vertex *) malloc (sizeof(vertex) * 2 * nVtx);
		memset (ordered_comp, -1, sizeof(vertex) * 2 * nVtx);
		vertex* component = (vertex *) malloc (sizeof(vertex) * 2 * nVtx);
		vertex* reversecomp = (vertex *) malloc (sizeof(vertex) * 2 * nVtx);
		vertex* reverse_ordered_comp = (vertex *) malloc (sizeof(vertex) * 2 * nVtx);
		memset (reverse_ordered_comp, -1, sizeof(vertex) * 2 * nVtx);
		memset (reversecomp, -1, sizeof(vertex) * 2 * nVtx);
		double* ordered_weight = (double *) malloc (sizeof(double) * 2 * nVtx);

		#pragma region CC
		double* ordered_ff = (double*)malloc(sizeof(double) * 2 * nVtx);
		#pragma endregion //CC

		card_info* ordered_cardinality = (card_info *)malloc(sizeof(card_info) * 2 * nVtx);
		int* all_graphs_xadj = (int *) malloc (sizeof(int) * 2 * nVtx);
		int* all_graphs_len = (int *) malloc (sizeof(int) * 2 * nVtx);
		int* labels = (int *) malloc (sizeof(int) * 2 * nVtx);
		vertex* bridges = (vertex *) malloc (sizeof(vertex) * 2 *(*pxadj)[nVtx]);
		int* nd = (int *) malloc (sizeof(int) * 2 * nVtx); //紀錄BFS traverse tree的node下面有幾個node
		int* l = (int *) malloc (sizeof(int) * 2 * nVtx); 
		int* h = (int *) malloc (sizeof(int) * 2 * nVtx);
		vertex* art_points = (vertex *) malloc (sizeof(vertex) * 2 * nVtx);
		vertex* art_track = (vertex *) malloc (sizeof(vertex) * 2 * nVtx);
		memset (art_track, -1, sizeof(vertex) * 2 * nVtx);

		int idv_sets_size = INIT_SIZE_FOR_IDV_SETS; // initial estimate for total number of idv sets
		int one_set_size = INIT_SIZE_FOR_IDV_SETS; // initial estimate for number of (vertices, weights) in an idv set
		int idv_track_size = 2 * nVtx;
		int* idv_track;

		// -1 if a vertex is not in an idv set, otherwise id for idv_set
		idv_track = (int *) malloc (sizeof(int) * idv_track_size);
		for (int i = 0; i < idv_track_size; i++) {
			idv_track[i] = -1;
		}

		int* identical_sets_sz;
		identical_sets_sz = (int*) malloc(idv_sets_size * sizeof(int));
		for (int i = 0; i < idv_sets_size; i++) {
			identical_sets_sz[i] = INIT_SIZE_FOR_IDV_SETS;
		}

		int* identical_sets_c;
		identical_sets_c = (int*) calloc(idv_sets_size, sizeof(int));
		idv_info** identical_sets;
		identical_sets = (idv_info**) malloc(sizeof(idv_info*) * idv_sets_size);
		for (int i = 0; i < idv_sets_size; i++)
			identical_sets[i] = (idv_info*) malloc(sizeof(idv_info) * one_set_size);
		int next_idvset_id = 0;

		int artc = 0;
		int bridges_c;
		all_graphs_xadj[0] = 0;
		int agx = 1;
		int k = 0, c = 0;
		int numofremovedbridges = 0;
		int size1 = 2 * initnVtx;
		int size2 = 2 * (initnVtx + 1);
		int size3 = 2 * initnedges;
		int size4 = 2 * initnVtx;

		int startingnVtx = nVtx;

#ifdef REDUCTION_STEPS
		int comp_id = 0;
#endif

#ifdef JUST_PREPROC
		int counter_to_count_components_before_reduction = 1;
#endif


		double* total_weights_of_each_comp = (double *) malloc(sizeof(double) * nVtx);
		int* comp_ids_of_each_v = (int *) calloc(2 * nVtx, sizeof(int));
		int comp_no = 1;

		for (vertex i = 0; i < startingnVtx; i++) {
			if(mark[i] == 0) {
				memset(component, -1, sizeof(vertex) * 2 * startingnVtx);
				// detecting connected component
				int compcounter = 0;
				component[compcounter++] = i;
				mark[i] = 1;
				int cur = compcounter - 1;
				double total_weight = 0;
				while (cur != compcounter) {
					vertex v = component[cur];
					weight[v] = 1;
					comp_ids_of_each_v[v] = comp_no;
					total_weight += weight[v];
					for (myindex j=(*pxadj)[v]; j < (*pxadj)[v+1]; j++) {
						vertex w =(*padj)[j];
						if (mark[w] == 0) {
							mark[w] = 1;
							component[compcounter] = w;
							compcounter++;
						}
					}
					cur++;
				}
				total_weights_of_each_comp[comp_no] = total_weight;
				comp_no++;
			}
		}

		memset (mark, 0, sizeof(int) * 2 * nVtx);

		for (vertex i = 0; i < startingnVtx; i++) {
			if(mark[i] == 0) {
#ifdef REDUCTION_STEPS
				printf("\nFOR COMPONENT-ID: %d\n", ++comp_id);
#endif
				memset(component, -1, sizeof(vertex) * 2 * startingnVtx);
				// detecting connected component
				int compcounter = 0;
				component[compcounter++] = i;
				mark[i] = 1;
				int cur = compcounter - 1;
				while (cur != compcounter) {
					vertex v = component[cur];
					for (myindex j=(*pxadj)[v]; j < (*pxadj)[v+1]; j++) {
						vertex w =(*padj)[j];
						if (mark[w] == 0) {
							mark[w] = 1;
							component[compcounter] = w;
							compcounter++;
						}
					}
					cur++;
				}

				if (*biggest_cc_before < compcounter)
					*biggest_cc_before = compcounter;

#ifdef JUST_PREPROC
				printf("Comp-%4d: %10d\n", counter_to_count_components_before_reduction++, compcounter);
#endif


				int len = compcounter;
				Bucket bs;
				int max_degree = 0;
				for (int i = 0; i < len; i++) {
					int j = component[i];
					if (((*pxadj)[j+1] -(*pxadj)[j]) > max_degree)
						max_degree =(*pxadj)[j+1] -(*pxadj)[j];
				}
				bs = Zoltan_Bucket_Initialize(max_degree+1, len);

				for (vertex i = 0; i < len; i++) {
					int j = component[i];
					int degree_of_i =(*pxadj)[j+1] -(*pxadj)[j];
					reversecomp[component[i]] = i;
					Zoltan_Bucket_Insert(&bs, i, degree_of_i);
				}

				int isReduced = 1;
				int nextvid = nVtx;
				int round = 0;

				//[Preprocess]
		
				while (isReduced) {
#ifdef REDUCTION_STEPS
					printf("**********ROUND: %d for COMPONENT: %d**********\n",round++,comp_id);
#endif

#ifdef DEBUG
					//to be removed
					graph_check((*pxadj), (*padj), nVtx);
#endif

					int numofremovededgesatbeginning = *numof_removed_edges;
					int numofartpointsatbeginning = *numof_art_points;

					//degree-1 removal
					{
#ifdef REDUCTION_STEPS
						int a = *numof_removed_edges;
#endif	
						printf("[D1]\n");
						util::timestamp ta;
						if (!dd) {
							remove_degree_1s (nVtx, component, reversecomp, &bs, weight, (*pxadj), (*padj),
									(*ptadj), idv_track, identical_sets_c, identical_sets, (*bc),
									numof_removed_edges, tmark, tbfsorder, total_weights_of_each_comp, comp_ids_of_each_v, CCs, ff);
						}
						util::timestamp tb;
						if (Try >= THROW_AWAY) {
							*deg1_rem += (tb - ta);
						}

#ifdef REDUCTION_STEPS
						int b = *numof_removed_edges;
						printf("at DEG-1 REMOVAL, %d num of edges are removed\n", b-a);
#endif
						printf("\n");
					}


					//bridge removal (also includes newly created degree-1 removal)
					{
#ifdef REDUCTION_STEPS
						int a = *numof_removed_edges;
#endif
						updatedlen = 0;
						for (int i = 0; i < len; i++) {
							if (component[i] != -1)
								updatedlen++;
						}
						numofremovedbridges -= *numof_removed_edges;

#ifdef DEBUG
						graph_check((*pxadj), (*padj), nVtx);
#endif
						printf("[Bridge]\n");
						if (!db && (updatedlen > 1)) {
							util::timestamp tc;
							bridges_c = 0;

							bridge_detection_multiple_cc (len, nVtx, labels, nd, l, h, component, (*pxadj), (*padj),
														bridges, &bridges_c);
							util::timestamp td;
							if (Try >= THROW_AWAY) {
								*bridge_det += (td - tc);
							}
							
							bridge_removal (nVtx, len, bridges, bridges_c, numof_removed_edges, &bs, component,
									reversecomp, (*padj), (*pxadj), (*ptadj), weight, (*bc), idv_track, identical_sets_c,
									identical_sets, tmark, tbfsorder, dd, total_weights_of_each_comp, comp_ids_of_each_v,
									&comp_no, CCs, ff);
							util::timestamp te;
							if (Try >= THROW_AWAY) {
								*bridge_rem += (te - td);
							}
						}
						numofremovedbridges += *numof_removed_edges;
#ifdef REDUCTION_STEPS
						int b = *numof_removed_edges;
						printf("at BRIDGE REMOVAL, %d num of edges are removed\n", b-a);
#endif
						printf("\n");
						// exit(1);
					}

					/**
					 * 用下面這個註解掉的 for迴圈，檢查 ff 中是否有 nan
					 * @brief
					 * 在還沒有identical vertex的情況下，目前沒有 ff 是 nan
					 * 
					*/
					// for(int i = 0 ; i < nVtx ; i ++){
					// 	if(std::isnan(ff[i])){
					// 		printf("[ERROR] ff[%d] = %f\n", i, ff[i]);
					// 		exit(1);
					// 	}
					// }


					//articulation point separation
					{
#ifdef REDUCTION_STEPS
						int a = *numof_art_points;
						int c = *numof_newly_created_vertices;
						int e = *numof_removed_edges;
#endif		
						printf("[AP]\n");
						updatedlen = 0;
						for (int i = 0; i < len; i++) {
							if (component[i] != -1)
								updatedlen++;
						}

#ifdef DEBUG
						graph_check((*pxadj), (*padj), nVtx);
#endif
						if (!da && (updatedlen > 1)) {
							util::timestamp ti;
							artc = 0;
							vertex* stack = new vertex[updatedlen];
							int* dfn = new int[nVtx];
							int* l = new int[nVtx];
							int* parent = new int[nVtx];
							int* already_art = new int[nVtx];
							memset (already_art, 0, sizeof(int) * nVtx);
							vertex* mark__comp = new vertex[nVtx];

							articulation_point_detection_multiple_cc (updatedlen, nVtx, art_points, &artc, stack, dfn, l,
									parent, already_art, mark__comp, component, (*padj), (*pxadj));
							
							printf("[AP : detection][Done]\n");
							*numof_art_points += artc;

							delete[] stack;
							delete[] l;
							delete[] dfn;
							delete[] parent;
							delete[] mark__comp;
							delete[] already_art;

							util::timestamp tj;
							if (Try >= THROW_AWAY) {
								*artp_det += (tj - ti);
							}
#ifdef DEBUG
							printf("%d num of artps\n",artc);
#endif
							//articulation point copying
							if (artc > 0) {
								int begin = nVtx;

								articulation_point_copy (&nVtx, artc, &nextvid, &len, &size1, &size2, &size4,
										art_points, art_track, pxadj, padj,
										bc, &num_edges, initnVtx, bfsorder, endpred,
										level, sigma, delta,
										mark, weight, component, ordered_comp, reversecomp,
										reverse_ordered_comp,
										ordered_weight, ordered_cardinality, all_graphs_xadj, all_graphs_len, labels, nd, l,
										h, &bs, &next_idvset_id, &idv_track_size, &idv_track, &identical_sets_c,
										&identical_sets_sz, &identical_sets, &idv_sets_size, one_set_size, tmark, tbfsorder,
										total_weights_of_each_comp, comp_ids_of_each_v, &comp_no, &CCs, &ff);

								int end = nVtx;
								*numof_newly_created_vertices += end - begin;
								if (num_edges > size3) {
									size3 = num_edges;
									bridges = (vertex *) myre_alloc (bridges, (size3) * sizeof(vertex));
								}

#ifdef VLIST
								printf("imm after artp func\n");
								// Graph Check : to be removed
								for (int i = 0; i < nVtx; i++) {
									printf("vertex %d : ", i+1);
									for (myindex j =(*pxadj)[i]; j <(*pxadj)[i+1]; j++) {
										int v = (*padj)[j];
										if (v != -1) {
											printf("%d ", v+1);
											int flag = 0;
											for (myindex k = (*pxadj)[v]; k < (*pxadj)[v+1]; k++) {
												if ((*padj)[k] == i) {
													flag = 1;
													break;
												}
											}
											if (flag == 0) {
												printf("\n%d and %d do not match each other, "
														"remember that nVtx was %d\n",i+1,v+1,startingnVtx);
												break;
											}
										}
									}
									printf("\n");
								}
#endif
								//degree-1 removal
								{
									printf("[AP : D1]\n");
									util::timestamp ta;
									if (!dd) {
										remove_degree_1s (nVtx, component, reversecomp, &bs, weight, (*pxadj), (*padj),
												(*ptadj), idv_track, identical_sets_c, identical_sets, (*bc),
												numof_removed_edges, tmark, tbfsorder,
												total_weights_of_each_comp, comp_ids_of_each_v, CCs, ff);
									}
									util::timestamp tb;
									if (Try >= THROW_AWAY) {
										*deg1_rem += (tb - ta);
										*artp_rem -= (tb - ta);
									}
								}
							}

							util::timestamp tk;
							if (Try >= THROW_AWAY) {
								*artp_rem += (tk - tj);
							}
							// exit(1);
						}
#ifdef REDUCTION_STEPS
						int b = *numof_art_points;
						int d = *numof_newly_created_vertices;
						int f = *numof_removed_edges;
						printf("at ARTP COPY, %d num of art points are found\n", b-a);
						printf("at ARTP COPY, %d num of vertices are created\n", d-c);
						printf("at ARTP COPY, %d num of edges are removed (as a result of deg1)\n", f-e);
#endif
					}

#ifdef VLIST
					// Graph Check : to be removed
					for (int i = 0; i < nVtx; i++) {
						printf("vertex %d : ", i+1);
						for (myindex j =(*pxadj)[i]; j <(*pxadj)[i+1]; j++) {
							int v = (*padj)[j];
							if (v != -1) {
								printf("%d ", v+1);
								int flag = 0;
								for (myindex k = (*pxadj)[v]; k < (*pxadj)[v+1]; k++) {
									if ((*padj)[k] == i) {
										flag = 1;
										break;
									}
								}
								if (flag == 0) {
									printf("%d and %d do not match each other, remember that nVtx was %d\n",
											i+1,v+1,startingnVtx);
									break;
								}
							}
						}
						printf("\n");
					}
#endif

					// identical vertices
					{
#ifdef REDUCTION_STEPS
						int a = *numof_identical_vertices;
						int c = *numof_removed_edges;
#endif
						updatedlen = 0;
						for (int i = 0; i < len; i++) {
							if (component[i] != -1)
								updatedlen++;
						}

						//detection and removal
						if (!di && (updatedlen > 1)) {
#ifdef BCCOMP_DBG
							printf("numof_removed_edges:%d\n",*numof_removed_edges);
#endif
							util::timestamp idv_rem1(0,0);
							util::timestamp idv_rem2(0,0);
							util::timestamp idv_det1(0,0);
							util::timestamp idv_det2(0,0);
							int begin = *numof_removed_edges;
							
							// TYPE-2 detection and removal, sum of neigs and itself is hashed to "hash_size" sized array
							printf("[Idv_detection_and_merge : type 2]\n");
							idv_detection_and_merge(2, len, component, reversecomp, (*pxadj),
									(*padj), nVtx, &idv_sets_size, &one_set_size, &idv_track_size, &next_idvset_id,
									&idv_track, &identical_sets_c, &identical_sets_sz,
									&identical_sets, weight, (*bc), &bs, numof_removed_edges, numof_identical_vertices,
									idv_det2, idv_rem2, CCs, ff);



#ifdef BCCOMP_DBG
							printf("type-2 removals\n");

							for (int i = 0; i < next_idvset_id; i++) {
								if (identical_sets_c[i] > 0) {
									printf("%d-th idv set: \n",i);
									for (int j = 0; j < identical_sets_c[i]; j++) {
										printf("identical_sets[%d][%d] : %d, %lf\n", i, j, identical_sets[i][j].id+1,
												identical_sets[i][j].weight);
									}
								}
							}
#endif

#ifdef VLIST
							for (int i = 0; i < nVtx; i++) {
								printf("vertex %d : ", i+1);
								for (myindex j =(*pxadj)[i]; j <(*pxadj)[i+1]; j++) {
									int v = (*padj)[j];
									if (v != -1) {
										printf("%d ", v+1);
										int flag = 0;
										for (myindex k = (*pxadj)[v]; k < (*pxadj)[v+1]; k++) {
											if ((*padj)[k] == i) {
												flag = 1;
												break;
											}
										}
										if (flag == 0) {
											printf("%d and %d do not match each other, remember that nVtx was %d\n",
													i+1, v+1, startingnVtx);
											break;
										}
									}
								}
								printf("\n");
							}
#endif

							// TYPE-1 detection and removal, sum of neigs is hashed to hash_size sized array
							printf("[Idv_detection_and_merge : type 1]\n");
							idv_detection_and_merge(1, len, component, reversecomp, (*pxadj),
									(*padj), nVtx, &idv_sets_size, &one_set_size, &idv_track_size, &next_idvset_id,
									&idv_track, &identical_sets_c, &identical_sets_sz,
									&identical_sets, weight, (*bc), &bs, numof_removed_edges, numof_identical_vertices,
									idv_det1, idv_rem1, CCs, ff);
							int end = *numof_removed_edges;

#ifdef BCCOMP_DBG
							printf("type-1 removals\n");
							for (int i = 0; i < next_idvset_id; i++) {
								if (identical_sets_c[i] > 0) {
									printf("%d-th idv set: \n",i);
									for (int j = 0; j < identical_sets_c[i]; j++) {
										printf("identical_sets[%d][%d] : %d, %lf\n", i, j, identical_sets[i][j].id+1,
												identical_sets[i][j].weight);
									}
								}
							}
#endif

#ifdef VLIST
							for (int i = 0; i < nVtx; i++) {
								printf("vertex %d : ", i+1);
								for (myindex j =(*pxadj)[i]; j <(*pxadj)[i+1]; j++) {
									int v = (*padj)[j];
									if (v != -1) {
										printf("%d ", v+1);
										int flag = 0;
										for (myindex k = (*pxadj)[v]; k < (*pxadj)[v+1]; k++) {
											if ((*padj)[k] == i) {
												flag = 1;
												break;
											}
										}
										if (flag == 0) {
											printf("%d and %d do not match each other, remember that nVtx was %d\n",
													i+1,v+1,startingnVtx);
											break;
										}
									}
								}
								printf("\n");
							}
#endif
							//degree-1 removal
							if (end > begin) {
								printf("\t[D1][idv D1]\n");
								util::timestamp ta;
								if (!dd) {
									
									remove_degree_1s (nVtx, component, reversecomp, &bs, weight, (*pxadj), (*padj),
											(*ptadj), idv_track, identical_sets_c, identical_sets, (*bc),
											numof_removed_edges, tmark, tbfsorder, total_weights_of_each_comp,
											comp_ids_of_each_v, CCs, ff);
								}
								util::timestamp tb;
								if (Try >= THROW_AWAY) {
									*deg1_rem += (tb - ta);
									*idv_rem += (tb-ta);
								}
								printf("\t[D1][idv D1 done]\n");
							}

							util::timestamp tk;
							if (Try >= THROW_AWAY) {
								*idv_det += idv_det1 + idv_det2;
								*idv_rem += idv_rem1 + idv_rem2;
							}
						}

						// exit(1);
#ifdef REDUCTION_STEPS
						int b = *numof_identical_vertices;
						int d = *numof_removed_edges;
						printf("at IDV REMOVAL, %d num of identical vertices are removed\n", b-a);
						printf("at IDV REMOVAL, %d num of edges are removed (deg1 rem included)\n", d-c);
#endif
					}

					
					int numofremovededges_uptonow = *numof_removed_edges;
					int numofartpoints_uptonow = *numof_art_points;

					isReduced = (numofremovededges_uptonow - numofremovededgesatbeginning) +
								(numofartpoints_uptonow - numofartpointsatbeginning);


					// clique only vertex detection and removal
					if (isReduced == 0) {
#ifdef REDUCTION_STEPS
						int a = *numof_removed_edges;
#endif
						//把 -1 從 padj(csrE) 中 移除
						remove_minus_ones_in_graph (nVtx, &num_edges, pxadj, padj);


						//side vertex 
						if (!dc && (updatedlen > 1)) {
							util::timestamp tf;
							int cov_i = 0;
							vertex* clique_only_v = new vertex [len];
							memset (clique_only_v, -1, sizeof(vertex) * len);
							
							// clique-only detection
							find_clique_only_vertices (clique_only_v, &cov_i, component, len,(*pxadj),(*padj));
							util::timestamp tg;
							if (Try >= THROW_AWAY) {
								*clique_det += (tg - tf);
							}


							printf("[Side vertex][remove]\n");
							double totalw_of_covs = 0;
							// clique-only vertex removal
							remove_covs (len, component, cov_i, clique_only_v, reversecomp, &bs, nVtx, weight,
									(*pxadj), (*padj), bfsorder, endpred, level, sigma, Pred, delta, (*bc),
									numof_removed_edges, idv_track, identical_sets_c, identical_sets, next_idvset_id,
									&totalw_of_covs, total_weights_of_each_comp, comp_ids_of_each_v, CCs, ff);
//							printf("totalw_of_covs: %lf\n", totalw_of_covs);
							printf("[Side vertex][remove][done]\n");

							util::timestamp th;
							if (Try >= THROW_AWAY) {
								*clique_rem += (th - tg);
							}

							//degree-1 removal
							{
								util::timestamp ta;
								if (!dd) {
									remove_degree_1s (nVtx, component, reversecomp, &bs, weight, (*pxadj), (*padj),
											(*ptadj), idv_track, identical_sets_c, identical_sets, (*bc),
											numof_removed_edges, tmark, tbfsorder, total_weights_of_each_comp,
											comp_ids_of_each_v, CCs, ff);
								}
								util::timestamp tb;
								if (Try >= THROW_AWAY) {
									*deg1_rem += (tb - ta);
								}
							}
							delete[] clique_only_v;
						}
#ifdef REDUCTION_STEPS
						int b = *numof_removed_edges;
						printf("at COV REMOVAL, %d num of edges are removed (as a result of deg1)\n", b-a);
#endif
					}

#ifdef DEBUG
					graph_check((*pxadj), (*padj), nVtx);
					for (int i = 0; i < nVtx; i++) {
						printf("%d [w:%lf] : ",i+1, weight[i]);
						for (myindex j =(*pxadj)[i]; j <(*pxadj)[i+1]; j++) {
							int v = (*padj)[j];
							printf("%d ",v+1);
						}
						printf("\n");
					}
#endif
					updatedlen = 0;
					for (int i = 0; i < len; i++) {
						if (component[i] != -1)
							updatedlen++;
					}
#ifdef BCCOMP_DBG
					printf("updatedlen:%d\n",updatedlen);
#endif

					int numofremovededgesatend = *numof_removed_edges;
					int numofartpointsatend = *numof_art_points;

					isReduced = (numofremovededgesatend - numofremovededgesatbeginning) +
							(numofartpointsatend - numofartpointsatbeginning);

				}
#ifdef DEBUG
				// Graph Check : to be removed
				for (int i = 0; i < nVtx; i++) {
					printf("vertex %d : ", i+1);
					for (myindex j =(*pxadj)[i]; j <(*pxadj)[i+1]; j++) {
						int v = (*padj)[j];
						if (v != -1) {
							printf("%d ", v+1);
							int flag = 0;
							for (myindex k = (*pxadj)[v]; k < (*pxadj)[v+1]; k++) {
								if ((*padj)[k] == i) {
									flag = 1;
									break;
								}
							}
							if (flag == 0) {
								printf("%d and %d do not match each other, remember that nVtx was %d\n",
										i+1,v+1,startingnVtx);
								break;
							}
						}
					}
					printf("\n");
				}
#endif

#ifdef BCCOMP_DBG
				for (int i = 0; i < idv_track_size; i++) {
					if (idv_track[i] != -1)
						printf("%d is a idv repr\n",i+1);
				}

				for (int i = 0; i < next_idvset_id; i++) {
					if (identical_sets_c[i] > 0) {
						printf("%d-th idv set: \n",i);
						for (int j = 0; j < identical_sets_c[i]; j++) {
							printf("identical_sets[%d][%d] : %d, %lf\n", i, j,
									identical_sets[i][j].id+1, identical_sets[i][j].weight);
						}
					}
				}
#endif

				if (num_edges > sizeofnewadjorpred) {
					Pred = (vertex *) myre_alloc (Pred, (num_edges) * sizeof(vertex));
					newadj = (vertex *) myre_alloc (newadj, (num_edges) * sizeof(vertex));
					sizeofnewadjorpred = num_edges;
				}

				if (nVtx + 1 > sizeofnewxadj) {
					newxadj = (vertex *) myre_alloc (newxadj, (nVtx + 1) * sizeof(vertex));
					sizeofnewxadj = nVtx + 1;
				}

				/**
				 * @brief
				 * Reorder : 只是把相近的node都放在一起，放在pxadj(csrV), padj(csrE)
				*/
				// Reordering wrt BFS order
				util::timestamp tl;
				
				reorder_graph (nVtx, len, component, ordered_comp, reverse_ordered_comp,
							weight, ordered_weight, pxadj, padj, newxadj, newadj,
							all_graphs_xadj, &agx, num_comp, &k, &c, biggest_cc_after, ff, ordered_ff);

				util::timestamp tm;
				if (Try >= THROW_AWAY) {
					*bfsorder_t += (tm - tl);
				}

				Zoltan_Bucket_Free(&bs);
			}
		}

		free (total_weights_of_each_comp);
		free (comp_ids_of_each_v);

		util::timestamp tn;
		/**
		 * @todo
		 * 還不知道 cardinality 是要幹嘛用的
		 * cardinality 	: 是這個 identical set 總共有多少個 identical nodes
		 * total_weight : 是這個 identical set 總共代表了多少 nodes，因為某些 identical nodes 在被移除時，可能已有先被當成 representative (AP, d1, side)
		*/
		// construct ordered_cardinality
		for (int i = 0; i < size4; i++) {
			if (ordered_comp[i] != -1) {
				if (idv_track[ordered_comp[i]] == -1) {
					ordered_cardinality[i].cardinality = 1;
					ordered_cardinality[i].total_weight = ordered_weight[i];
					ordered_cardinality[i].total_ff = ordered_ff[i];
				}
				else {
					assert((identical_sets_c[idv_track[ordered_comp[i]]] - 1) > 1);
					ordered_cardinality[i].cardinality = (identical_sets_c[idv_track[ordered_comp[i]]] - 1);
					ordered_cardinality[i].total_weight = identical_sets[idv_track[ordered_comp[i]]][0].weight;
					ordered_cardinality[i].total_ff = identical_sets[idv_track[ordered_comp[i]]][0].idv_ff;
				}
			}
		}
		// reordered adjacency lists are sorted
		for (int ii = 0; ii < *num_comp; ii++) {
			for (myindex i = all_graphs_xadj[ii]; i < all_graphs_xadj[ii+1]; i++) {
				qsort(newadj + newxadj[i], newxadj[i+1] - newxadj[i], sizeof(vertex), compare);
			}
		}
		util::timestamp to;
		if (Try >= THROW_AWAY) {
			*bfsorder_t += (to - tn);
		}

		util::timestamp t2;

#ifdef BCCOMP_DBG
		printf("BEFORE ALL SOURCES\n");

#ifdef DEBUG
		for (vertex source=0; source < nVtx; source++) {
#else
		for (vertex source=0; source < initnVtx; source++) {
#endif
			printf("bc[%d] : %lf\n",source+1,(*bc)[source]);

		}
#endif

#ifndef JUST_PREPROC
		// REAL BC COMPUTATION
		{
			/**
			 * @brief
			 * 每個component去計算BC
			*/
			for (int ii = 0; ii < (*num_comp); ii++) {
				int start = all_graphs_xadj[ii];
				int end = all_graphs_xadj[ii+1];

#ifdef PARTIAL_BC
				end = start + (int)(partial_bc_factor * (double)(end - start));
#endif

				if ((end - start) > 1) {
					int kernel = select_kernel (start, end, ordered_weight, ordered_cardinality);
					if (kernel == 0) // no weight, no cardinality
						compute_bc_base(start, end, ordered_comp, newxadj, newadj, bfsorder,
								endpred, level, sigma, Pred, delta, (*bc), phase1time, phase2time, CCs, ff);
					else if (kernel == 1) // only weight //[Not Done]
						compute_bc_weight (start, end, ordered_comp, ordered_weight, newxadj,
								newadj, bfsorder, endpred, level, sigma,
								Pred, delta, (*bc), phase1time, phase2time);
					else if (kernel == 2) // only cardinality
						compute_bc_card (start, end, ordered_comp, newxadj, newadj, bfsorder, endpred, level, sigma,
								Pred, delta, (*bc), ordered_cardinality, phase1time, phase2time, CCs);
					else if (kernel == 3) // weight and cardinality //[Not Done]
						compute_bc_weight_card(start, end, ordered_comp, ordered_weight, newxadj, newadj, bfsorder,
								endpred, level, sigma, Pred, delta, (*bc), ordered_cardinality, phase1time, phase2time);
				}
			}
		}


#ifdef BCCOMP_DBG
		printf("AFTER ALL SOURCES\n");
#ifdef DEBUG
		for (vertex source=0; source<nVtx; source++) {
#else
		for (vertex source=0; source<initnVtx; source++) {
#endif
			printf("bc[%d] : %lf\n",source+1,(*bc)[source]);

		}
#endif
		util::timestamp t_4190;

		// bc's of idv vertices are adjusted, update AB
		double* total_comp_weights_of_each_comp = (double *) malloc(sizeof(double) * nVtx);
		int* comp_ids = (int *) calloc(nVtx, sizeof(int));
		comp_no = 1;

		/**
		 * 跑過每個 identical set，
		 * next_idvset_id : idvset的長度
		 * 
		 * 下面這個 for 迴圈 是把每個 component 的 total weight 算好
		*/
		for (int i = 0; i < next_idvset_id; i++) {
			if (identical_sets_c[i] > 0) {
				double total_weight = 0;

				// u 是這個 identical set 還在 graph 裡面的 node，其他的 identical node 都被壓掉了
				vertex u = identical_sets[i][1].id;

				if (comp_ids[u] != 0) {
					continue;
				}

				memset (tmark, 0, sizeof(int) * nVtx);
				int endofbfsorder = 1;
				tbfsorder[0] = u;
				int cur = 0;
				tmark[u] = 1;
				while (cur != endofbfsorder) {
					vertex v = tbfsorder[cur];
					comp_ids[v] = comp_no;
					if (idv_track[v] == -1)
						total_weight += weight[v];
					else {
						int idx_of_v = idv_track[v];
						total_weight += identical_sets[idx_of_v][0].weight;
					}
					for (myindex j = (*pxadj)[v]; j < (*pxadj)[v + 1]; j++) {
						vertex w = (*padj)[j];
						if ((w != -1) && (tmark[w] == 0)) {
							tbfsorder[endofbfsorder++] = w;
							tmark[w] = 1;
						}
					}
					cur++;
				}
				total_comp_weights_of_each_comp[comp_no] = total_weight;
				comp_no++;
			}
		}


		/**
		 * total_comp_weights_of_each_comp : 現在已取得每個 component 的 weight
		 * 
		 * 下面這個 for迴圈 在處理 被壓掉的 identical vertex，用 identical vertex 的代表點 把 bc 共享給 identical vertex
		*/
		for (int i = 0; i < next_idvset_id; i++) {
			if (identical_sets_c[i] > 0) {

				double total_weight_of_the_graph = total_comp_weights_of_each_comp[comp_ids[identical_sets[i][1].id]];

				Betweenness bc_of_idv_repr = (*bc)[identical_sets[i][1].id];
#ifdef BCCOMP_DBG
				printf("bc of idv-%d is %lf in update AB\n",identical_sets[i][1].id + 1, bc_of_idv_repr );
#endif
				double w_of_idv_repr = identical_sets[i][1].weight;
				for (int j = 2; j < identical_sets_c[i]; j++) {
					double tobe_added = 0;
					if ((total_weight_of_the_graph - identical_sets[i][0].weight) <= 0)
						tobe_added = bc_of_idv_repr;
					else
						tobe_added = bc_of_idv_repr + ((identical_sets[i][j].weight - w_of_idv_repr) *
								(total_weight_of_the_graph - identical_sets[i][0].weight));
					(*bc)[identical_sets[i][j].id] += tobe_added;
#ifdef BCCOMP_DBG
					printf("%lf is added to bc[%d] at update AB, btw total_weight_of_the_graph is %lf\n",
							tobe_added, identical_sets[i][j].id + 1,
							(total_weight_of_the_graph - identical_sets[i][0].weight));
#endif
				}
			}
		}

		free (total_comp_weights_of_each_comp);
		free (comp_ids);

		// bc's of copied articulation vertices are transferred to associated real articulation vertex,
		// must be done at the end. Notice the reverse traversal for the sake of correct bc calculation
		/**
		 * 對那些新的 AP分身
		*/
		for (int i = nVtx - 1; i >= startingnVtx; i--) {
			int u = art_track[i - startingnVtx];
			if (u != -1) {
				(*bc)[u] += (*bc)[i];
#ifdef BCCOMP_DBG
				printf("%lf from copy artp: %d is transferred to orig artp: %d\n", (*bc)[i], i+1, u+1);
#endif
			}
		}

		util::timestamp t_4191;
		if (Try >= THROW_AWAY) {
			*preproc += (t_4191 - t_4190);
		}

#endif

#ifdef BCCOMP_DBG
		printf("WEIGHTS:\n");
		for (vertex source=0; source<initnVtx; source++) {
			printf("weight[%d] : %lf\n",source+1, weight[source]);
		}
#endif

		free(newxadj);
		free(newadj);
		free(bfsorder);
		free(Pred);
		free(endpred);
		free(level);
		free(mark);
		free(tmark);
		free(tbfsorder);
		free(weight);
		free(sigma);
		free(delta);
		free(ordered_comp );
		free(component);
		free(reverse_ordered_comp);
		free(all_graphs_xadj);
		free(all_graphs_len);
		free(bridges);
		free(nd);
		free(l);
		free(h);
		free(reversecomp);
		free(labels);
		free(ordered_weight);
		free(ordered_cardinality);
		free(art_points);
		free(art_track);
		free(idv_track);
		free(identical_sets_c);
		free(identical_sets_sz);
		for (int i = 0; i < idv_sets_size; i++)
			free(identical_sets[i]);
		free(identical_sets);

		util::timestamp t3;
		int count = 0;

#ifdef CC_PROFILE
		int numof_CCs_ate= 0;
		int maxCC_nvtx_ate = 0;
		int maxCC_nedge_ate = 0;
		int total_nvtx_ate = 0;
		int total_nedge_ate = 0;
		printf("---------CC INFORMATION AT THE END---------\n");
		extract_cc_info (nVtx, (*pxadj), (*padj), &numof_CCs_ate, &maxCC_nvtx_ate, &maxCC_nedge_ate,
				&total_nvtx_ate, &total_nedge_ate);
		printf("-------------------------------------------------\n\n\n\n\n");
#endif

#ifndef PARTIAL_BC
#ifndef JUST_PREPROC
#ifdef DEBUG
		for (vertex source=0; source<nVtx; source++) {
#else
		for (vertex source=0; source<initnVtx; source++) {
#endif

#ifdef PRINT_BC_VALUES
			printf("bc[%d] : %lf\n",source+1,(*bc)[source]);
#endif

				if ((*bc)[source] == 0)
					count++;
		}
#endif
#endif

		if (Try >= THROW_AWAY) {
			*preproc += (t2-t1);
			*pout += ((t3-t1));
		}

#ifdef CC_PROFILE
		printf("\nAt the beginning,\tnumOfCCs: %d,\tmaxCC_nvtx: %d,\tmaxCC_nedge: %d,\ttotal_nvtx: %d,\ttotal_nedge: %d\n",
				numof_CCs_atb, maxCC_nvtx_atb, maxCC_nedge_atb, total_nvtx_atb, total_nedge_atb);
		printf  ("At the end,      \tnumOfCCs: %d,\tmaxCC_nvtx: %d,\tmaxCC_nedge: %d,\ttotal_nvtx: %d,\ttotal_nedge: %d\n\n",
				numof_CCs_ate, maxCC_nvtx_ate, maxCC_nedge_ate, total_nvtx_ate, total_nedge_ate);
#endif

	}

	return cnt;
}

