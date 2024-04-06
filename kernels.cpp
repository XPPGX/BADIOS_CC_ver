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

/**
 * @brief
 * anyweighted	: 如果這個 component 的某些 node 被當成 representative (AP, d1, side)，則這個變數需要被設置成 true
 * anycarded 	: 如果這個 component 的某些 node 是 identical node
*/
int select_kernel (int start, int end, double* ordered_weight, card_info* ordered_cardinality) {
	bool anyweighted = false;
	bool anycarded = false;
	for (vertex i = start; i < end; i++) {
		if (!anyweighted && (ordered_weight[i] > 1))
			anyweighted = true;
		if (!anycarded && (ordered_cardinality[i].cardinality > 1))
			anycarded = true;
		if (anyweighted && anycarded)
			break;
	}

	if (!anyweighted && !anycarded) // base
		return 0;
	else if (anyweighted && !anycarded) // only weight
		return 1;
	else if (!anyweighted && anycarded) // only cardinality
		return 2;
	else if (anyweighted && anycarded) // weight and cardinality
		return 3;

	return -1;
}


void compute_bc_weight_card (int start, int end, vertex*  ordered_comp, double*  ordered_weight, vertex*  newxadj,
		vertex*  newadj, vertex*  bfsorder, int*  endpred, int*  level, pathnumber*  sigma,
		vertex*  Pred, Betweenness*  delta, Betweenness*  bc, card_info* ordered_cardinality,
		util::timestamp& phase1time, util::timestamp& phase2time) {

	char *p1 = (char*) &phase1time;
	char *p2 = (char*) &phase2time;

	util::timestamp *phase1 = (util::timestamp*) p1;
	util::timestamp *phase2 = (util::timestamp*) p2;

	int len = end - start;
	for (vertex source = 0; source < len; source++) {

		assert(ordered_cardinality[start + source].cardinality >= 1);
		// if source is not idv
		if (ordered_cardinality[start + source].cardinality == 1) {
			util::timestamp t1;
			int endofbfsorder = 1;
			bfsorder[0] = source;

			for (int i = start; i < end; i++)
				endpred[i - start] = newxadj[i];

			for (int i = 0; i < len; i++)
				level[i] = -2;
			level[source] = 0;

			for (int i = 0; i < len; i++)
				sigma[i] = 0;
			sigma[source] = 1;

			//step 1: build shortest path graph
			int cur = 0;
			while (cur != endofbfsorder) {
				vertex v = bfsorder[cur];
				assert (level[v] >= 0);
				for (myindex j = newxadj[start + v]; j < newxadj[start + v + 1]; j++) {
					vertex w = newadj[j];
					if (level[w] < 0) {
						level[w] = level[v]+1;
						bfsorder[endofbfsorder++] = w;
					}
					if (level[w] == level[v]+1) {
						sigma[w] += sigma[v] * (ordered_cardinality[start+v].cardinality);
					}
					else if (level[w] == level[v] - 1) {
						Pred[endpred[v]++] = w;
					}
				}
				cur++;
			}

			for (int i = 0; i <  len; i++) {
				delta[i] = 0.;
			}

			for (int i = 0; i < len; i++) {
				delta[i] += (ordered_weight[start + i] - 1); // effect of dependents on real node
			}

			//step 2: compute betweenness
			util::timestamp t2;
			for (int i = endofbfsorder - 1; i > 0; i--) {
				vertex w = bfsorder[i];
				if (ordered_cardinality[start + w].cardinality == 1) {
					for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
						vertex v = Pred[j];
						delta[v] += sigma[v] * (1+delta[w]) / sigma[w];
					}
				}
				else { // phase-2 hack for handling idv vertices
					for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
						vertex v = Pred[j];
						delta[v] += sigma[v] * ((ordered_cardinality[start + w].cardinality * (delta[w] + 1 - ordered_weight[start + w]))
								+ ordered_cardinality[start + w].total_weight) / sigma[w];
					}
				}

				bc[ordered_comp[start + w]] += ordered_weight[start + source] * delta[w]; // effect of source and source's dependents on w

#ifdef BCCOMP_DBG
				// printf("source %d adds bc[%d]: %lf\n",ordered_comp[source]+1,ordered_comp[w]+1,ordered_weight[start + source] * delta[w]);
#endif

			}
			util::timestamp t3;
			*phase1 += (t2-t1);
			*phase2 += (t3-t2);
		}
		else {
			util::timestamp t1;
			int endofbfsorder = 1;
			bfsorder[0] = source;

			for (int i = start; i < end; i++)
				endpred[i - start] = newxadj[i];

			for (int i = 0; i < len; i++)
				level[i] = -2;
			level[source] = 0;

			for (int i = 0; i < len; i++)
				sigma[i] = 0;
			sigma[source] = 1;

			//step 1: build shortest path graph
			int cur = 0;
			while (cur != endofbfsorder) {
				vertex v = bfsorder[cur];
				assert (level[v] >= 0);
				for (myindex j = newxadj[start + v]; j < newxadj[start + v + 1]; j++) {
					vertex w = newadj[j];
					if (level[w] < 0) {
						level[w] = level[v]+1;
						bfsorder[endofbfsorder++] = w;
					}
					if (level[w] == level[v]+1) { // phase-1 hack for handling idv vertices
						if ((ordered_cardinality[start + v].cardinality > 1) && (v != source)) {
							sigma[w] += sigma[v] * ordered_cardinality[start + v].cardinality;
						}
						else
							sigma[w] += sigma[v];
					}
					else if (level[w] == level[v] - 1) {
						Pred[endpred[v]++] = w;
					}
				}
				cur++;
			}

			for (int i = 0; i < len; i++) {
				delta[i] = 0.;
			}

			for (int i = 0; i < len; i++) {
				delta[i] += (ordered_weight[start + i] - 1); // effect of dependents on real node
			}

			//step 2: compute betweenness
			util::timestamp t2;
			for (int i = endofbfsorder - 1; i > 0; i--) {
				vertex w = bfsorder[i];
				if ((ordered_cardinality[start + w].cardinality > 1) && (w != source)) { // phase-2 hack for handling idv vertices
					for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
						vertex v = Pred[j];
						delta[v] += sigma[v] * ((ordered_cardinality[start + w].cardinality * (delta[w] + 1 - ordered_weight[start + w]))
								+ ordered_cardinality[start + w].total_weight) / sigma[w];
					}
				}
				else {
					for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
						vertex v = Pred[j];
						delta[v] += sigma[v] * (1+delta[w]) / sigma[w];
					}
				}

				// phase-2 hack for handling idv vertices
				bc[ordered_comp[start + w]] += ordered_cardinality[start + source].total_weight * delta[w]; // effect of source and source's dependents on w
#ifdef BCCOMP_DBG
				 printf("source %d adds bc[%d]: %lf\n",ordered_comp[source]+1,ordered_comp[w]+1,ordered_weight[start + source] * delta[w]);
#endif

			}
			util::timestamp t3;
			*phase1 += (t2-t1);
			*phase2 += (t3-t2);
		}
	}
}


void compute_bc_weight (int start, int end, vertex* ordered_comp, double* ordered_weight, vertex* newxadj,
		vertex* newadj, vertex* bfsorder, int* endpred, int* level, pathnumber* sigma,
		vertex* Pred, Betweenness* delta, Betweenness* bc, util::timestamp& phase1time, util::timestamp& phase2time) {

	char *p1 = (char*) &phase1time;
	char *p2 = (char*) &phase2time;

	util::timestamp *phase1 = (util::timestamp*) p1;
	util::timestamp *phase2 = (util::timestamp*) p2;

	int len = end - start;
	for (vertex source = 0; source < len; source++) {

		// if source is not idv
		util::timestamp t1;
		int endofbfsorder = 1;
		bfsorder[0] = source;

		for (int i = start; i < end; i++)
			endpred[i - start] = newxadj[i];

		for (int i = 0; i < len; i++)
			level[i] = -2;
		level[source] = 0;

		for (int i = 0; i < len; i++)
			sigma[i] = 0;
		sigma[source] = 1;

		//step 1: build shortest path graph
		int cur = 0;
		while (cur != endofbfsorder) {
			vertex v = bfsorder[cur];
			assert (level[v] >= 0);
			for (myindex j = newxadj[start + v]; j < newxadj[start + v + 1]; j++) {
				vertex w = newadj[j];
				if (level[w] < 0) {
					level[w] = level[v]+1;
					bfsorder[endofbfsorder++] = w;
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

		for (int i = 0; i <  len; i++) {
			delta[i] = 0.;
		}

		for (int i = 0; i < len; i++) {
			delta[i] += (ordered_weight[start + i] - 1); // effect of dependents on real node
		}

		//step 2: compute betweenness
		util::timestamp t2;
		for (int i = endofbfsorder - 1; i > 0; i--) {
			vertex w = bfsorder[i];
			for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
				vertex v = Pred[j];
				delta[v] += sigma[v] * (1+delta[w]) / sigma[w];
			}

			bc[ordered_comp[start + w]] += ordered_weight[start + source] * delta[w]; // effect of source and source's dependents on w

#ifdef BCCOMP_DBG
			 printf("source %d adds bc[%d]: %lf\n",ordered_comp[source]+1,ordered_comp[w]+1,ordered_weight[start + source] * delta[w]);
#endif

		}
		util::timestamp t3;
		*phase1 += (t2-t1);
		*phase2 += (t3-t2);
	}
}


void compute_bc_card (int start, int end, vertex*  ordered_comp, vertex*  newxadj, vertex*  newadj,
		vertex*  bfsorder, int*  endpred, int*  level, pathnumber*  sigma, vertex*  Pred, Betweenness*  delta,
		Betweenness*  bc, card_info* ordered_cardinality, util::timestamp& phase1time, util::timestamp& phase2time, double* CCs) {

	char *p1 = (char*) &phase1time;
	char *p2 = (char*) &phase2time;

	util::timestamp *phase1 = (util::timestamp*) p1;
	util::timestamp *phase2 = (util::timestamp*) p2;

	int len = end - start;
	for (vertex source = 0; source < len; source++) {

		// if source is not idv
		if (ordered_cardinality[start + source].cardinality == 1) {
			util::timestamp t1;
			int endofbfsorder = 1;
			bfsorder[0] = source;

			// for (int i = start; i < end; i++)
			// 	endpred[i - start] = newxadj[i];

			for (int i = 0; i < len; i++)
				level[i] = -2;
			level[source] = 0;

			// for (int i = 0; i < len; i++)
			// 	sigma[i] = 0;
			// sigma[source] = 1;

			//step 1: build shortest path graph
			int cur = 0;
			while (cur != endofbfsorder) {
				vertex v = bfsorder[cur];
				assert (level[v] >= 0);
				for (myindex j = newxadj[start + v]; j < newxadj[start + v + 1]; j++) {
					vertex w = newadj[j];
					if (level[w] < 0) {
						level[w] = level[v]+1;
						bfsorder[endofbfsorder++] = w;
						
						#pragma region CC

						if(ordered_cardinality[start + w].cardinality == 1){
							CCs[source] += ordered_cardinality[start + w].total_ff + level[w];
						}
						else if(ordered_cardinality[start + w].cardinality > 1){
							//wait reorder ff
							CCs[source] += ordered_cardinality[start + w].total_ff + level[w] * ordered_cardinality[start + w].total_weight;
						}

						#pragma endregion //CC
					}

					// if (level[w] == level[v]+1) {
					// 	sigma[w] += sigma[v] * ordered_cardinality[start + v].cardinality;
					// }
					// else if (level[w] == level[v] - 1) {
					// 	Pred[endpred[v]++] = w;
					// }
				}
				cur++;
			}

			// for (int i = 0; i <  len; i++) {
			// 	delta[i] = 0.;
			// }

			//step 2: compute betweenness
			util::timestamp t2;
// 			for (int i = endofbfsorder - 1; i > 0; i--) {
// 				vertex w = bfsorder[i];
// 				if (ordered_cardinality[start + w].cardinality == 1) {
// 					for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
// 						vertex v = Pred[j];
// 						delta[v] += sigma[v] * (1+delta[w])/sigma[w];
// 					}
// 				}
// 				else { // phase-2 hack for handling idv vertices
// 					for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
// 						vertex v = Pred[j];
// 						delta[v] += sigma[v] * ((ordered_cardinality[start + w].cardinality * delta[w]) +
// 								(ordered_cardinality[start + w].total_weight)) / sigma[w];
// 					}
// 				}

// 				bc[ordered_comp[start + w]] += delta[w];
// #ifdef BCCOMP_DBG
// 				printf("source %d adds bc[%d]: %lf\n",ordered_comp[source]+1,ordered_comp[start+w]+1, delta[w]);
// #endif

// 			}
			util::timestamp t3;
			*phase1 += (t2-t1);
			*phase2 += (t3-t2);
		}
		else {
			util::timestamp t1;
			int endofbfsorder = 1;
			bfsorder[0] = source;

			// for (int i = start; i < end; i++)
			// 	endpred[i - start] = newxadj[i];

			for (int i = 0; i < len; i++)
				level[i] = -2;
			level[source] = 0;

			// for (int i = 0; i < len; i++)
			// 	sigma[i] = 0;
			// sigma[source] = 1;

			//step 1: build shortest path graph
			int cur = 0;
			while (cur != endofbfsorder) {
				vertex v = bfsorder[cur];
				assert (level[v] >= 0);
				for (myindex j = newxadj[start + v]; j < newxadj[start + v + 1]; j++) {
					vertex w = newadj[j];
					if (level[w] < 0) {
						level[w] = level[v]+1;
						bfsorder[endofbfsorder++] = w;

						#pragma region CC
						
						if(ordered_cardinality[start + w].cardinality == 1){
							CCs[source] += ordered_cardinality[start + w].total_ff + level[w];
							// CCs[source] += ff[w] + level[w];
						}
						else if(ordered_cardinality[start + w].cardinality > 1){
							//wait reorder ff
							CCs[source] += ordered_cardinality[start + w].total_ff + level[w] * ordered_cardinality[start + w].total_weight;
						}

						#pragma endregion //CC
					}


					// if (level[w] == level[v]+1) { // phase-1 hack for handling idv vertices
					// 	if ((ordered_cardinality[start + v].cardinality > 1) && (v != source)) {
					// 		sigma[w] += sigma[v] * ordered_cardinality[start + v].cardinality;
					// 	}
					// 	else
					// 		sigma[w] += sigma[v];
					// }
					// else if (level[w] == level[v] - 1) {
					// 	Pred[endpred[v]++] = w;
					// }
				}
				cur++;
			}

			// for (int i = 0; i < len; i++) {
			// 	delta[i] = 0.;
			// }

			//step 2: compute betweenness
			util::timestamp t2;
// 			for (int i = endofbfsorder - 1; i > 0; i--) {
// 				vertex w = bfsorder[i];
// 				if ((ordered_cardinality[start + w].cardinality > 1) && (w != source)) { // phase-2 hack for handling idv vertices
// 					for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
// 						vertex v = Pred[j];
// 						delta[v] += sigma[v] * ((ordered_cardinality[start + w].cardinality * delta[w]) +
// 								ordered_cardinality[start + w].total_weight) / sigma[w];
// 					}
// 				}
// 				else {
// 					for (myindex j = newxadj[w + start]; j < endpred[w]; j++) {
// 						vertex v = Pred[j];
// 						delta[v] += sigma[v] * (1+delta[w])/sigma[w];
// 					}
// 				}

// 				// phase-2 hack for handling idv vertices
// 				bc[ordered_comp[start + w]] += ordered_cardinality[start + source].total_weight * delta[w]; // effect of source and source's dependents on w
// #ifdef BCCOMP_DBG
// 				printf("source %d adds bc[%d]: %lf\n", ordered_comp[source]+1, ordered_comp[start+w]+1, (identical_sets_c[idx_of_source] - 1));
// #endif

// 			}
			util::timestamp t3;
			*phase1 += (t2-t1);
			*phase2 += (t3-t2);
		}
	}
}


void compute_bc_base (int start, int end, vertex* ordered_comp, vertex* newxadj,
		vertex* newadj, vertex* bfsorder, int* endpred, int* level, pathnumber* sigma,
		vertex* Pred, Betweenness* delta, Betweenness* bc, util::timestamp& phase1time, util::timestamp& phase2time, double* CCs, double* ff) {

	char *p1 = (char*) &phase1time;
	char *p2 = (char*) &phase2time;

	util::timestamp *phase1 = (util::timestamp*) p1;
	util::timestamp *phase2 = (util::timestamp*) p2;

	int len = end- start;
	for (vertex source = 0; source < len; source++) {
		util::timestamp t1;
		int endofbfsorder = 1;
		bfsorder[0] = source;

		// for (int i = start; i < end; i++)
		// 	endpred[i - start] = newxadj[i];

		for (int i = 0; i < len; i++)
			level[i] = -2;
		level[source] = 0;

		// for (int i = 0; i < len; i++)
		// 	sigma[i] = 0;
		// sigma[source] = 1;

		//step 1: build shortest path graph
		int cur = 0;
		while (cur != endofbfsorder) {
			vertex v = bfsorder[cur];
			assert (level[v] >= 0);
			for (myindex j = newxadj[start + v]; j < newxadj[start + v + 1]; j++) {
				vertex w = newadj[j];
				if (level[w] < 0) {
					level[w] = level[v]+1;
					bfsorder[endofbfsorder++] = w;
					CCs[source] += level[w];
				}
				// if (level[w] == level[v]+1) {
				// 	sigma[w] += sigma[v];
				// }
				// else if (level[w] == level[v] - 1) {
				// 	Pred[endpred[v]++] = w;
				// }
			}
			cur++;
		}

		// for (int i = 0; i <  len; i++) {
		// 	delta[i] = 0.;
		// }

		//step 2: compute betweenness
		util::timestamp t2;
		// for (int i = endofbfsorder - 1; i > 0; i--) {
		// 	vertex w = bfsorder[i];
		// 	for (int j = newxadj[w + start]; j < endpred[w]; j++) {
		// 		vertex v = Pred[j];
		// 		delta[v] += sigma[v] * (1 + delta[w])/sigma[w];
		// 	}
		// 	bc[ordered_comp[start + w]] += delta[w];
		// }
		util::timestamp t3;
		*phase1 += (t2-t1);
		*phase2 += (t3-t2);
	}
}
