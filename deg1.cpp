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


void remove_degree_1s (int nVtx, vertex* component, vertex* reversecomp, Bucket* bs, double* weight,
		int* xadj, vertex* adj, vertex *tadj, int* idv_track, int* identical_sets_c, idv_info** identical_sets,
		Betweenness* bc, int* numof_removed_edges, int* tmark, int *tbfsorder, double* total_weights_of_each_comp,
		int* comp_ids_of_each_v, double* CCs, double* ff) {
	
	double totalremw, total_w_of_comp;
	int count = 0;
	while(1) {
		// u 只是 index 而已
		int u = Zoltan_Bucket_PopMin(bs);
		if (u == -1)
			break;

		if (bs->values[u] == 0) {
			component[u] = -1;
			continue;
		}

		if (bs->values[u] > 1)
			break;

		bs->values[u] = INT_MAX;
		
		/**
		 * realu : D1 node
		*/
		int realu = component[u];
		component[u] = -1;

		/**
		 * @brief
		 * [Note]
		 * v : realu(d1Node) 的 parent
		 * 
		 * 把 realu 的 parent v 從 realu 這裡的 adj 拆掉
		 * 再去 parent v 的 adj 把 realu 拆掉
		*/
		int v = -1;
		for (myindex j = xadj[realu]; j < xadj[realu+1]; j++) {

			v = adj[j];
			if (v != -1) {
				adj[j] = -1;
				(*numof_removed_edges)++;

				for (myindex j = xadj[v]; j < xadj[v+1]; j++) {
					if (adj[j] == realu) {
						adj[j] = -1;
						break;
					}
				}
				break;
			}
		}

#ifdef BCCOMP_DBG
		printf("%d is deg1 vertex\n",realu+1);
#endif
		/**
		 * 有可能 ap 切割完之後， 生成的 AP 分身是 d1 node，
		 * 這時候 AP 分身就有可能是 identical node
		*/
		int idx_of_u = idv_track[realu];
		bool is_u_idv = (idx_of_u == -1) ? false : true;
		int idx_of_v = idv_track[v];
		bool is_v_idv = (idx_of_v == -1) ? false : true;

		totalremw = total_weights_of_each_comp[comp_ids_of_each_v[v]] -
				totalw_idv(realu, nVtx, adj, xadj, weight, idv_track, identical_sets, tmark, tbfsorder);

#ifdef BCCOMP_DBG
		printf("totalremw: %lf\n",totalremw);
#endif
		/**
		 * @todo [Wait]
		 * if, else if, else if 的 statement 都還沒看是怎麼處理的，等看完idv之後再回來看這裡
		*/
		if (is_u_idv && is_v_idv) { //case-c in papers
			// printf("\t[D1 : u_idv, v_idv]\r");
			// double added_to_repr_u = (identical_sets[idx_of_u][1].weight - 1) * (totalremw);
			// bc[identical_sets[idx_of_u][1].id] += added_to_repr_u;
#ifdef BCCOMP_DBG
			printf("%lf IS ADDED TO BC[%d]\n", added_to_repr_u, identical_sets[idx_of_u][1].id+1);
#endif
// 			for (int i = 2; i < identical_sets_c[idx_of_u]; i++) {
// 				int u_id     = identical_sets[idx_of_u][i].id;
// 				double u_weight = identical_sets[idx_of_u][i].weight;
// 				bc[u_id] += (u_weight - 1) * (totalremw);
// 				bc[u_id] -= added_to_repr_u;
// #ifdef BCCOMP_DBG
// 				printf("%lf IS ADDED TO BC[%d]\n", ((u_weight - 1) * (totalremw)) - added_to_repr_u, u_id+1);
// #endif
// 			}

			// bc[identical_sets[idx_of_v][1].id] +=
			// 		(identical_sets[idx_of_u][0].weight * (totalremw - identical_sets[idx_of_v][0].weight)) /
			// 		(identical_sets_c[idx_of_v] - 1);
#ifdef BCCOMP_DBG
			printf("%lf IS ADDED TO BC[%d]\n",
					(identical_sets[idx_of_u][0].weight * (totalremw - identical_sets[idx_of_v][0].weight)) /
					(identical_sets_c[idx_of_v] - 1), identical_sets[idx_of_v][1].id+1);
#endif
			// bc[identical_sets[idx_of_v][1].id] +=
			// 		(identical_sets[idx_of_v][1].weight - 1) * identical_sets[idx_of_u][0].weight;
#ifdef BCCOMP_DBG
			printf("%lf IS ADDED TO BC[%d]\n",
					(identical_sets[idx_of_v][1].weight - 1) * identical_sets[idx_of_u][0].weight,
					identical_sets[idx_of_v][1].id+1);
#endif
// 			for (int i = 2; i < identical_sets_c[idx_of_v]; i++) {
// 				int v_id     = identical_sets[idx_of_v][i].id;
// 				double v_weight = identical_sets[idx_of_v][i].weight;
// 				bc[v_id] += (v_weight - identical_sets[idx_of_v][1].weight) * identical_sets[idx_of_u][0].weight;
// #ifdef BCCOMP_DBG
// 				printf("%lf IS ADDED TO BC[%d]\n", (v_weight - identical_sets[idx_of_v][1].weight) *
// 						identical_sets[idx_of_u][0].weight, v_id+1);
// #endif
// 			}

			


			for (int i = 1; i < identical_sets_c[idx_of_v]; i++) {
				double tobe_added = identical_sets[idx_of_u][0].weight / (identical_sets_c[idx_of_v] - 1);
#ifdef BCCOMP_DBG
				printf("weight[%d]: %lf and then %lf is added to weight[%d]\n", identical_sets[idx_of_v][i].id+1,
						identical_sets[idx_of_v][i].weight, tobe_added, identical_sets[idx_of_v][i].id+1);
#endif
				identical_sets[idx_of_v][i].weight += tobe_added;
				// reflect..
				weight[identical_sets[idx_of_v][i].id] += tobe_added;
			}

			#pragma region CC
			ff[identical_sets[idx_of_v][0].id] += ff[identical_sets[idx_of_u][0].id] + identical_sets[idx_of_u][0].weight;
			identical_sets[idx_of_u][0].idv_ff += ff[identical_sets[idx_of_u][0].id] + identical_sets[idx_of_u][0].weight;
			#pragma endregion //CC 
			
			identical_sets[idx_of_v][0].weight += identical_sets[idx_of_u][0].weight;
#ifdef BCCOMP_DBG
			printf("total weight of idv set of %d becomes %lf\n", identical_sets[idx_of_v][0].id+1,
					identical_sets[idx_of_v][0].weight);
#endif

			// cancel the idv set list of u, so first add bc's of idv deps back..
// 			for (int i = 2; i < identical_sets_c[idx_of_u]; i++) {
// 				bc[identical_sets[idx_of_u][i].id] += bc[identical_sets[idx_of_u][1].id];
// #ifdef BCCOMP_DBG
// 				printf("%lf (coming from bc of %d) is added to bc of %d and so bc[%d]: %lf\n",
// 						bc[identical_sets[idx_of_u][1].id], identical_sets[idx_of_u][1].id+1,
// 						identical_sets[idx_of_u][i].id+1,
// 						identical_sets[idx_of_u][i].id+1, bc[identical_sets[idx_of_u][i].id]);
// #endif
// 			}

			identical_sets_c[idx_of_u] = 0;
			identical_sets[idx_of_u][0].id = -1;
			// cancel idv_track of u
			idv_track[realu] = -1;
#ifdef BCCOMP_DBG
			double totalremw = totalw_idv(v, nVtx, adj, xadj, weight, idv_track, identical_sets, tmark, tbfsorder);
			printf("totalremw becomes %lf after weight addition\n",totalremw);
#endif

		}
		else if (is_u_idv && !is_v_idv) { //case-a in papers
			// printf("\t[D1 : u_idv, v_nidv]\r");
			// bc[v] += identical_sets[idx_of_u][0].weight * (totalremw - 1);
#ifdef BCCOMP_DBG
			printf("%lf IS ADDED TO BC[%d]\n", identical_sets[idx_of_u][0].weight * (totalremw - 1), v+1);
#endif

			// cancel the idv set list of u, so first add bc's of idv deps back..
// 			for (int i = 2; i < identical_sets_c[idx_of_u]; i++) {
// 				bc[identical_sets[idx_of_u][i].id] += bc[identical_sets[idx_of_u][1].id];
// #ifdef BCCOMP_DBG
// 				printf("%lf (coming from bc of %d) is added to bc of %d and so bc[%d]: %lf\n",
// 						bc[identical_sets[idx_of_u][1].id], identical_sets[idx_of_u][1].id+1,
// 						identical_sets[idx_of_u][i].id+1,
// 						identical_sets[idx_of_u][i].id+1, bc[identical_sets[idx_of_u][i].id]);
// #endif
// 			}

// 			for (int i = 1; i < identical_sets_c[idx_of_u]; i++) {
// 				int u_id = identical_sets[idx_of_u][i].id;
// 				double u_weight = identical_sets[idx_of_u][i].weight;
// 				bc[u_id] += (u_weight - 1) * (totalremw);
// #ifdef BCCOMP_DBG
// 				printf("%lf IS ADDED TO BC[%d], BTW u_weight:%lf\n",
// 						((u_weight - 1) * (totalremw/* + identical_sets[idx_of_u][0].weight - u_weight*/))
// 						/*- added_to_repr*/, u_id+1, u_weight);
// #endif
// 			}

#ifdef BCCOMP_DBG
			printf("%lf is added to weight[%d]\n", identical_sets[idx_of_u][0].weight, v+1);
#endif

			#pragma region CC
			//只要 ff 跟 weight，bc的都可以不要
			ff[v] += ff[identical_sets[idx_of_u][0].id] + identical_sets[idx_of_u][0].weight;
			#pragma endregion //CC

			weight[v] += identical_sets[idx_of_u][0].weight;

			identical_sets_c[idx_of_u] = 0;
			identical_sets[idx_of_u][0].id = -1;
			// cancel idv_track of u
			idv_track[realu] = -1;
#ifdef BCCOMP_DBG
			double totalremw = totalw_idv(v, nVtx, adj, xadj, weight, idv_track, identical_sets, tmark, tbfsorder);
			printf("totalremw becomes %lf after weight addition\n",totalremw);
#endif
		}
		else if (!is_u_idv && is_v_idv) { //case-b in papers
			// printf("\t[D1 : u_nidv, v_idv]\r");
			// bc[realu] += (weight[realu] - 1) * (totalremw);
#ifdef BCCOMP_DBG
			printf("%lf IS ADDED TO BC[%d]\n", (weight[realu] - 1) * (totalremw), realu+1);
#endif
			// bc[identical_sets[idx_of_v][1].id] += weight[realu] * (identical_sets[idx_of_v][1].weight - 1);
#ifdef BCCOMP_DBG
			printf("%lf IS ADDED TO BC[%d]\n", weight[realu] * (identical_sets[idx_of_v][1].weight - 1),
					identical_sets[idx_of_v][1].id + 1);
#endif
			// bc[identical_sets[idx_of_v][1].id] += (weight[realu] * (totalremw - identical_sets[idx_of_v][0].weight)) /
			// 		(identical_sets_c[idx_of_v] - 1);
#ifdef BCCOMP_DBG
			printf("%lf IS ADDED TO BC[%d]\n", (weight[realu] * (totalremw - identical_sets[idx_of_v][0].weight)) /
					(identical_sets_c[idx_of_v] - 1), identical_sets[idx_of_v][1].id+1);
#endif
// 			for (int i = 2; i < identical_sets_c[idx_of_v]; i++) {
// 				int v_id = identical_sets[idx_of_v][i].id;
// 				double v_weight = identical_sets[idx_of_v][i].weight;
// 				bc[v_id] += weight[realu] * (v_weight - identical_sets[idx_of_v][1].weight);
// #ifdef BCCOMP_DBG
// 				printf("%lf IS ADDED TO BC[%d]\n", weight[realu] * (v_weight - identical_sets[idx_of_v][1].weight), v_id+1);
// #endif
// 			}
			for (int i = 1; i < identical_sets_c[idx_of_v]; i++) {
				double tobe_added = weight[realu] / (identical_sets_c[idx_of_v] - 1);
#ifdef BCCOMP_DBG
				printf("%lf is added to weight[%d]\n", tobe_added, identical_sets[idx_of_v][i].id+1);
#endif
				identical_sets[idx_of_v][i].weight += tobe_added;
				// reflect..
				weight[identical_sets[idx_of_v][i].id] += tobe_added;
			}
			
			#pragma region CC
			//只要 ff 跟 weight，bc的都可以不要
			ff[identical_sets[idx_of_v][0].id] += ff[realu] + weight[realu];
			identical_sets[idx_of_v][0].idv_ff += ff[realu] + weight[realu];
			#pragma endregion //CC

			//             weight[identical_sets[idx_of_v][1].id] = identical_sets[idx_of_v][1].weight;
			identical_sets[idx_of_v][0].weight += weight[realu];
#ifdef BCCOMP_DBG
			double totalremw = totalw_idv(v, nVtx, adj, xadj, weight, idv_track, identical_sets, tmark, tbfsorder);
			printf("totalremw becomes %lf after weight addition\n",totalremw);
#endif
		}
		else { //usual case
			// printf("\t[D1 : usual]\r");
			// bc[realu] += (weight[realu] - 1) * totalremw; // effect of u's dependents on v's component, for bc of u
#ifdef BCCOMP_DBG
			printf("source d1 adds bc[%d]: %lf\n",realu+1, (weight[realu] - 1) * totalremw);
#endif
			// bc[v] +=  weight[realu] * (totalremw - 1); // effect of u on v's component except v, for bc of v
#ifdef BCCOMP_DBG
			printf("source d1 adds bc[%d]: %lf\n",v+1, weight[realu] * (totalremw - 1));
#endif
#ifdef BCCOMP_DBG
			printf("%lf is added to weight[%d]\n", weight[realu], v+1);
#endif
			#pragma region CC

			//只要 ff 跟 weight，bc的都可以不要
			ff[v] += ff[realu] + weight[realu];
			// printf("ff[%d] = %f\n", v, ff[v]);
			// if(std::isnan(ff[v])){
			// 	printf("ff[%d] = %f\n", v, ff[v]);
			// }

			#pragma endregion //CC

			weight[v] += weight[realu]; // effect of u on v's components, for bc's of v's component
			
#ifdef BCCOMP_DBG
			double totalremw = totalw_idv(v, nVtx, adj, xadj, weight, idv_track, identical_sets, tmark, tbfsorder);
			printf("totalremw becomes %lf after weight addition\n",totalremw);
#endif
		}

		int rev_v = reversecomp[v];
		Zoltan_Bucket_DecVal(bs, rev_v);
	}

	return;
}
