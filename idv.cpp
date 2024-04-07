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


void idv_detection_and_merge(int type, int len, int* component, int* reversecomp, int* xadj,
		int* adj, int nVtx, int* idv_sets_size, int* one_set_size,
		int* idv_track_size, int* next_idvset_id, int** idv_track, int** identical_sets_c, int** identical_sets_sz,
		idv_info*** identical_sets, double* weight, Betweenness* bc,
		Bucket* bs, int* numof_removed_edges, int* numof_identical_vertices,
		util::timestamp& idvdet, util::timestamp& idvrem, double* CCs, double* ff) {

	
	/**
	 * identical_sets_c : 這個 2維陣列，被初始化為 0
	*/
	vertex hash_size = nVtx;
	char *p1 = (char*) &idvdet;
	char *p2 = (char*) &idvrem;

	util::timestamp *det = (util::timestamp*) p1;
	util::timestamp *rem = (util::timestamp*) p2;

	util::timestamp t1;

	/**
	 * neigsums 	: 	This is a hash table that contains each node u which is hashed into this bin
	 * 					by the hash value we get from the summation of neighborsID of each node.
	*/
	vector<vector<int> > neigsums;
	neigsums.resize(hash_size);

	/**
	 * 把每個 node u 的 neighborID 都加起來當作準備要hash的值 nsum，
	 * 用 nsum 把這個 node u 塞進對應的 hash table中。
	*/
	for (int i = 0; i < len; i++) {
		int nsum = 0;
		int u = component[i];
		if (u != -1) {
			for (myindex j = xadj[u]; j < xadj[u+1]; j++) {
				if (adj[j] != -1)
					nsum += adj[j];
			}
			if (nsum > 0) {
				if (type == 2)
					nsum += u; //for type-2 idv detection
				int idx = nsum % hash_size;
				assert(idx >= 0);
				assert(idx < hash_size);
				assert(neigsums[idx].size() < len);
				neigsums[idx].push_back(u);
			}
		}
	}

#ifdef VLIST
	for (int i = 0; i < hash_size; i++) {
		if (neigsums[i].size() > 0) {
			printf("type:%d, neigsums[%d]: \n",type,i);
			for (int j = 0; j < neigsums[i].size(); j++)
				printf("neigsums[%d][%d]: %d\n", i, j, neigsums[i][j]+1);
		}
	}
#endif

	/**
	 * 檢查 idv_track_size 的大小，如果 idv_track_size < nVtx(新的node size (因為有新的AP分身))，
	 * 則 realloc idv_track, 並把新的空間都 assign -1
	*/
	//myre_allocations for idv_track
	{
		
		if (nVtx > *idv_track_size) {
			(*idv_track) = (int*) myre_alloc ((*idv_track), sizeof(int) * nVtx);
			for (int i = *idv_track_size; i < nVtx; i++)
				(*idv_track)[i] = -1;
			*idv_track_size = nVtx;
		}
	}

	util::timestamp t2;
	*det += (t2 - t1);

	/**
	 * 掃過 hash table 的所有 bin
	 * 如果 neigsnums[i].size() > 1，代表這個 hash table 的 bin 裡面有2個以上的 nodes，他們的 neighbor 可能相同
	 * 
	*/
	// idv sets construction, bc adjustments
	for (int i = 0; i < hash_size; i++) {
		if (neigsums[i].size() > 1) {
			//// all insertions to already existing idv sets are done here
			for (int j = 0; j < neigsums[i].size(); j++) {
				int v1 = neigsums[i][j];
				if ((v1 != -1) && ((*idv_track)[v1] != -1)) { // 如果 v1 是 identical vertex representative
					int neig_num = 0;// used in bc calculation below
					for (myindex k = xadj[v1]; k < xadj[v1+1]; k++) {
						
						//計數 neighbor 的 數量，如果neighbor是idv，則一次就會加一堆
						int vtx_k = adj[k];
						if (vtx_k != -1) {
							
							//一開始沒有node是 identical vertex，所以都會先做 if statement
							if ((*idv_track)[vtx_k] == -1)
								neig_num++; //如果 vtx_k 不是 identical vertex representative，neig_num 就 +1
							else
								neig_num += (*identical_sets_c)[(*idv_track)[vtx_k]] - 1; //不知道這邊為啥要 - 1
						}
					}

					//idx_of_v1 : 是 identical set 的 id，不是 node ID
					int idx_of_v1 = (*idv_track)[v1];
					for (int k = 0; k < neigsums[i].size(); k++) {
						int v2 = neigsums[i][k];
						//如果 v2 跟 v1 是 identical 的
						if ((v2 != -1) && (v2 != v1) && equal_vertices(type, v1, v2, xadj, adj)) {
							(*numof_identical_vertices)++;
#ifdef BCCOMP_DBG
							printf("%d and %d are equal vertices\n",v1+1,v2+1);
#endif
							int idx_of_v2 = (*idv_track)[v2];
							bool is_v2_idv = (idx_of_v2 == -1) ? false : true;
							
							//new_weight : 準備要被壓掉的 weight 就記錄在這裡

							double new_weight;
							
							#pragma region CC
							
							//new_ff : 準備要被壓掉的 ff 就記錄在這裡
							double new_ff;

							#pragma endregion //CC
							
							/**
							 * [BC update] 因為identical vertex而更新
							 * 如果 v2 是 identical vertex :
							 * 
							 * @todo [Wait]
							 * identical_sets[idx_of_v2][1] 跟 identical_sets[idx_of_v1][0] 是甚麼
							 * [0] : 是對外，給原本不在這個identical_set的人看的 identical_set 的 總weight
							 * [1] : 有可能是紀錄 idv代表點 原本還沒成為idv_node之前，自己的weight
							*/
							if (is_v2_idv) {
								// update C is moved here, since it must be added when idv created
								// bc[(*identical_sets)[idx_of_v2][1].id] +=
								// 		((*identical_sets)[idx_of_v2][1].weight - 1) * (*identical_sets)[idx_of_v1][0].weight;
#ifdef BCCOMP_DBG
								printf("%lf IS ADDED TO bc[%d], UPDATE C\n",
										((*identical_sets)[idx_of_v2][1].weight - 1) * (*identical_sets)[idx_of_v1][0].weight,
										(*identical_sets)[idx_of_v2][1].id+1);
#endif
// 								for (int l = 2; l < (*identical_sets_c)[idx_of_v2]; l++) {
// 									bc[(*identical_sets)[idx_of_v2][l].id] +=
// 											((*identical_sets)[idx_of_v2][l].weight - (*identical_sets)[idx_of_v2][1].weight)
// 											* (*identical_sets)[idx_of_v1][0].weight;
// #ifdef BCCOMP_DBG
// 									printf("%lf IS ADDED TO bc[%d], UPDATE C\n",
// 											((*identical_sets)[idx_of_v2][l].weight - (*identical_sets)[idx_of_v2][1].weight)
// 											* (*identical_sets)[idx_of_v1][0].weight, (*identical_sets)[idx_of_v2][l].id+1);
// #endif
// 								}

								// v1 updates
								// bc[(*identical_sets)[idx_of_v1][1].id] += ((*identical_sets)[idx_of_v1][1].weight - 1) *
								// 		(*identical_sets)[idx_of_v2][0].weight;
#ifdef BCCOMP_DBG
								printf("%lf IS ADDED TO bc[%d], UPDATE C\n", ((*identical_sets)[idx_of_v1][1].weight - 1) *
										(*identical_sets)[idx_of_v2][0].weight, (*identical_sets)[idx_of_v1][1].id+1);
#endif

// 								for (int l = 2; l < (*identical_sets_c)[idx_of_v1]; l++) {
// 									bc[(*identical_sets)[idx_of_v1][l].id] +=
// 											((*identical_sets)[idx_of_v1][l].weight - (*identical_sets)[idx_of_v1][1].weight)
// 											* (*identical_sets)[idx_of_v2][0].weight;
// #ifdef BCCOMP_DBG
// 									printf("%lf IS ADDED TO bc[%d], UPDATE C\n",
// 											((*identical_sets)[idx_of_v1][l].weight - (*identical_sets)[idx_of_v1][1].weight)
// 											* (*identical_sets)[idx_of_v2][0].weight, (*identical_sets)[idx_of_v1][l].id+1);
// #endif
// 								}

								// v2 is newly binded to v1. So, 'bc of v1' amount of centrality at this point is subtracted
								// from bc of v2 and its idv array at this point. It will be compensated
								// at the end.
								// However, when v2 (which is to be merged to v1) is repr of another idv, then 'bc of v2'
								// amount of centrality is added back
								// to other elements of that idv set (i think), so do it that first.
// 								for (int l = 2; l < (*identical_sets_c)[idx_of_v2]; l++) {
// 									bc[(*identical_sets)[idx_of_v2][l].id] += bc[(*identical_sets)[idx_of_v2][1].id];
// #ifdef BCCOMP_DBG
// 									printf("%lf IS ADDED TO bc[%d]\n", bc[(*identical_sets)[idx_of_v2][1].id],
// 											(*identical_sets)[idx_of_v2][l].id+1);
// #endif
// 								}

// 								for (int l = 1; l < (*identical_sets_c)[idx_of_v2]; l++) {
// 									bc[(*identical_sets)[idx_of_v2][l].id] -= bc[(*identical_sets)[idx_of_v1][1].id];
// #ifdef BCCOMP_DBG
// 									printf("%lf IS REMOVED FROM bc[%d]\n", bc[(*identical_sets)[idx_of_v1][1].id],
// 											(*identical_sets)[idx_of_v2][l].id+1);
// #endif
// 								}


								/**
								 * @brief
								 * 更新 identical_sets 把 v2 的 idv list 都搬到 v1 的 idv list 裡面
								 * 
								 * Note :
								 * (*identical_sets_c)[idx_of_v1] 	= v1 的 identical set 的 node count(nodeNum)
								 * (*identical_sets)[idx_of_v1][0] 	= v1 的 identical set 的 第0個node
								*/
								// move v2's idv list to v1's idv list (firstly, check for realloc)
								int newsize = (*identical_sets_c)[idx_of_v1] + (*identical_sets_c)[idx_of_v2];
								if (newsize > (*identical_sets_sz)[idx_of_v1]) {
									(*identical_sets_sz)[idx_of_v1] = ENHANCE_FACTOR * newsize;
									(*identical_sets)[idx_of_v1] = (idv_info*) myre_alloc((*identical_sets)[idx_of_v1],
																		sizeof(idv_info) * (*identical_sets_sz)[idx_of_v1]);
								}

								//把identical_set[idx_of_v2]的每一個 node 的 id, weight 都搬到 identical_sets[idx_of_v1]的末端
								for (int l = 1; l < (*identical_sets_c)[idx_of_v2]; l++) {
									(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]].id = (*identical_sets)[idx_of_v2][l].id;
									
									#pragma region CC
									(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]].idv_ff = (*identical_sets)[idx_of_v2][l].idv_ff;
									#pragma endregion //CC

									(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]++].weight = (*identical_sets)[idx_of_v2][l].weight;

								}

								/**
								 * (*identical_sets)[idx_of_v2][0].id = -1，代表v2的identical set被erase
								 * (*identical_sets_c)[idx_of_v2] = 0，代表v2的identical set的長度也被設為0
								*/
								// cancel the v2's old idv set list
								(*identical_sets)[idx_of_v2][0].id = -1; 
								(*identical_sets_c)[idx_of_v2] = 0;
								// cancel idv_track
								(*idv_track)[v2] = -1;
								new_weight = (*identical_sets)[idx_of_v2][0].weight;
								
								#pragma region CC
								new_ff = (*identical_sets)[idx_of_v2][0].idv_ff;
								#pragma endregion //CC

							}
							else {
								// update C is moved here, since it must be added when idv created
								// bc[v2] += (weight[v2] - 1) * (*identical_sets)[idx_of_v1][0].weight;
#ifdef BCCOMP_DBG
								printf("%lf IS ADDED TO bc[%d], UPDATE C\n", (weight[v2] - 1) *
										(*identical_sets)[idx_of_v1][0].weight, v2+1);
#endif

								// bc[(*identical_sets)[idx_of_v1][1].id] += ((*identical_sets)[idx_of_v1][1].weight - 1) *
								// 		weight[v2];
#ifdef BCCOMP_DBG
								printf("%lf IS ADDED TO bc[%d], UPDATE C\n", ((*identical_sets)[idx_of_v1][1].weight - 1) *
										weight[v2], (*identical_sets)[idx_of_v1][1].id+1);
#endif
// 								for (int l = 2; l < (*identical_sets_c)[idx_of_v1]; l++) {
// 									bc[(*identical_sets)[idx_of_v1][l].id] +=
// 											((*identical_sets)[idx_of_v1][l].weight - (*identical_sets)[idx_of_v1][1].weight)
// 											* weight[v2];

// #ifdef BCCOMP_DBG
// 									printf("%lf IS ADDED TO bc[%d], UPDATE C\n",
// 											((*identical_sets)[idx_of_v1][l].weight - (*identical_sets)[idx_of_v1][1].weight)
// 											* weight[v2], (*identical_sets)[idx_of_v1][l].id+1);
// #endif
// 								}
								// v2 is newly binded to v1. So, 'bc of v1' amount of centrality at this point is subtracted
								// from bc of v2 at this point.
								// bc[v2] -= bc[(*identical_sets)[idx_of_v1][1].id];
#ifdef BCCOMP_DBG
								printf("%lf IS REMOVED FROM bc[%d]\n", bc[(*identical_sets)[idx_of_v1][1].id], v2+1);
#endif

								// move v2 to v1's idv list (firstly, check for realloc)
								/**
								 * 因為 v2 不是 identical vertex 所以只需要把 v2 這個點 搬到 v1 的 identical_sets[idx_of_v1]就好
								*/
								int newsize = (*identical_sets_c)[idx_of_v1] + 1;
								if (newsize > (*identical_sets_sz)[idx_of_v1]) {
									(*identical_sets_sz)[idx_of_v1] = ENHANCE_FACTOR * newsize;
//									printf("5th: will reallocate %ld bytes b\n",
//																	sizeof(idv_info) * (*identical_sets_sz)[idx_of_v1]);
									(*identical_sets)[idx_of_v1] = (idv_info*) myre_alloc((*identical_sets)[idx_of_v1],
																		sizeof(idv_info) * (*identical_sets_sz)[idx_of_v1]);
								}

								(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]].id = v2;
								
								#pragma region CC
								(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]].idv_ff = ff[v2];
								#pragma endregion //CC
								
								(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]++].weight = weight[v2];
								new_weight = weight[v2];

								#pragma region CC
								new_ff = ff[v2];
								#pragma endregion //CC
							}

							neigsums[i][k] = -1; // removed from neigsums list

							/**
							 * [BC update]
							 * @todo [Wait]
							 * 這邊也有BC更新的部分
							*/
							// adjust bc's of neig's of idvs if it's TYPE-1
// 							if (type == 1) {
// 								double existing_weight = (*identical_sets)[idx_of_v1][0].weight;
// 								double total_weight = new_weight + existing_weight;
// 								Betweenness tobe_added = ((total_weight * total_weight) -
// 										((existing_weight * existing_weight) + (new_weight * new_weight))) / neig_num;
// 								for (myindex k = xadj[v1]; k < xadj[v1+1]; k++) {
// 									int neig = adj[k];
// 									if (neig != -1) {
// 										if ((*idv_track)[neig] != -1) {
// 											int idx_of_neig = (*idv_track)[neig];
// 											bc[(*identical_sets)[idx_of_neig][1].id] += tobe_added;
// #ifdef BCCOMP_DBG
// 											printf("%lf is added to bc[%d]\n", tobe_added,
// 													(*identical_sets)[idx_of_neig][1].id+1);
// #endif
// 										}
// 										else {
// 											bc[neig] += tobe_added;
// #ifdef BCCOMP_DBG
// 											printf("%lf is added to bc[%d] at idv_merge of %d to %d\n", tobe_added,
// 													neig+1, v2+1, v1+1);
// #endif
// 										}
// 									}
// 								}
// 							}


							/**
							 * [Note]
							 * (*identical_sets)[idx_of_v1][0].weight : 可能是這個 identical set 的 總weight
							*/
							// total weight is updated
							(*identical_sets)[idx_of_v1][0].weight += new_weight;

							#pragma region CC

							if(type == 1){
								(*identical_sets)[idx_of_v1][0].idv_ff += new_ff + 2 * new_weight;
							}
							else{ //type == 2
								(*identical_sets)[idx_of_v1][0].idv_ff += new_ff + new_weight;
							}

							#pragma endregion //CC

							// REMOVE v2 and its adj edges from the graph
							int rev_v2 = reversecomp[v2];

							// edge removals
							for (myindex k = xadj[v2]; k < xadj[v2+1]; k++) {
								int y = adj[k];
								if (y != -1) {
									// edge from identical vertex to its neig is removed
									adj[k] = -1;
									// bucket value is decremented
									Zoltan_Bucket_DecVal(bs, rev_v2);
									// edge from its neig to identical vertex is removed
									int rev_y = reversecomp[y];
									for (myindex l = xadj[y]; l < xadj[y+1]; l++) {
										if (adj[l] == v2) {
											adj[l] = -1;
											// bucket value is decremented
											Zoltan_Bucket_DecVal(bs, rev_y);
											break;
										}
									}
									(*numof_removed_edges)++;
								}
							}

							// vertex removal
							component[rev_v2] = -1;
						}
					}
				}
			}

			//// new idv sets are created here. At each (*identical_sets)[i], we store pairs of (id, weight). Also
			//// (id_of_repr, total_weight) is stored as first element for calculation avoidance in future
			for (int j = 0; j < neigsums[i].size(); j++) {
				int v1 = neigsums[i][j];
				bool already_added = false;
				int idx_of_v1 = -1;
				if ((v1 != -1) && ((*idv_track)[v1] == -1)) { // if v1 is not an idv representative, it becomes so
					int neig_num = 0;
					for (myindex k = xadj[v1]; k < xadj[v1+1]; k++) {
						int vtx_k = adj[k];
						if (vtx_k != -1) {
							if ((*idv_track)[vtx_k] == -1)
								neig_num++;
							else
								neig_num += (*identical_sets_c)[(*idv_track)[vtx_k]] - 1;
						}
					}
					for (int k = j+1; k < neigsums[i].size(); k++) {
						int v2 = neigsums[i][k];
						if ((v2 != -1) && (v1 != v2) && equal_vertices(type, v1, v2, xadj, adj)) {
							(*numof_identical_vertices)++;
#ifdef BCCOMP_DBG
							printf("at type:%d, %d and %d are equal vertices\n",type, v1+1,v2+1);
#endif
							/**
							 * [Note]
							 * (*next_idvset_id) : 最新的idvset_id
							 * (*idv_track)[v1] : v1 所在的 identical set 的 id，這個 id 不是 真實的 nodeID 而是從 0 開始的 setID
							*/
							if (!already_added) {
								(*idv_track)[v1] = (*next_idvset_id)++;
								idx_of_v1 = (*idv_track)[v1];
								
								/**
								 * 這邊只是在看空間是否需要延長
								*/
								// enhance the size of identical_sets, identical_sets_c and identical_sets_sz, if needed
								int newsize = idx_of_v1;
								if (newsize >= (*idv_sets_size)) {
									int sz = (newsize * ENHANCE_FACTOR);
									*identical_sets = (idv_info**) myre_alloc((*identical_sets), sizeof(idv_info*) * sz);

									for (int i = *idv_sets_size; i < sz; i++)
										(*identical_sets)[i] = (idv_info*) malloc(sizeof(idv_info) * INIT_SIZE_FOR_IDV_SETS);

									(*identical_sets_c) = (int*) myre_alloc((*identical_sets_c), sizeof(int) * sz);
									for (int i = *idv_sets_size; i < sz; i++)
										(*identical_sets_c)[i] = 0;

									(*identical_sets_sz) = (int*) myre_alloc((*identical_sets_sz), sizeof(int) * sz);
									for (int i = *idv_sets_size; i < sz; i++)
										(*identical_sets_sz)[i] = INIT_SIZE_FOR_IDV_SETS;

									(*idv_sets_size) = sz;

								}


								/**
								 * (*identical_sets)[idx_of_v1][0].id 		: 紀錄 identical sets 的代表點ID
								 * (*identical_sets)[idx_of_v1][0].weight 	: 紀錄整個 identical sets 的 weight
								 * 
								 * (*identical_sets)[idx_of_v1][1].id		: 紀錄 identical sets 的代表點ID
								 * (*identical_sets)[idx_of_v1][1].weight	: 紀錄 identical node 的 weight
								 * 
								 * (*idnetical_sets)[idx_of_v1][2].id		: 紀錄被吃進此 identical sets 的 第一個 nodeID
								 * (*identical_sets)[idx_of_v1][2].weight	: 紀錄被吃進來的 第一個 node 的 weight
								*/
								int index = (*identical_sets_c)[idx_of_v1]; //如果是新創的 idv_set, 那 (*identical_sets_c)[idx_of_v1] 會是 0
								(*identical_sets)[idx_of_v1][index].id = v1;
								(*identical_sets)[idx_of_v1][index].weight = weight[v1];
								
								#pragma region CC
								(*identical_sets)[idx_of_v1][index].idv_ff = ff[v1];
								#pragma endregion //CC
								
								(*identical_sets_c)[idx_of_v1]++;


								////////////////////////////////////////////////////////////

								index = (*identical_sets_c)[idx_of_v1];
								(*identical_sets)[idx_of_v1][index].id = v1;
								(*identical_sets)[idx_of_v1][index].weight = weight[v1];

								#pragma region CC
								(*identical_sets)[idx_of_v1][index].idv_ff = ff[v1];
								#pragma endregion //CC

								(*identical_sets_c)[idx_of_v1]++;

								already_added = true;
							}
							// update C is moved here, since it must be added when idv created
							// bc[v2] += (weight[v2] - 1) * (*identical_sets)[idx_of_v1][0].weight;
#ifdef BCCOMP_DBG
							printf("%lf IS ADDED TO bc[%d], UPDATE C\n", (weight[v2] - 1) *
									(*identical_sets)[idx_of_v1][0].weight, v2+1);
#endif

							// bc[(*identical_sets)[idx_of_v1][1].id] += ((*identical_sets)[idx_of_v1][1].weight - 1) * weight[v2];
#ifdef BCCOMP_DBG
							printf("%lf IS ADDED TO bc[%d], UPDATE C\n", ((*identical_sets)[idx_of_v1][1].weight - 1) *
									weight[v2], (*identical_sets)[idx_of_v1][1].id+1);
#endif

// 							for (int l = 2; l < (*identical_sets_c)[idx_of_v1]; l++) {
// 								bc[(*identical_sets)[idx_of_v1][l].id] += ((*identical_sets)[idx_of_v1][l].weight -
// 										(*identical_sets)[idx_of_v1][1].weight) * weight[v2];

// #ifdef BCCOMP_DBG
// 								printf("%lf IS ADDED TO bc[%d], UPDATE C\n",
// 										((*identical_sets)[idx_of_v1][l].weight - (*identical_sets)[idx_of_v1][1].weight)
// 										* weight[v2], (*identical_sets)[idx_of_v1][l].id+1);
// #endif
// 							}

							// v2 is newly binded to v1. So, 'bc of v1' amount of centrality at this point is subtracted
							// from bc of v2 at this point.
							// bc[v2] -= bc[v1];
#ifdef BCCOMP_DBG
							printf("%lf IS REMOVED FROM bc[%d]\n", bc[v1], v2+1);
#endif

							// move v2 to v1's idv list (firstly, check for realloc)

							int newsize = (*identical_sets_c)[idx_of_v1] + 1;
							if (newsize > (*identical_sets_sz)[idx_of_v1]) {
								(*identical_sets_sz)[idx_of_v1] = ENHANCE_FACTOR * newsize;
								(*identical_sets)[idx_of_v1] = (idv_info*) myre_alloc((*identical_sets)[idx_of_v1],
																	sizeof(idv_info) * (*identical_sets_sz)[idx_of_v1]);
							}


							(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]].id = v2;
							
							#pragma region CC
							
							(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]].idv_ff = ff[v2];

							#pragma endregion //CC

							(*identical_sets)[idx_of_v1][(*identical_sets_c)[idx_of_v1]++].weight = weight[v2];

							




							neigsums[i][k] = -1; //removed from neigsums list
							// adjust bc's of neig's of idvs if it's TYPE-1
// 							if (type == 1) {
								
// 								//這裡的內容可以整個註解掉，之後再改成 CC的版本
// 								double new_weight = weight[v2];
// 								double existing_weight = (*identical_sets)[idx_of_v1][0].weight;
// 								double total_weight = new_weight + existing_weight;
// 								Betweenness tobe_added = ((total_weight * total_weight) -
// 										((existing_weight * existing_weight) + (new_weight * new_weight))) / neig_num;
// 								for (myindex k = xadj[v1]; k < xadj[v1+1]; k++) {
// 									int neig = adj[k];
// 									if (neig != -1) {
// 										if ((*idv_track)[neig] != -1) {
// 											int idx_of_neig = (*idv_track)[neig];
// 											bc[(*identical_sets)[idx_of_neig][1].id] += tobe_added;
// #ifdef BCCOMP_DBG
// 											printf("%lf IS ADDED TO bc[%d]\n", tobe_added,
// 													(*identical_sets)[idx_of_neig][1].id+1);
// #endif
// 										}
// 										else {
// 											bc[neig] += tobe_added;
// #ifdef BCCOMP_DBG
// 											printf("%lf is added to bc[%d] at new idv creation of id-repr:%d and %d\n",
// 													tobe_added, neig+1, v1+1, v2+1);
// #endif
// 										}
// 									}
// 								}
// 							}
							// total weight is updated
							(*identical_sets)[idx_of_v1][0].weight += weight[v2]; // total weight is updated
							
							#pragma region CC

							if(type == 1){
								(*identical_sets)[idx_of_v1][0].idv_ff += ff[v2] + 2 * weight[v2];
							}
							else{ //type == 2
								(*identical_sets)[idx_of_v1][0].idv_ff += ff[v2] + weight[v2];
							}
							
							#pragma endregion //CC

							// REMOVE v2 and its adj edges from the graph
							int rev_v2 = reversecomp[v2];
							// edge removals
							for (myindex k = xadj[v2]; k < xadj[v2+1]; k++) {
								int y = adj[k];
								if (y != -1) {
									// edge from identical vertex to its neig is removed
									adj[k] = -1;
									// bucket value is decremented
									Zoltan_Bucket_DecVal(bs, rev_v2);
									// edge from its neig to identical vertex is removed
									int rev_y = reversecomp[y];
									for (myindex l = xadj[y]; l < xadj[y+1]; l++) {
										if (adj[l] == v2) {
											adj[l] = -1;
											// bucket value is decremented
											Zoltan_Bucket_DecVal(bs, rev_y);
											break;
										}
									}
									(*numof_removed_edges)++;
								}
							}
							// vertex removals
							component[rev_v2] = -1;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < hash_size; i++)
		neigsums[i].clear();
	neigsums.clear();

	util::timestamp t3;
	*rem += (t3 - t2);
}



bool equal_vertices(int type, vertex v1, vertex v2, int* xadj, int* adj) {

	int len_v1 = 0, len_v2 = 0;
	bool flag = true;
	for (myindex i = xadj[v1]; i < xadj[v1+1]; i++) {
		if (adj[i] != -1)
			len_v1++;
	}
	for (myindex i = xadj[v2]; i < xadj[v2+1]; i++) {
		if (adj[i] != -1)
			len_v2++;
	}

	if (len_v1 != len_v2)
		return false;

	if (type == 1) {
		/**
		 * 直接檢查兩個nodes的neighbors是否一樣(不包含對方 : v1不包含v2, v2不包含v1)
		*/
		int equalities = 0;
		for (myindex i = xadj[v1]; i < xadj[v1+1]; i++) {
			if (adj[i] != -1) {
				bool sflag = false;
				for (myindex j = xadj[v2]; j < xadj[v2+1]; j++) {
					if (adj[i] == adj[j]) {
						equalities++;
						sflag = true;
						break;
					}
				}
				if (sflag)
					continue;
				else {
					flag = false;
					break;
				}
			}
		}
		if (flag)
			assert (equalities == len_v1);
	}
	else if (type == 2) {
		/**
		 * 先判斷 v1 的 neighbor中 是否有 v2
		 * 再判斷 v2 的 neighbor中 是否有 v1
		*/
		bool tiny_flag = false;
		for (myindex i = xadj[v1]; i < xadj[v1+1]; i++) {
			if (adj[i] == v2) {
				tiny_flag = true;
				break;
			}
		}
		tiny_flag = false;
		for (myindex i = xadj[v2]; i < xadj[v2+1]; i++) {
			if (adj[i] == v1) {
				tiny_flag = true;
				break;
			}
		}
		if (!tiny_flag)
			return false;
		
		/**
		 * 這邊才真的開雙層迴圈在檢查兩個nodes v1, v2的neighbor是否完全一樣，但跳過 adj[i] == v2的情況
		*/
		int equalities = 0;
		for (myindex i = xadj[v1]; i < xadj[v1+1]; i++) {
			if (adj[i] != -1) {
				bool sflag = false;
				for (myindex j = xadj[v2]; j < xadj[v2+1]; j++) {
					if (adj[i] == adj[j]) {
						sflag = true;
						equalities++;
						break;
					}
				}
				if (sflag || (adj[i] == v2))
					continue;
				else {
					flag = false;
					break;
				}
			}
		}
		if (flag)
			assert (equalities == (len_v1 - 1));
	}
	return flag;
}
