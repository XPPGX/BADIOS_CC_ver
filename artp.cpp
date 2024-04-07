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
 * @todo
 * 遇到有可能會新增 node 的 function，某些會改變長度的 array 需要使用 "pointer of pointer" 去傳遞 address
 * CCs 跟 ff 記得要改雙指標
 * 檢查 是否 ff array 中還有 Nan : 還有
 * 檢查 每個值的寫法是否正確
*/

#include "bc-seq-brandes.h"

/**
 * @todo [Wait]
 * 等看完identical vertex要回來看這個
*/
void articulation_point_copy (int* nVtx, int artc, int* nextvid, int* len, int* size1, int* size2, int* size4,
		vertex* art_points, vertex* art_track, int** pxadj,
		int** padj, Betweenness** bc, int* num_edges, int initnVtx, vertex* bfsorder,
		int* endpred, int* level, pathnumber* sigma, Betweenness* delta,
		int* mark, double* weight, vertex* component, vertex* ordered_comp, vertex* reversecomp,
		vertex* reverse_ordered_comp, double* ordered_weight, card_info* ordered_cardinality,
		int* all_graphs_xadj, int* all_graphs_len,
		int* labels, int* nd, int* l, int* h, Bucket* bs, int* next_idvset_id,
		int* idv_track_size, int** idv_track, int** identical_sets_c, int** identical_sets_sz,
		idv_info*** identical_sets, int* idv_sets_size, int one_set_size, int* tmark, int* tbfsorder,
		double* total_weights_of_each_comp, int* comp_ids_of_each_v, int* comp_no, double** CCs, double** ff) {

	vertex* newoldneigs = (vertex *)malloc(sizeof(vertex) * 2 * (*nVtx));// for each neig of art point, it is list of new
																		 // and old neigs for them
	int* componentid = (int*) malloc(sizeof(int) * (*nVtx));// conn component ids for all vertices around an art point

	//還不知道這個是幹嘛的
	int veryoldnVtx = (*nVtx);
	// art points are processed one-by-one
	for (int i = 0; i < artc; i++) {
		//新創建的 node 的 nodeID 從 nVtx 開始，nVtx 代表當前整個 graph 的 node 數量
		int startnextvid = *nextvid;// id for newcomer vertices
		memset (newoldneigs, -1, sizeof(int) * 2 * (*nVtx));
		for (int j = 0; j < (*nVtx); j++)
			componentid[j] = -2;
		
		//cid 為啥一開始 是 -1  = =
		int cid = -1;
		// vertices connected to an art point are labeled with a comp id
		int u = art_points[i]; //取得 AP nodeID
		double total_weight_at_the_beginning = totalw_idv(u, (*nVtx), (*padj), (*pxadj), weight, (*idv_track),
															(*identical_sets), tmark, tbfsorder);
#ifdef BCCOMP_DBG
		printf("art_p: %d\n",u+1);
#endif
		/**
		 * 這個 for 迴圈可以知道該 AP 連接到多少個 component
		 * 並且 assign 不同 component 內的 node 一個 cid
		 * 這些 cid 會從 "-1" 到 "component數量"
		 * 
		 * @todo
		 * 這裡要改成在 assign_component_ids 的時候，順便也要 return 該 component 的每個 nodes 對 AP本尊 的總距離
		 * 1. 宣告一個 array "comp_dist_from_u" : 去紀錄每個component對u的距離
		 * 
		 * 2. 最後要free(comp_dist_from_u)
		*/


		// printf("\t[AP : assign_component_ids]\n");
		#pragma region CC
		//comp_dist_from_u 的 index 索引方式 : cid + 1 就是 cid 對應的 comp_dist_from_u
		/**
		 * @todo 要記得free comp_dist_from_u
		*/
		// printf("current_count_AP = %d\n", i + 1);
		double* comp_dist_from_u = (double*)malloc(sizeof(double) * (*nVtx));
		int* dist_arr = (int*)malloc(sizeof(int) * (*nVtx));
		memset(comp_dist_from_u, 0, sizeof(double) * (*nVtx));
		for (myindex j =(*pxadj)[u]; j <(*pxadj)[u+1]; j++) {
			int v =(*padj)[j];
			if (v != -1) {
				if (componentid[v] == -2) {
					assign_component_ids (cid, v, u, componentid, (*padj), (*pxadj), (*nVtx), comp_dist_from_u, dist_arr, weight, *ff);
					cid++;
				}
			}
		}
		free(dist_arr);
		// printf("\t[AP : assign_component_ids][Done]\n");
		/**
		 * 取得 u 所在的 component 對 u 的 ff
		 * 之後要取得在comp[i] 的AP分身 看外面所有的ff，只要用 total_comp_dist_from_u - comp_dist_from_u[i]就好，
		*/
		//
		
		double total_comp_dist_from_u = 0;
		for(int i = 0 ; i < cid + 1 ; i ++){
			// printf("comp_dist_from_u[%d] = %f\n", i, comp_dist_from_u[i]);
			total_comp_dist_from_u += comp_dist_from_u[i];
		}
		// printf("total_comp_dist_from_u = %f\n", total_comp_dist_from_u);

		/**
		 * @brief
		 * 到這裡已取得
		 * 1. AP本尊 對切開後的每個 component 的 ff
		 * 2. AP本尊 對切開前整個 component 的 ff
		*/
		#pragma endregion //CC

		/**
		 * cid	== -1 的 nodes 代表在原本的 AP 所在的component 
		 * cid	!= -1 的 nodes 代表在原本的 AP 不在的component
		 * 
		 * 用 newoldneigs[] 紀錄每個node對應到的 "新AP(*nextvid + componentid[v])" 與 "原始AP(u)"
		 * 每個 vertex 用 newoldneigs 用連續的兩格記錄上述
		*/
		// new and old neigs of the vertices (neigs of art points) are recorded
		for (myindex j =(*pxadj)[u]; j < (*pxadj)[u+1]; j++) {
			int v =(*padj)[j]; // v 是 nodeID
			if (v != -1) {
				if (componentid[v] > -1) {
					newoldneigs[2*v] = *nextvid + componentid[v]; //u 的 AP分身ID
					newoldneigs[2*v + 1] = u; //u 的 AP本尊 ID
				}
			}
		}
		*nextvid += cid;// next vertex id to assign vertices is incremented by the number of total conn components
						// connected to an art point

		int maxdeg = 0;// max degree
		if (maxdeg < ((*pxadj)[u+1] - (*pxadj)[u]))
			maxdeg = (*pxadj)[u+1] - (*pxadj)[u];

		int endnextvid = *nextvid;
		int oldnVtx = (*nVtx);
		int oldnum_edges = *num_edges;
		int* counter_extxadj;
		int incr = 0; //新創建的 nodes 數量

		//用incr 紀錄新創建的 nodes 數量
		if ((*idv_track)[u] != -1) {
			//如果 AP本尊 u 是 identical vertex 的話，新增的 AP分身 會很多
			int idx_of_u = (*idv_track)[u];
			incr = cid * ((*identical_sets_c)[idx_of_u] - 1);
		}
		else {
			incr = cid;
		}

		/**
		 * 重新紀錄 pxadj(csrV), padj(csrE)
		 * 
		 * Note : 
		 * realloc之後，原本的記憶體如果可以延長，會直接延長，原本的值還會在裡面。
		 * 如果原本的記憶體不能延長，會把裡面的值複製到新的地方，然後free掉原本的空間。
		 * 
		 * pxadj 的改動 => 為 pxadj 重新申請一個長度 = nVtx + 1 的空間，代表新創立的AP，而 pxadj[AP] 值為 padj 中新的 offset
		 * 
		*/
		*nVtx += incr;// num of vertices are increased
		*num_edges += (incr) * maxdeg;// num of edges are increased //這裡新增的值不知道是不是對
		*pxadj = (int *) myre_alloc (*pxadj, sizeof(int) * ((*nVtx) + 1));//myre_allocations for xadj are done
		counter_extxadj = (int *)malloc(sizeof(int) * (incr));// counter for xadj array of newly created vertices
		
		/**
		 * 為啥這邊 offset 要記錄兩次
		*/
		for (int i = 0; i < incr; i++) {
			counter_extxadj[i] = (i) * maxdeg + oldnum_edges;
		}

		for (int i = 0; i < incr; i++) {
			//這邊 + 1 只是他nodeID的記法
			(*pxadj)[oldnVtx+i+1] = (i+1) * maxdeg + oldnum_edges;
		}

		*padj = (vertex *) myre_alloc (*padj, sizeof(vertex) * (*num_edges));//myre_allocations for adj are done
		for (myindex i = oldnum_edges; i < *num_edges; i++)
			(*padj)[i] = -1;

		//myre_allocations for used arrays, if needed of course
		#pragma region CC

		// printf("\t[AP][CC, ff realloc]\n");
		*CCs = (double*)myre_alloc(*CCs, sizeof(double) * (*nVtx));
		*ff = (double*)myre_alloc(*ff, sizeof(double) * (*nVtx));
		for(int i = startnextvid ; i < endnextvid ; i ++){
			(*CCs)[i] = 0;
			(*ff)[i] = 0;
			// printf("CC[%d] = %f, ff[%d] = %f\n", i, CCs[i], i, ff[i]);
		}
		// printf("\t[AP][CC, ff realloc][Done]\n");
		
		#pragma endregion //CC


		if ((*nVtx) > (*size1)) { //(*size1) was (2 * initnVtx) at the very beginning
			*size1 = *nVtx;
			// printf("rellaoc ff\n");
#ifdef DEBUG
			printf("some myre_alloc\n");
#endif
			mark = (int *) myre_alloc (mark, sizeof(int) * (*nVtx));
			weight = (double *) myre_alloc (weight, sizeof(double) * (*nVtx));
			component = (vertex *) myre_alloc (component, sizeof(vertex) * (*nVtx));
			reversecomp = (vertex *) myre_alloc (reversecomp, sizeof(vertex) * (*nVtx));
			art_track = (vertex *)myre_alloc (art_track, sizeof(vertex) * (*nVtx));
#ifdef DEBUG
			printf("size of art_track, new nVtx:%d\n",*nVtx);
#endif
			//因為原本的 nodes 不會是新創建的 AP, 只有新創建的 AP分身 才需要用art_track去追哪個node是 AP本尊
			for (int i = oldnVtx; i < *nVtx; i++)
				art_track[i - initnVtx] = -1;
		}

		// myre_allocations
		if ((*nVtx) > (*idv_track_size)) {
			(*idv_track) = (int *) myre_alloc ((*idv_track), (*nVtx) * sizeof(int));
			for (int i = *idv_track_size; i < *nVtx; i++)
				(*idv_track)[i] = -1;
			*idv_track_size = *nVtx;
		}

		int newsizeofidv_sets = (*next_idvset_id) + (endnextvid - startnextvid) + 1;
		if (newsizeofidv_sets >= (*idv_sets_size)) {

			int sz = (newsizeofidv_sets * ENHANCE_FACTOR);
			(*identical_sets) = (idv_info**) myre_alloc((*identical_sets), sz * sizeof(idv_info*));
			(*identical_sets_c) = (int *) myre_alloc ((*identical_sets_c), sz * sizeof(int));
			(*identical_sets_sz) = (int *) myre_alloc ((*identical_sets_sz), sz * sizeof(int));
			for (int i = *idv_sets_size; i < sz; i++) {
				(*identical_sets)[i] = (idv_info*) malloc(INIT_SIZE_FOR_IDV_SETS * sizeof(idv_info));
				(*identical_sets_c)[i] = 0;
				(*identical_sets_sz)[i] = INIT_SIZE_FOR_IDV_SETS;
			}
			*idv_sets_size = sz;
		}

		/**
		 * 這邊在更新 padj(csrE)，(原本有連到AP本尊 && 之後要連到AP分身) 的 nodes，更新他們的 padj(csrE)
		 * AP 本尊 => 斷開與 neighbor vid 的 edge
		 * AP 分身 => 建立與 neighbor vid 的 edge
		 * vid => 代表每個node的ID，每個 node 都要看自己是否需要把原本有連線的 AP本尊ID 換成 AP分身ID
		 * art_track => 紀錄AP分身的本尊是誰
		*/
		// neigbor lists are updated according to newly created vertices
		for (int ii = 0; ii < (2 * oldnVtx); ii+=2) {
			int vid = ii/2;
			int newneig = newoldneigs[ii]; //AP 分身 ID
			if (newneig != -1) {
				int oldneig = newoldneigs[ii+1]; //AP 本尊 ID
				for (myindex j =(*pxadj)[vid]; j <(*pxadj)[vid+1]; j++) {
					if ((*padj)[j] == oldneig) {
						//把原本padj(csrE)中原本是vid鄰居的AP本尊ID，換成AP分身ID
						(*padj)[j] = newneig;//fix the edge from neigbor_of_AP to new_AP_copy

						//這裡在建立 AP 分身的 padj(csrE)
						(*padj)[counter_extxadj[newneig - startnextvid]++] = vid;
						//紀錄 AP 分身 的 本尊 是誰
						art_track[newneig - initnVtx] = oldneig;//keeping track between newly created vertices and
																// their old art point ids
#ifdef DEBUG
						printf("art_track[%d]:%d\n", newneig, art_track[newneig - initnVtx]);
#endif
						break;
					}
				}
				
				//AP本尊 斷開與 vid 之間的 edge
				int u = oldneig; //AP 本尊 ID
				for (myindex j =(*pxadj)[u]; j <(*pxadj)[u+1]; j++) {
					int v = (*padj)[j];
					if (v == vid) {
						(*padj)[j] = -1;//fix the edge from oldartpoint to neigofart
						
						break;
					}
				}
			}
		}

		/**
		 * @todo [Wait]
		 * 如果 AP本尊 是idv，則這邊會先把每個 AP分身 都當成有自己的一個 identical set，
		 * 而且這個 identical set 的內容跟 AP本尊 一樣。
		*/
		// if art_point is also an idv
		if ((*idv_track)[u] != -1) {
			int idx_of_u = (*idv_track)[u];
			int nextnum = *nextvid;
			
			//處理這次的 AP本尊分割的時候，新產生的幾個 AP 分身的 weight(identical vertex version)
			for (int i = startnextvid; i < endnextvid; i++) {
				(*idv_track)[i] = (*next_idvset_id)++;
				int idx_of_i = (*idv_track)[i];
				// idv_set for newly created vertex is created
				(*identical_sets)[idx_of_i][0].id = i;
				(*identical_sets)[idx_of_i][0].weight = (*identical_sets)[idx_of_u][0].weight;
				(*identical_sets)[idx_of_i][1].id = i;
				(*identical_sets)[idx_of_i][1].weight = (*identical_sets)[idx_of_u][1].weight;

				// reflect weights
				weight[(*identical_sets)[idx_of_i][1].id] = (*identical_sets)[idx_of_u][1].weight;
				for (int j = 2; j < (*identical_sets_c)[idx_of_u]; j++) {
					// idv_set of newly created real vertex also consists of newly created vertices
					(*identical_sets)[idx_of_i][j].id = nextnum;
					// weights are adjusted
					(*identical_sets)[idx_of_i][j].weight = (*identical_sets)[idx_of_u][j].weight;
					// reflect
					weight[(*identical_sets)[idx_of_i][j].id] = (*identical_sets)[idx_of_u][j].weight;
					// new vertices in the idv_set of newly created real vertex will be mapped to their
					// correspondent in idv_set of orig artp
					art_track[nextnum - initnVtx] = (*identical_sets)[idx_of_u][j].id;
					nextnum++;
				}
				(*identical_sets_c)[idx_of_i] = (*identical_sets_c)[idx_of_u];
			}
			*nextvid = nextnum;
		}
		else {
			//處理在這次AP本尊分割的時候，新產生的幾個 AP 分身的weight (普通的點的version)
			for (int i = startnextvid; i < endnextvid; i++) {
				weight[i] = weight[u];
			}
		}


		/**
		 * @todo [Wait]
		 * 等看完 identical vertex 的處理，要回來看這個
		 * u_weight = u 所在的 component 的所有 nodes 的 weight 加總
		*/
		double u_weight = totalw_idv(u, (*nVtx), (*padj), (*pxadj), weight, (*idv_track), (*identical_sets), tmark, tbfsorder);
		if ((*idv_track)[u] == -1)
			u_weight -= weight[u]; //整個u所在的component的weight - u自己的weight
		else
			u_weight -= (*identical_sets)[(*idv_track)[u]][0].weight; //整個u所在的component的weight - u這個點所代表的 weight
		
		//AP本尊所在的 component 的 weight(不包含 AP本尊自身的weight)
		double old_comp_weight =  u_weight; //old_comp_weight is initialized to weight of old art point component
											// (except old art point)

		/**
		 * 對每個新的 AP 分身，以 AP分身 當 source 進行 BFS traverse 取得該 AP分身 所在的 component 的 weight
		 * 並加總到 old_comp_weight
		 * 
		 * @brief
		 * comp_
		 * 
		*/
		
		//comp_du_index : component dist from u index
		int comp_du_index = 1;
		for (int i = startnextvid; i < endnextvid; i++) {
			if (art_track[i - initnVtx] == u) {

				double i_weight;
				double total_weight = 0;
				memset (tmark, 0, sizeof(int) * (*nVtx));
				memset (tbfsorder, 0, sizeof(int) * (*nVtx));
				int endofbfsorder = 1;
				tbfsorder[0] = i;
				int cur = 0;
				tmark[i] = 1;

				//計算 AP分身 所在的component的 total weight(包含 AP分身)
				while (cur != endofbfsorder) {
					vertex v = tbfsorder[cur];
					comp_ids_of_each_v[v] = (*comp_no); //賦予每個traverse到的node一個component號碼
					if ((*idv_track)[v] == -1)
						total_weight += weight[v];
					else {
						int idx_of_v = (*idv_track)[v];
						total_weight += (*identical_sets)[idx_of_v][0].weight;
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
				i_weight = total_weight;

				(*comp_no)++;

				//i_weight = 該AP分身所在的component的total_weight - AP分身自身的weight
				if ((*idv_track)[i] == -1)
					i_weight -= weight[i];
				else
					i_weight -= (*identical_sets)[(*idv_track)[i]][0].weight;
				


				#pragma region CC
				
				//不該用 i 去索引 comp_dist_from_u，應該用別的
				(*ff)[i] += (*ff)[u] + (total_comp_dist_from_u - comp_dist_from_u[comp_du_index]);
				comp_du_index ++;
				// if(std::isnan((*ff)[i])){
				// 	printf("new ff[%d] = %f, ff[u] = %f, total_comp_dist_from_u = %f, comp_dist_from_u[%d] = %f\n", i, (*ff)[i], (*ff)[u], total_comp_dist_from_u, i, comp_dist_from_u[comp_du_index]);
				// 	exit(1);
				// }
				#pragma endregion //CC
				

				old_comp_weight += i_weight;
			}
		}
		//AP本尊 u 的 ff 要加上除了自己以外，外部的所有ff。
		(*ff)[u] += total_comp_dist_from_u - comp_dist_from_u[0];
		// printf("ff[%d] = %f\n", u, (*ff)[u]);

		/**
		 * 到這裡為止
		 * old_comp_weight = 原本AP所在的component的total weight - 原本的AP的weight
		 * 
		 * @brief
		 * 要把 AP本尊 (u)的 weight 改成，代表在AP本尊所在的component之外有多少點是AP本尊代表的
		*/
		// ALSO SET THE WEIGHT OF U and other idvs IN IDENTICAL_SETS THING, IF THEY ARE SO
		if ((*idv_track)[u] != -1) {
			int idx_of_u = (*idv_track)[u];
			double tobe_added = (old_comp_weight - u_weight) / ((*identical_sets_c)[idx_of_u] - 1);
			(*identical_sets)[idx_of_u][0].weight += (old_comp_weight - u_weight);
			for (int i = 1; i < (*identical_sets_c)[idx_of_u]; i++) {
				(*identical_sets)[idx_of_u][i].weight += tobe_added;
				// reflect weights
				weight[(*identical_sets)[idx_of_u][i].id] += tobe_added;
			}
		}
		else {
			//讓 AP本尊 (u) 紀錄他代表外部的多少點
			weight[u] += old_comp_weight - u_weight;
		}


		/**
		 * @brief 把新創建的 AP分身 的 weight 也計算好
		 * 作法 :
		 * 1. 取得每個AP分身所在的component的weight總和 i_weight
		 * 2. i_weight = i_weight - weight[i]，代表 i_weight 扣掉 AP分身 的 weight
		 * 3. 判斷 AP分身 是否是 identical vertex，然後再用不同做法去 assign weight
		 * @todo [Wait]
		 * 在看完 identical vertex之後，要回來看這個
		*/
		for (int i = startnextvid; i < endnextvid; i++) {
			if (art_track[i - initnVtx] == u) {
				double i_weight = totalw_idv(i, (*nVtx), (*padj), (*pxadj), weight, (*idv_track), (*identical_sets),
											tmark, tbfsorder);
				if ((*idv_track)[i] == -1)
					i_weight -= weight[i];
				else
					i_weight -= (*identical_sets)[(*idv_track)[i]][0].weight;

				if ((*idv_track)[i] == -1) {
					weight[i] += old_comp_weight - i_weight;
				}
				else {
					int idx_of_i = (*idv_track)[i];
					double tobe_added = (old_comp_weight - i_weight) / ((*identical_sets_c)[idx_of_i] - 1);
					(*identical_sets)[idx_of_i][0].weight += (old_comp_weight - i_weight);
					for (int j = 1; j < (*identical_sets_c)[idx_of_i]; j++) {
						(*identical_sets)[idx_of_i][j].weight += tobe_added;
						// reflect..
						weight[(*identical_sets)[idx_of_i][j].id] += tobe_added;
					}
				}
			}
		}


		/**
		 * component 	: 紀錄這個component 原本的所有nodeID之外，也在這裡紀錄 新加入的AP分身的ID
		 * reversecomp 	: 記錄這個component的nodeID對應到  component_arr 中的 index 
		*/
		for (int i = startnextvid; i < endnextvid; i++) {
			component[*len] = i;
			reversecomp[i] = *len;
#ifdef DEBUG
			int assoc_artpoint = art_track[i - initnVtx];
			if (assoc_artpoint == -1)
				printf("sth is wrong!!\n");
#endif
			(*len)++;
			mark[i] = 1;
			total_weights_of_each_comp[comp_ids_of_each_v[i]] = total_weight_at_the_beginning; //這裡不知道是否有記錯，或是他就是要這樣記
		}

		free(counter_extxadj);
		newoldneigs = (vertex *) myre_alloc (newoldneigs, 2 * (*nVtx) * sizeof(vertex));
		componentid = (int *) myre_alloc (componentid, (*nVtx) * sizeof(int));
	}

	*bc = (Betweenness *) myre_alloc (*bc, sizeof(Betweenness) * (*nVtx)); //myre_allocations for bc are done

	for (int i = veryoldnVtx; i < (*nVtx); i++)
		(*bc)[i] = 0.;
		

	if ((*nVtx) > (*size4)) {
		*size4 = *nVtx;
		bfsorder = (vertex *) myre_alloc (bfsorder, (*nVtx) * sizeof(vertex));
		endpred = (vertex *) myre_alloc (endpred, (*nVtx) * sizeof(vertex));
		level = (vertex *) myre_alloc (level, (*nVtx) * sizeof(vertex));
		sigma = (pathnumber *) myre_alloc (sigma, (*nVtx) * sizeof(pathnumber));
		delta = (Betweenness *) myre_alloc (delta, (*nVtx) * sizeof(Betweenness));
		ordered_comp = (vertex *) myre_alloc (ordered_comp, sizeof(vertex) * (*nVtx));
		reverse_ordered_comp = (vertex *) myre_alloc (reverse_ordered_comp, sizeof(vertex) * (*nVtx));
		ordered_weight = (double *)myre_alloc (ordered_weight, sizeof(double) * (*nVtx));
		ordered_cardinality = (card_info*)myre_alloc (ordered_cardinality, sizeof(card_info) * (*nVtx));
		all_graphs_xadj = (int *)myre_alloc (all_graphs_xadj, sizeof(int) * (*nVtx));
		all_graphs_len = (int *) myre_alloc (all_graphs_len, sizeof(int) * (*nVtx));
		labels = (int *) myre_alloc (labels, sizeof(int) * (*nVtx));
		nd = (int *) myre_alloc (nd, sizeof(int) * (*nVtx));
		l = (int *) myre_alloc (l, sizeof(int) * (*nVtx));
		h = (int *) myre_alloc (h, sizeof(int) * (*nVtx));
		art_points = (vertex *) myre_alloc (art_points, sizeof(vertex) * (*nVtx));
	}

	if (*nVtx > (*size2)) { //(*size2) was (2 * (initnVtx + 1)) at the beginning
		*size2 = *nVtx;
	}

	// allocated memories are freed
	free(componentid);
	free(newoldneigs);
	Zoltan_Bucket_Free(bs);
	
#ifdef BCCOMP_DBG
	// Graph Check : to be removed
	for (int i = 0; i < *nVtx; i++) {
		for (myindex j =(*pxadj)[i]; j <(*pxadj)[i+1]; j++) {
			int v = (*padj)[j];
			if (v != -1) {
				int flag = 0;
				for (myindex k = (*pxadj)[v]; k < (*pxadj)[v+1]; k++) {
					if ((*padj)[k] == i) {
						flag = 1;
						break;
					}
				}
				if (flag == 0) {
					printf("\n%d and %d do not match each other\n",i+1,v+1);
					break;
				}
			}
		}
		printf("\n");
	}
#endif

	// for(int i = 0 ; i < *nVtx ; i ++){
	// 	// printf("ff[%d] = %f\033[K\r", i, ff[i]);
	// 	printf("ff[%d] = %f\n", i, ff[i]);
	// }
	// printf("\n");
	/**
	 * 這裡還看不懂，大致上是要利用 bucket 來加速尋找 clique
	*/
	// Bucket is repopulated since it seems easier and cheaper
	int max_degree = 0;
	for (int i = 0; i < *nVtx; i++) {
		if (((*pxadj)[i+1] - (*pxadj)[i]) > max_degree)
			max_degree =(*pxadj)[i+1] -(*pxadj)[i];
	}
	*bs = Zoltan_Bucket_Initialize(max_degree+1, *len);
	for (vertex i = 0; i < *len; i++) {
		int j = component[i];
		if (j != -1) {
			int degree_of_i = 0;
			for (myindex k =(*pxadj)[j]; k <(*pxadj)[j+1]; k++) {
				if ((*padj)[k] != -1)
					degree_of_i++;
			}
			Zoltan_Bucket_Insert(bs, i, degree_of_i);
		}
	}
	return;
}


/**
 * 使用 bridge 切完之後，原本的component可能變成多個component，
 * 但是 len (Component nodeNum)不變
*/
void articulation_point_detection_multiple_cc (int len, int nVtx, vertex* art_points, int* artc, vertex* stack, int* dfn,
		int* l, int* parent, int* already_art, vertex* markcomp, vertex* component,
		vertex* adj, int* xadj) {

	memset(markcomp, 0, sizeof(vertex) * nVtx);
	for (int i = 0; i < len; i++) {
		int u = component[i];
		if ((u != -1) && (markcomp[u] == 0)) {
			memset(dfn, 0, sizeof(int) * nVtx);
			memset(l, 0, sizeof(int) * nVtx);
			memset(parent, 0, sizeof(int) * nVtx);
			articulation_point_detection (u, art_points, artc, stack, dfn, l, parent, already_art, markcomp,
					adj, xadj);
		}
	}
}



/**
 * Note : 
 * 1. vertex* art_points 	= AP 都存在這個 array 裡面
 * 2. int* already_art 		= 照 nodeID 紀錄誰是 AP
*/
void articulation_point_detection (vertex u, vertex* art_points, int* artc,
		vertex* stack, int* dfn, int* l, int* parent, int* already_art, vertex* markcomp,
		vertex* adj, int* xadj) {

	int num = 1;

	int root = u;
	markcomp[u] = 1;
	int stackc = 0;
	stack[stackc] = u;
	dfn[u] = num++;
	l[u] = dfn[u];

	/**
	 * 透過 dfn 跟 l 來判斷誰是 AP :
	 * 不過這種做法有可能會有多個連續的 AP
	*/
	while (stackc >= 0) {
		vertex v = stack[stackc];
		bool flag = false;
		for (myindex j = xadj[v]; j < xadj[v+1]; j++) {
			vertex w = adj[j];
			if ((w != -1) && (markcomp[w] == 0)) {
				markcomp[w] = 1;
				stack[++stackc] = w;
				dfn[w] = num++;
				parent[w] = v;
				l[w] = dfn[w];
				flag = true;
				break;
			}
		}
		if (!flag) {
			stackc--;
			for (myindex j = xadj[v]; j < xadj[v+1]; j++) {
				vertex w = adj[j];
				if (w != -1) {
					if ((w != parent[v]) && (dfn[w] < dfn[v]))
						l[v] = (l[v]<dfn[w])?l[v]:dfn[w];
					else if (v == parent[w]) {
						l[v] = (l[v]<l[w])?l[v]:l[w];
						if ((l[w] >= dfn[v]) && (v != root) && (already_art[v] == 0)) {
							art_points[(*artc)++] = v;
							already_art[v] = 1;
						}
					}
				}
			}
		}
	}

	//檢查 DFS 的 source 是否也是 AP
	int chno = 0;
	for (myindex j = xadj[u]; j < xadj[u+1]; j++) {
		vertex w = adj[j];
		if ((w != -1) && (parent[w] == u))
			chno++;
	}
	if ((chno > 1) && (already_art[u] == 0)) {
		art_points[(*artc)++] = u;
		already_art[u] = 1;
	}
	return;
}
