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
 * @brief [Hanging]
 * @todo 
 * 1. bridge_removal的 if, else if, else if 的 statement 都還沒看
 * 
 * 2. 這邊用來 BFS traverse 的 for迴圈 可以減少
*/

#include "bc-seq-brandes.h"


void bridge_removal (int nVtx, int len, vertex* bridges, int bridges_c, int* numof_removed_edges,
		Bucket* bs, vertex* component, vertex* reversecomp, vertex* adj, int* xadj, int* tadj,
		double* weight, Betweenness*bc, int* idv_track, int* identical_sets_c, idv_info** identical_sets,
		int* tmark, int* tbfsorder, bool dd, double* total_weights_of_each_comp, int* comp_ids_of_each_v,
		int* comp_no, double* CCs, double* ff) {

	double total_w_of_comp;
	int count = 0;

	//i += 2 是因為他存bridge的時候，是用bridges[i] = u, bridges[i + 1] = v, 表示(u, v)這條edge是bridge
	for (int i = 0; i < bridges_c; i+=2) {
		int u = bridges[i];
		int v = bridges[i+1];

		int bridgestillthere = 0;
		for (myindex j = xadj[u]; j < xadj[u+1]; j++) {
			if (v == adj[j]) {
				bridgestillthere = 1;
				break;
			}
		}

		if (bridgestillthere == 1) {
			total_w_of_comp = total_weights_of_each_comp[comp_ids_of_each_v[u]];

			for (myindex j = xadj[u]; j < xadj[u+1]; j++) {
				//這裡就把(u, v)這條 bridge 從 u 與 v 的 adj(csrE) 斷開
				if (v == adj[j]) {
					adj[j] = -1;
					(*numof_removed_edges)++;
					for (myindex k = xadj[v]; k < xadj[v+1]; k++) {
						if (adj[k] == u) {
							adj[k] = -1;
							break;
						}
					}
					break;
				}
			}


			/**
			 * 計算整個 component 有多少nodes，然後把 u 能探到的 nodes 的 comp_id 改成新的comp_no，代表這是最新的component
			 * 
			 * total_weight : 計算整個 component 有多少 nodes
			 * 
			*/
			double wu = 0;
			double wv = 0;
			double total_weight = 0;
			memset (tmark, 0, sizeof(int) * nVtx);
			int endofbfsorder = 1;
			tbfsorder[0] = u;
			int cur = 0;
			tmark[u] = 1;
			while (cur != endofbfsorder) {
				vertex v = tbfsorder[cur];
				comp_ids_of_each_v[v] = (*comp_no);
				if (idv_track[v] == -1)
					total_weight += weight[v];
				else {
					int idx_of_v = idv_track[v];
					total_weight += identical_sets[idx_of_v][0].weight;
				}
				for (myindex j = xadj[v]; j < xadj[v + 1]; j++) {
					vertex w = adj[j];
					if ((w != -1) && (tmark[w] == 0)) {
						tbfsorder[endofbfsorder++] = w;
						tmark[w] = 1;
					}
				}
				cur++;
			}
			/**
			 * wu : 不是u所在的component的所有node數量，而是 u 能夠探訪到的所有nodes數量，因為某些 node 已經開始有 weight，所以上面的 traverse 可以記數到除了 v 的 component之外的 所有 nodes 數量
			 * wv : v所在的component的所有node數量
			 * 這下面看不懂就畫圖
			*/
			wu = total_weight;
			wv = total_w_of_comp - wu; //check過了，這個 total_w_of_comp - wu，算出來，跟 v 自己去 traverse 自己的 component 是一樣的


			(*comp_no)++;

			total_weights_of_each_comp[comp_ids_of_each_v[u]] =	total_w_of_comp;  //這是對的，因為每次都要從整張圖一開始的component去看
			total_weights_of_each_comp[comp_ids_of_each_v[v]] = total_w_of_comp;  //這是對的，因為每次都要從整張圖一開始的component去看

			if(std::isnan(wu) || std::isnan(wu) || std::isnan(total_w_of_comp)){
				printf("[ERROR] wu = %f, wv = %f, total_weight_of_each_comp = %f\n", wu, wv, total_w_of_comp);
				exit(1);
			}
			else{
				// printf("cid_u = %d, wu = %f, cid_v = %d, wv = %f, total_weight_of_each_comp = %f\n", comp_ids_of_each_v[u], wu, comp_ids_of_each_v[v], wv, total_w_of_comp);
			}

			int idx_of_u = idv_track[u];
			bool is_u_idv = (idx_of_u == -1) ? false : true;
			int idx_of_v = idv_track[v];
			bool is_v_idv = (idx_of_v == -1) ? false : true;


			#pragma region CC
			/**
			 * @brief 取得以下兩個值
			 * comp_dist_from_u : u 所在的 component 的所有 nodes 對 u 的距離總和
			 * comp_dist_from_v : v 所在的 component 的所有 nodes 對 v 的距離總和
			 * 
			 * CC的 bridge removal 會比 BC久，因為 BC 在 bridge removal 只要算最短路徑個數(大約是兩個component的 nodes數相乘)，
			 * 而 CC 的 bridge removal 除了要算 weight(代表的nodes個數)，還要算 ff(代表的距離)
			*/

			double comp_dist_from_u = 0;
			double comp_dist_from_v = 0;

			double* dist_arr = (double*)malloc(sizeof(double) * nVtx);

			//先取得 comp_dist_from_u
			memset(tmark, 0, sizeof(int) * nVtx);
			memset(dist_arr, 0, sizeof(double) * nVtx);
			endofbfsorder = 1;
			tbfsorder[0] = u;
			cur = 0;
			tmark[u] = 1;
			dist_arr[u] = 0;
			while (cur != endofbfsorder) {
				vertex v = tbfsorder[cur];
				
				for (myindex j = xadj[v]; j < xadj[v + 1]; j++) {
					vertex w = adj[j];
					if ((w != -1) && (tmark[w] == 0)) {
						tbfsorder[endofbfsorder++] = w;
						tmark[w] = 1;

						dist_arr[w] = dist_arr[v] + 1;
						if(idv_track[w] == -1){
							comp_dist_from_u += dist_arr[w] * weight[w] + ff[w];
						}
						else{
							int idx_of_w = idv_track[w];
							comp_dist_from_u += dist_arr[w] * identical_sets[idx_of_w][0].weight + ff[w];
						}
					}
				}
				cur++;
			}

			//再取得 comp_dist_from_v
			memset(tmark, 0, sizeof(int) * nVtx);
			memset(dist_arr, 0, sizeof(double) * nVtx);
			endofbfsorder = 1;
			tbfsorder[0] = v;
			cur = 0;
			tmark[v] = 1;
			dist_arr[v] = 0;
			while (cur != endofbfsorder) {
				vertex currentNodeID = tbfsorder[cur];
				
				for (myindex j = xadj[currentNodeID]; j < xadj[currentNodeID + 1]; j++) {
					vertex w = adj[j];
					if ((w != -1) && (tmark[w] == 0)) {
						tbfsorder[endofbfsorder++] = w;
						tmark[w] = 1;

						dist_arr[w] = dist_arr[currentNodeID] + 1;
						if(idv_track[w] == -1){
							comp_dist_from_v += dist_arr[w] * weight[w] + ff[w];
						}
						else{
							int idx_of_w = idv_track[w];
							comp_dist_from_v += dist_arr[w] * identical_sets[idx_of_w][0].weight + ff[w];
						}
					}
				}
				cur++;
			}

			// printf("u = %d, comp_dist_from_u = %f\n", u, comp_dist_from_u);
			// printf("v = %d, comp_dist_from_v = %f\n", v, comp_dist_from_v);
			free(dist_arr);
			#pragma endregion //CC




			/**
			 * @todo [Wait]
			 * if, else if, else if 的 statement都還沒看是怎麼處理的，等idv的處理看完之後再回來看這裡
			*/
			if (is_u_idv && is_v_idv) {

				// bc updates for u's idv set
				bc[identical_sets[idx_of_u][1].id] += (wv * (wu - identical_sets[idx_of_u][0].weight)) / (identical_sets_c[idx_of_u] - 1);
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", (wv * (wu - identical_sets[idx_of_u][0].weight)) / (identical_sets_c[idx_of_u] - 1), identical_sets[idx_of_u][1].id+1);
#endif

				bc[identical_sets[idx_of_u][1].id] += (identical_sets[idx_of_u][1].weight - 1) * wv;
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", (identical_sets[idx_of_u][1].weight - 1) * wv, identical_sets[idx_of_u][1].id+1);
#endif

				for (int i = 2; i < identical_sets_c[idx_of_u]; i++) {
#ifdef BCCOMP_DBG
					printf("%lf IS ADDED TO BC[%d]\n", (identical_sets[idx_of_u][i].weight - identical_sets[idx_of_u][1].weight) * wv, identical_sets[idx_of_u][i].id+1);
#endif
					bc[identical_sets[idx_of_u][i].id] += (identical_sets[idx_of_u][i].weight - identical_sets[idx_of_u][1].weight) * wv;
				}

				// bc updates for v's idv set
				bc[identical_sets[idx_of_v][1].id] += (wu * (wv - identical_sets[idx_of_v][0].weight)) / (identical_sets_c[idx_of_v] - 1);
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", (wu * (wv - identical_sets[idx_of_v][0].weight)) / (identical_sets_c[idx_of_v] - 1), identical_sets[idx_of_v][1].id+1);
#endif

				bc[identical_sets[idx_of_v][1].id] += (identical_sets[idx_of_v][1].weight - 1) * wu;
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", (identical_sets[idx_of_v][1].weight - 1) * wu, identical_sets[idx_of_v][1].id+1);
#endif

				for (int i = 2; i < identical_sets_c[idx_of_v]; i++) {
#ifdef BCCOMP_DBG
					printf("%lf IS ADDED TO BC[%d]\n", (identical_sets[idx_of_v][i].weight - identical_sets[idx_of_v][1].weight) * wu, identical_sets[idx_of_v][i].id+1);
#endif
					bc[identical_sets[idx_of_v][i].id] += (identical_sets[idx_of_v][i].weight - identical_sets[idx_of_v][1].weight) * wu;
				}


				#pragma region CC

				double temp_u_ff = ff[u];
				double temp_v_ff = ff[v];

				ff[u] += temp_v_ff + comp_dist_from_v + wv;
				ff[v] += temp_u_ff + comp_dist_from_u + wu;

				#pragma endregion //CC


				double tobe_added_to_u = wv;
				double tobe_added_to_v = wu;

				// weight updates for u's idv set
				for (int i = 1; i < identical_sets_c[idx_of_u]; i++) {
					double tobe_added = tobe_added_to_u / (identical_sets_c[idx_of_u] - 1);
					identical_sets[idx_of_u][i].weight += tobe_added;
					//reflect..
					weight[identical_sets[idx_of_u][i].id] += tobe_added;
				}
				//                 // reflect it to real weight array (repr weight)
				//                 weight[identical_sets[idx_of_u][1].id] = identical_sets[idx_of_u][1].weight;
				// total idv weight is set
				identical_sets[idx_of_u][0].weight += tobe_added_to_u;

				// weight updates for v's idv set
				for (int i = 1; i < identical_sets_c[idx_of_v]; i++) {
					double tobe_added = tobe_added_to_v / (identical_sets_c[idx_of_v] - 1);
					identical_sets[idx_of_v][i].weight += tobe_added;
					// reflect..
					weight[identical_sets[idx_of_v][i].id] += tobe_added;
				}
				//                 // reflect it to real weight array (repr weight)
				//                 weight[identical_sets[idx_of_v][1].id] = identical_sets[idx_of_v][1].weight;
				// total idv weight is set
				identical_sets[idx_of_v][0].weight += tobe_added_to_v;
			}
			else if (is_u_idv && !is_v_idv) {

				bc[v] += wu * (wv - 1);
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", wu * (wv - 1), v+1);
#endif

				bc[identical_sets[idx_of_u][1].id] += (wv * (wu - identical_sets[idx_of_u][0].weight)) / (identical_sets_c[idx_of_u] - 1);
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", (wv * (wu - identical_sets[idx_of_u][0].weight))  / (identical_sets_c[idx_of_u] - 1), identical_sets[idx_of_u][1].id+1);
#endif

				bc[identical_sets[idx_of_u][1].id] += (identical_sets[idx_of_u][1].weight - 1) * wv;
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", (identical_sets[idx_of_u][1].weight - 1) * wv, identical_sets[idx_of_u][1].id+1);
#endif

				for (int i = 2; i < identical_sets_c[idx_of_u]; i++) {
					bc[identical_sets[idx_of_u][i].id] += (identical_sets[idx_of_u][i].weight - identical_sets[idx_of_u][1].weight) * wv;
#ifdef BCCOMP_DBG
					printf("%lf IS ADDED TO BC[%d]\n", (identical_sets[idx_of_u][i].weight - identical_sets[idx_of_u][1].weight) * wv, identical_sets[idx_of_u][i].id+1);
#endif
				}



				#pragma region CC

				double temp_u_ff = ff[u];
				double temp_v_ff = ff[v];

				ff[u] += temp_v_ff + comp_dist_from_v + wv;
				ff[v] += temp_u_ff + comp_dist_from_u + wu;

				#pragma endregion //CC
				


				weight[v] += wu;
				for (int i = 1; i < identical_sets_c[idx_of_u]; i++) {
					double tobe_added = wv / (identical_sets_c[idx_of_u] - 1);
					identical_sets[idx_of_u][i].weight += tobe_added;
					//reflect..
					weight[identical_sets[idx_of_u][i].id] += tobe_added;
				}
				//                 //reflect the weight change
				//                 weight[identical_sets[idx_of_u][1].id] = identical_sets[idx_of_u][1].weight;
				// set total weights of idv set
				identical_sets[idx_of_u][0].weight += wv;
			}
			else if (!is_u_idv && is_v_idv) {

				bc[u] += wv * (wu - 1);
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", wv * (wu - 1), u+1);
#endif

				bc[identical_sets[idx_of_v][1].id] += (identical_sets[idx_of_v][1].weight - 1) * wu;
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", (identical_sets[idx_of_v][1].weight - 1) * wu, identical_sets[idx_of_v][1].id+1);
#endif

				bc[identical_sets[idx_of_v][1].id] += (wu * (wv - identical_sets[idx_of_v][0].weight)) / (identical_sets_c[idx_of_v] - 1);
#ifdef BCCOMP_DBG
				printf("%lf IS ADDED TO BC[%d]\n", (wu * (wv - identical_sets[idx_of_v][0].weight)) / (identical_sets_c[idx_of_v] - 1),
						identical_sets[idx_of_v][1].id+1);
#endif



				for (int i = 2; i < identical_sets_c[idx_of_v]; i++) {
					bc[identical_sets[idx_of_v][i].id] += (identical_sets[idx_of_v][i].weight - identical_sets[idx_of_v][1].weight) * wu;
#ifdef BCCOMP_DBG
					printf("%lf IS ADDED TO BC[%d]\n", (identical_sets[idx_of_v][i].weight - identical_sets[idx_of_v][1].weight) * wu, identical_sets[idx_of_v][i].id+1);
#endif
				}

				#pragma region CC

				double temp_u_ff = ff[u];
				double temp_v_ff = ff[v];
				
				ff[u] += temp_v_ff + comp_dist_from_v + wv;
				ff[v] += temp_u_ff + comp_dist_from_u + wu;

				#pragma endregion //CC



				weight[u] += wv;
				for (int i = 1; i < identical_sets_c[idx_of_v]; i++) {
					double tobe_added = wu / (identical_sets_c[idx_of_v] - 1);
					identical_sets[idx_of_v][i].weight += tobe_added;
					// reflect..
					weight[identical_sets[idx_of_v][i].id] += tobe_added;
				}
				//                 //reflect the weight change
				//                 weight[identical_sets[idx_of_v][1].id] = identical_sets[idx_of_v][1].weight;
				// set total weights of idv set
				identical_sets[idx_of_v][0].weight += wu;


			}
			else {

				bc [u] += (wu - 1) * wv;
#ifdef BCCOMP_DBG
				printf("source br adds bc[%d]: %lf\n",u+1,(wu - 1) * wv);
#endif
				bc [v] += (wv - 1) * wu;
#ifdef BCCOMP_DBG
				printf("source br adds bc[%d]: %lf\n",v+1,(wv - 1) * wu);
#endif



				#pragma region CC

				double temp_u_ff = ff[u];
				double temp_v_ff = ff[v];

				ff[u] += temp_v_ff + comp_dist_from_v + wv;
				ff[v] += temp_u_ff + comp_dist_from_u + wu; 
				// printf("[Bridge : usual] ff[u] = %f, ff[v] = %f\r", ff[u], ff[v]);

				#pragma endregion //CC
				

				/**
				 * weight[u] : 除了u自己原先代表的nodes數量，現在也要代表 v 所在的component的 nodes 數量
				 * weight[v] : 除了v自己原先代表的nodes數量，現在也要代表 v 所在的component的 nodes 數量
				*/
				weight[u] += wv;
				weight[v] += wu;
			}

			int rev_v = reversecomp[v];
			Zoltan_Bucket_DecVal(bs, rev_v);
			int rev_u = reversecomp[u];
			Zoltan_Bucket_DecVal(bs, rev_u);

			//degree-1 removal
			if (!dd)
				remove_degree_1s (nVtx, component, reversecomp, bs, weight, xadj, adj, tadj, idv_track,
						identical_sets_c, identical_sets, bc, numof_removed_edges, tmark, tbfsorder,
						total_weights_of_each_comp, comp_ids_of_each_v, CCs, ff);

		}
	}
	

	//類似檢查的機制
	for (int i = 0; i < len; i++) {
		int u = component[i];
		int flag = 0;
		if (u != -1) {
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

	return;
}


void bridge_detection_multiple_cc (int len, int nVtx, int* labels, int* nd, int*l, int* h, vertex* component, int* xadj, vertex* adj, vertex* bridges, int* bridges_c) {
	/**
	 * 取得在 u 在他的component中的nodes數量，如果 v 與 u 同個component則不會再次成為source，
	 * 重複上述直到所有node都已 visited
	*/
	int* mark = new int[nVtx];
	memset (mark, 0, sizeof(int) * nVtx);
	memset (labels, -1, sizeof(int) * nVtx);
	for (int i = 0; i < len; i++) {
		int u = component[i];
		if ((u != -1) && (mark[u] == 0)) {
			int* bfsorder = new int[nVtx];
			int endofbfsorder = 1;
			bfsorder[0] = u;
			int cur = 0;
			mark[u] = 1;
			while (cur != endofbfsorder) {
				vertex v = bfsorder[cur];
				for (myindex j = xadj[v]; j < xadj[v + 1]; j++) {
					vertex w = adj[j];
					if ((w != -1) && (mark[w] == 0)) {
						bfsorder[endofbfsorder++] = w;
						mark[w] = 1;
					}
				}
				cur++;
			}
			delete[] bfsorder;
			bridge_detection (u, endofbfsorder, nVtx, labels, nd, l, h, xadj, adj, bridges, bridges_c);
		}
	}
	delete[] mark;
	return;
}


void bridge_detection (vertex u, int remlen, int nVtx, int* labels, int* nd, int*l, int* h, int* xadj,
		vertex* adj, vertex* bridges, int* bridges_c) {

	vertex* stack = new vertex[remlen];
	vertex* newtree = new vertex[remlen];
	vertex* children = new vertex[xadj[nVtx]];
	vertex* nontree = new vertex[xadj[nVtx]];
	int* nxadj = new int[remlen+1]; // nontree edges
	int* cxadj = new int[remlen+1]; // children edges
	vertex* parents = new vertex[nVtx];
	vertex* markcomp = new vertex[nVtx];
	memset (parents, -1, sizeof(vertex) * nVtx);
	memset (nxadj, 0, sizeof(int) * (remlen + 1));
	memset (cxadj, 0, sizeof(int) * (remlen + 1));

	//labc : label count
	int labc = 0, c = remlen - 1;
	memset(markcomp, 0, sizeof(vertex) * nVtx);
	markcomp[u] = 1;
	int stackc = 0;
	stack[stackc] = u;

	/**
	 * traverse component, if one node visit a node that has label already, then stop traverse.
	 * @todo
	 * 1. 還看不懂labels[]在記啥，目前認為是記錄此次探訪中，探訪到同個component中的幾個點。
	*/
	while (stackc >= 0) {
		vertex v = stack[stackc];

#ifdef DEBUG
		if(c == -1) {

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
			exit(1);
		}
#endif

		newtree[c] = v;
		c--;
		stackc--;
		if (labels[v] == -1) {
			labels[v] = labc++;
		}
		else
			break;
		
		/**
		 * 1. xadj = csrV
		 * 2. adj = csrE
		*/
		for (myindex j = xadj[v]; j < xadj[v+1]; j++) {
			vertex w = adj[j];
			if (w != -1) {
				if (markcomp[w] == 0) {
					markcomp[w] = 1;
					stack[++stackc] = w;
				}
			}
		}
	}

	/**
	 * 建立BFS traverse tree：
	 * 1. 知道哪些 edge 	是 		BFS traverse tree 的構成edge (child edge 	: cxadj)
	 * 2. 知道哪些 edge 	不是 	BFS traverse tree 的構成edge (nontree edge 	: nxadj)
	 * 
	 * Variables :
	 * 1. cxadj = child edge 	的 csrV,	children 	= child edge 	的 csrE : (cxadj, children)由前往後是root到leaf
	 * 2. nxadj = nontree edge	的 csrV,	nontree 	= nontree edge 	的 csrE : (nxadj, nontree)由前往後是root到leaf
	*/
	memset(markcomp, 0, sizeof(vertex) * nVtx);
	for (int i = 0; i < remlen; i++) {
		cxadj[i + 1] = cxadj[i];
		nxadj[i + 1] = nxadj[i];
		int v = newtree[remlen - i - 1];
		markcomp[v] = 1;
		for (myindex j = xadj[v]; j < xadj[v+1]; j++) {
			vertex w = adj[j];
			if (w != -1) {
				if (markcomp[w] == 0) {
					parents[w] = v;
					children[cxadj[i + 1]++] = w;
					markcomp[w] = 1;
				}
				else if (w != parents[v]) {
					nontree[nxadj[i + 1]++] = w;
				}
			}
		}
	}

	/**
	 * @brief 紀錄每個node在一個component中連接多少child =>
	 * 從tree的底部開始往上trace (cxadj由前往後是 root 到 leaf)
	 * i = 0 : tree 的 底
	 * i = remlen - 1 是 root
	*/
	for (int i = 0; i < remlen; i++) {
		int u = newtree[i];
		nd[u] = 1;
		for (int j = cxadj[remlen - 1 - i]; j < cxadj[remlen - i]; j++) {
			nd[u] += nd[children[j]];
		}
	}

	/**
	 * 從 newtree 底部的 node 開始看 child neighbors 跟 non-tree neighbors 
	 * l => 記錄自己可以摸到最靠近 root 的 node 的 label
	 * h => 紀錄自己可以摸到最遠離 root 的 node 的 label
	*/
	for (int i = 0; i < remlen; i++) {
		int u = newtree[i];
		int minl = labels[u];
		int maxh = labels[u];
		for (int j = cxadj[remlen - 1 -i]; j < cxadj[remlen - i]; j++) {
			int lm = l[children[j]];
			int hm = h[children[j]];
			if (lm < minl)
				minl = lm;
			if (hm > maxh)
				maxh = hm;
		}
		for (int j = nxadj[remlen - 1 -i]; j < nxadj[remlen - i]; j++) {
			int lm = labels[nontree[j]];
			if (lm < minl)
				minl = lm;
			if (lm > maxh)
				maxh = lm;
		}
		l[u] = minl;
		h[u] = maxh;
	}


	for (int i = 0; i < remlen; i++) {
		int v = newtree[i];
		for (int j = cxadj[remlen - 1 - i]; j < cxadj[remlen - i]; j++) {
			int w = children[j];
			if ((l[w] == labels[w]) && (h[w] < labels[w] + nd[w])) {
				bridges[(*bridges_c)++] = v;
				bridges[(*bridges_c)++] = w;
#ifdef DEBUG
printf("%d-%d is a bridge!!\n",v+1,w+1);
#endif
			}
		}
	}
	delete[] stack;
	delete[] newtree;
	delete[] children;
	delete[] nontree;
	delete[] nxadj;
	delete[] cxadj;
	delete[] parents;
	delete[] markcomp;
	return;
}
