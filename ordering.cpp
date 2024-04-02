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
 * 
*/
void reorder_graph (int nVtx, int len, vertex* component, vertex* ordered_comp, vertex* reverse_ordered_comp,
			double* weight, double* ordered_weight, int** pxadj, int** padj, vertex* newxadj, vertex* newadj,
			int* all_graphs_xadj, int* agx, int* num_comp, int* k, int* c, int* biggest_cc_after) {

	int newlen = 0;
	vertex* markcomp = new vertex[nVtx];
	memset(markcomp, 0, sizeof(vertex) * nVtx);
	for (vertex i = 0; i < len; i++) {
		//從 0 ~ len 有可能會因為 AP切割，bridge切割，而造成在 0 ~ len之間 有不同 component 的 node
		int u = component[i];
		if ((u != -1) && (markcomp[u] == 0)) {
			(*num_comp)++;
			markcomp[u] = 1;
			int comp_counter = 0;
			ordered_comp[(*k) + comp_counter] = u;
			ordered_weight[(*k) + comp_counter] = weight[u];
			reverse_ordered_comp[u] = comp_counter;
			comp_counter++;
			int cur = comp_counter - 1;

			/**
			 * @brief
			 * [Note]
			 * Order 的意思是 依照同個 component 的 node 都放在一起，放在order_comp的 array裡面
			 * 以 u 為 source 進行BFS 
			 * 蒐集 u 所在的 component 的所有 node
			*/
			while (cur != comp_counter) {
				vertex v = ordered_comp[(*k) + cur];
				for (myindex j =(*pxadj)[v]; j <(*pxadj)[v+1]; j++) {
					vertex w =(*padj)[j];
					if (w != -1) {
						if (markcomp[w] == 0) {
							markcomp[w] = 1;
							ordered_comp[(*k) + comp_counter] = w;
							ordered_weight[(*k) + comp_counter] = weight[w];
							reverse_ordered_comp[w] = comp_counter;
							comp_counter++;
						}
					}
				}
				cur++;
			}

			if (*biggest_cc_after < comp_counter)
				*biggest_cc_after = comp_counter;

			/**
			 * 建立
			 * 1. newadj
			 * 2. newxadj
			*/
			int old_k = (*k);
			newlen = comp_counter;
			for (vertex i = 0; i < newlen; i++) {
				newxadj[(*k) + 1] = newxadj[(*k)];
				int u = ordered_comp[old_k + i];
				for (myindex j =(*pxadj)[u]; j <(*pxadj)[u + 1]; j++) {
					int e =(*padj)[j];
					if (e != -1) {
//						assert ((*c) < size3);
						newadj[(*c)++] = reverse_ordered_comp[e];
						newxadj[(*k) + 1]++;
					}
				}
				(*k)++;
			}

			//這個應該是每個component的第一個點的csrV offset
			all_graphs_xadj[(*agx)++] = (*k);
		}
	}

	delete[] markcomp;
	return;
}


/**
 * @brief [done]
 * 這個function 會把 padj(csrE)裡面的 -1 都移除掉，並且重新計算 pxadj
*/
void remove_minus_ones_in_graph (int nVtx, int* num_edges, int** pxadj, int** padj) {

	int nEdge = (*pxadj)[nVtx];
	int global_count = 0;
	int previous_index = 0;
	int i;
	for (i = 0; i < nVtx; i++) {
		for (int j = (*pxadj)[i]; j < (*pxadj)[i+1]; j++) {
			if ((*padj)[j] != -1) {
				(*padj)[global_count++] = (*padj)[j];
			}
		}

		(*pxadj)[i] = previous_index;
		previous_index = global_count;

	}
	(*pxadj)[nVtx] = previous_index;
	*num_edges = previous_index;

	

	return;
}
