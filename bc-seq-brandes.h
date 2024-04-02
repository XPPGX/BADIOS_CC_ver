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

#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>

#include "ulib.h"
#include "bucket.h"
#include "timestamp.hpp"


// #define DEBUG
// #define BCCOMP_DBG
// #define VLIST
// #define REDUCTION_STEPS
// #define CC_PROFILE
// #define PARTIAL_BC
// #define JUST_PREPROC
// #define DEBUG

#define THROW_AWAY 0
#define INIT_SIZE_FOR_IDV_SETS 200
#define ENHANCE_FACTOR 1.3

typedef int attr_id_t;
typedef double Betweenness;
typedef long long myindex;
typedef int vertex;
typedef double pathnumber;

struct idv_info {
	int id;
	double weight;
};

struct card_info {
	int cardinality;
	double total_weight;
};

using namespace std;

void extract_cc_info (int nVtx, int* xadj, int* adj, int* numof_CCs, int* maxCC_nvtx,
		int* maxCC_nedge, int* total_nvtx, int* total_nedge);
int select_kernel (int start, int end, double* ordered_weight, card_info* ordered_cardinality);
void* myre_alloc (void* ptr, size_t size);
void idv_detection_and_merge(int type, int len, int* component, int* reversecomp, int* xadj,
		int* adj, int nVtx, int* idv_sets_size, int* one_set_size,
		int* idv_track_size, int* next_idvset_id, int** idv_track, int** identical_sets_c, int** identical_sets_sz,
		idv_info*** identical_sets, double* weight, Betweenness* bc,
		Bucket* bs, int* numof_removed_edges, int* numof_identical_vertices, util::timestamp& idvdet,
		util::timestamp& idvrem);
bool equal_vertices(int type, vertex v1, vertex v2, int* xadj, int* adj);

void articulation_point_copy (int* nVtx, int artc, int* nextvid, int* len, int* size1, int* size2, int* size4,
		vertex* art_points, vertex* art_track, int** pxadj,
		int** padj, Betweenness** bc, int* num_edges, int initnVtx, vertex* bfsorder,
		int* endpred, int* level, pathnumber* sigma, Betweenness* delta,
		int* mark, double* weight, vertex* component, vertex* ordered_comp, vertex* reversecomp,
		vertex* reverse_ordered_comp, double* ordered_weight, card_info* ordered_cardinality,
		int* all_graphs_xadj, int* all_graphs_len,
		int* labels, int* nd, int* l, int* h, Bucket* bs, int* next_idvset_id,
		int* idv_track_size, int** idv_track, int** identical_sets_c, int** identical_sets_sz, idv_info*** identical_sets,
		int* idv_sets_size, int one_set_size, int* tmark, int* tbfsorder, double* total_weights_of_each_comp,
		int* comp_ids_of_each_v, int* comp_no, double* CCs, double* ff);

void graph_check (int* xadj, vertex* adj, int nVtx);

void assign_component_ids (int cid, vertex v, vertex u, int* componentid, vertex* adj, int* xadj, int nVtx, double* comp_dist_from_u);

void articulation_point_detection_multiple_cc (int len, int nVtx, vertex* art_points, int* artc, vertex* stack, int* dfn,
		int* l, int* parent, int* already_art, vertex* markcomp, vertex* component,
		vertex* adj, int* xadj);

void articulation_point_detection (vertex u, vertex* art_points, int* artc,
		vertex* stack, int* dfn, int* l, int* parent, int* already_art, vertex* markcomp,
		vertex* adj, int* xadj);

void remove_covs (int len, vertex* component, int cov_i, vertex* clique_only_v, vertex* reversecomp, Bucket* bs, int nVtx,
		double* weight, int* xadj, vertex* adj, vertex* bfsorder, int* endpred, int* level, pathnumber* sigma,
		vertex* Pred, Betweenness* delta, Betweenness* bc, int* numof_removed_edges, int* idv_track,
		int* identical_sets_c, idv_info** identical_sets, int next_idvset_id, double* totalw_of_covs,
		double* total_weights_of_each_comp, int* comp_ids_of_each_v);

int compare(const void *p1, const void *p2);
void bc_comp_for_cov (int source, int nvtx, double* weight, int* xadj,
		vertex* adj, vertex* bfsorder, int* endpred, int* level, pathnumber* sigma,
		vertex* Pred, Betweenness* delta, Betweenness* bc, int* idv_track,
		int* identical_sets_c, idv_info** identical_sets);
int is_part_of_clique (int u, int* xadj, vertex* adj, int clique_neig_num);
void find_clique_only_vertices (vertex* clique_only_v, int* cov_i, vertex* component, int len, int* xadj, vertex* adj);
int totalsize (vertex u, int nVtx, vertex* adj, int* xadj);
double totalw_idv(vertex u, int nVtx, vertex* adj, int* xadj, double* weight, int* idv_track, idv_info** identical_sets,
		int* tmark, int* tbfsorder);
double totalw(vertex u, int nVtx, vertex* adj, int* xadj, double* weight);

void remove_degree_1s (int nVtx, vertex* component, vertex* reversecomp, Bucket* bs, double* weight,
		int* xadj, vertex* adj, vertex *tadj, int* idv_track, int* identical_sets_c, idv_info** identical_sets,
		Betweenness* bc, int* numof_removed_edges, int* tmark, int *tbfsorder, double* total_weights_of_each_comp,
		int* comp_ids_of_each_v, double* CCs, double* ff);

void compute_bc_weight (int start, int end, vertex* ordered_comp, double* ordered_weight, vertex* newxadj,
		vertex* newadj, vertex* bfsorder, int* endpred, int* level, pathnumber* sigma,
		vertex* Pred, Betweenness* delta, Betweenness* bc, util::timestamp& phase1time, util::timestamp& phase2time);

void compute_bc_weight_card (int start, int end, vertex*  ordered_comp, double*  ordered_weight, vertex*  newxadj,
		vertex*  newadj, vertex*  bfsorder, int*  endpred, int*  level, pathnumber*  sigma,
		vertex*  Pred, Betweenness*  delta, Betweenness*  bc, card_info* ordered_cardinality,
		util::timestamp& phase1time, util::timestamp& phase2time);

void compute_bc_card (int start, int end, vertex*  ordered_comp, vertex*  newxadj, vertex*  newadj,
		vertex*  bfsorder, int*  endpred, int*  level, pathnumber*  sigma, vertex*  Pred, Betweenness*  delta,
		Betweenness*  bc, card_info* ordered_cardinality, util::timestamp& phase1time, util::timestamp& phase2time);

void compute_bc_base (int start, int end, vertex* ordered_comp, vertex* newxadj,
		vertex* newadj, vertex* bfsorder, int* endpred, int* level, pathnumber* sigma,
		vertex* Pred, Betweenness* delta, Betweenness* bc, util::timestamp& phase1time, util::timestamp& phase2time);

void bridge_detection (vertex u, int remlen, int nVtx, int* labels, int* nd, int*l, int* h, int* xadj,
		vertex* adj, vertex* bridges, int* bridges_c);


void bridge_detection_multiple_cc (int len, int nVtx, int* labels, int* nd, int*l, int* h, vertex* component, int* xadj,
		vertex* adj, vertex* bridges, int* bridges_c);

void bridge_removal (int nVtx, int len, vertex* bridges, int bridges_c, int* numof_removed_edges,
		Bucket* bs, vertex* component, vertex* reversecomp, vertex* adj, int* xadj, int* tadj,
		double* weight, Betweenness*bc, int* idv_track, int* identical_sets_c, idv_info** identical_sets,
		int* tmark, int* tbfsorder, bool dd, double* total_weights_of_each_comp, int* comp_ids_of_each_v, int* comp_no, double* CCs, double* ff);

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
);

int base_bc(int nVtx, int *xadj, int *adj, Betweenness *bc,
        int maxvertex, int nTry, //algo parameter
        util::timestamp& totaltime, util::timestamp& phase1time, util::timestamp& phase2time
#ifdef PARTIAL_BC
        , double partial_bc_factor
#endif
        );


void reorder_graph (int nVtx, int len, vertex* component, vertex* ordered_comp, vertex* reverse_ordered_comp,
			double* weight, double* ordered_weight, int** pxadj, int** padj, vertex* newxadj, vertex* newadj,
			int* all_graphs_xadj, int* agx, int* num_comp, int* k, int* c, int* biggest_cc_after);

void remove_minus_ones_in_graph (int nVtx, int* num_edges, int** pxadj, int** padj);
