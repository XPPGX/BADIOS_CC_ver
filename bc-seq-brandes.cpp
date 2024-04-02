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
#include "graph.hpp"


int main(int argc, char *argv[]) {

	int i, nVtx, *tadj, *xadj, *adj, nTry = 1;
	myindex nEdge;
	char *filename = argv[1];
	char timestr[20];
	char timeprestr[20];
	char time1str[20];
	char time2str[20];
	char deg1remstr[20];
	char bridgedetstr[20];
	char bridgeremstr[20];
	char cliquedetstr[20];
	char cliqueremstr[20];
	char artpdetstr[20];
	char artpremstr[20];
	char idvdetstr[20];
	char idvremstr[20];
	char bfsorderstr[20];

	int numof_removed_edges = 0;
	int numof_newly_created_vertices = 0;
	int numof_art_points = 0;
	int numof_identical_vertices = 0;
	int biggest_cc_before = 0;
	int biggest_cc_after = 0;
	int num_comp = 0;

	if (argc < 2) {
		fprintf(stderr, "Usage: %s\n"
			    "	<input filename>			// path of the input file (Chaco or Matrix Market format)\n"
			    "	<nTry>					// number of runs\n"
			    "	<0 for BADIOS, 1 for base-brandes>	// when BADIOS is selected, ordering is enabled as default\n"
			    "	<side-vertex_disable>               	// 1 to disable, 0 to enable side vertex removal\n"
			    "	<articulation-vertex_disable>       	// 1 to disable, 0 to enable articulation vertex copy\n"
			    "	<bridge_disable>                    	// 1 to disable, 0 to enable bridge removal\n"
			    "	<degree1-vertex_disable>            	// 1 to disable, 0 to enable degree-1 vertex removal\n"
			    "	<identical-vertex_disable>          	// 1 to disable, 0 to enable identical vertex removal (both type-I and type-II)\n"
				, argv[0]);
		exit(1);
	}

	if (argc >= 3) {
		nTry = atoi(argv[2]);
		if (!nTry) {
			fprintf(stderr, "nTry was 0, setting it to 1\n");
			nTry = 1;
		}
	}

	ReadGraph<int,int>(argv[1], &nVtx, &xadj, &adj, &tadj, NULL, NULL);
	nEdge = xadj[nVtx];

	Betweenness* bc = (Betweenness *) malloc(sizeof(Betweenness) * nVtx);
	for (i=0; i<nVtx; ++i)
		bc[i] = 0.;

	if (strrchr(argv[1], '/'))
		filename = strrchr(argv[1], '/') + 1;

	util::timestamp totaltime(0,0);
	util::timestamp preproctime(0,0);
	util::timestamp phase1time(0,0);
	util::timestamp phase2time(0,0);
	util::timestamp deg1remtime(0,0);
	util::timestamp bridgedettime(0,0);
	util::timestamp bridgeremtime(0,0);
	util::timestamp cliquedettime(0,0);
	util::timestamp cliqueremtime(0,0);
	util::timestamp artpdettime(0,0);
	util::timestamp artpremtime(0,0);
	util::timestamp idvdettime(0,0);
	util::timestamp idvremtime(0,0);
	util::timestamp bfsordertime(0,0);

	std::string algo_out;

	int maxvertex=nVtx;
	{
		char* mv;
		if ((mv = getenv("MAXVERTEX") ) != NULL)
		{
			maxvertex = atoi (mv); //TODO: boo unsafe!!
		}
	}

	int alg = atoi(argv[3]);

	bool disable_clique = (atoi(argv[4])==1)?true:false;
	bool disable_art = (atoi(argv[5])==1)?true:false;
	bool disable_bridge = (atoi(argv[6])==1)?true:false;
	bool disable_degree1 = (atoi(argv[7])==1)?true:false;
	bool disable_idv = (atoi(argv[8])==1)?true:false;

	double partial_bc_factor = -1;
#ifdef PARTIAL_BC
	partial_bc_factor = atof(argv[9]);
#endif

	if (alg == 0) {
	  main_bc(nVtx, &xadj, &adj, &tadj, &bc, nTry, totaltime, preproctime, phase1time, phase2time, deg1remtime, bridgedettime,
		bridgeremtime, cliquedettime, cliqueremtime, artpdettime, artpremtime, idvdettime, idvremtime, bfsordertime,
		&numof_removed_edges, &numof_art_points, &numof_newly_created_vertices, &numof_identical_vertices, &biggest_cc_before,
		&biggest_cc_after, &num_comp, disable_clique, disable_art, disable_bridge, disable_degree1, disable_idv
#ifdef PARTIAL_BC
		, partial_bc_factor
#endif
		);
	}
	else
		base_bc(nVtx, xadj, adj, bc, maxvertex, nTry, totaltime, phase1time, phase2time
#ifdef PARTIAL_BC
		, partial_bc_factor
#endif
		);

	int maxv = 0;
	for (int i=1; i< nVtx; i++) {
		if (bc[i] > bc[maxv])
			maxv=i;
	}

	totaltime /= nTry;
	totaltime.to_c_str(timestr, 20);
	preproctime /= nTry;
	preproctime.to_c_str(timeprestr, 20);
	phase1time /= nTry;
	phase1time.to_c_str(time1str, 20);
	phase2time /= nTry;
	phase2time.to_c_str(time2str, 20);

	deg1remtime /= nTry;
	deg1remtime.to_c_str(deg1remstr, 20);

	bridgedettime /= nTry;
	bridgedettime.to_c_str(bridgedetstr, 20);
	bridgeremtime /= nTry;
	bridgeremtime.to_c_str(bridgeremstr, 20);

	cliquedettime /= nTry;
	cliquedettime.to_c_str(cliquedetstr, 20);
	cliqueremtime /= nTry;
	cliqueremtime.to_c_str(cliqueremstr, 20);

	artpdettime /= nTry;
	artpdettime.to_c_str(artpdetstr, 20);
	artpremtime /= nTry;
	artpremtime.to_c_str(artpremstr, 20);

	idvdettime /= nTry;
	idvdettime.to_c_str(idvdetstr, 20);
	idvremtime /= nTry;
	idvremtime.to_c_str(idvremstr, 20);

	bfsordertime /= nTry;
	bfsordertime.to_c_str(bfsorderstr, 20);

	std::cout.precision(15);

	std::cout<<"filename:                      "<<filename<<endl
			<<" nVtx:                          "<<nVtx<<endl
		    <<" nEdge:                         "<<nEdge/2<<endl
		    <<" TotTime:                       "<<timestr<<endl
		    <<" PreprocTime:                   "<<timeprestr<<endl
		    <<" Ph1Time:                       "<<time1str<<endl
		    <<" Ph2Time:                       "<<time2str<<endl
		    <<" Deg1RemovalTime:               "<<deg1remstr<<endl
		    <<" BridgeDetectionTime:           "<<bridgedetstr<<endl
		    <<" BridgeRemovalTime:             "<<bridgeremstr<<endl
		    <<" CliqueDetectionTime:           "<<cliquedetstr<<endl
		    <<" CliqueRemovalTime:             "<<cliqueremstr<<endl
		    <<" ArtPDetectionTime:             "<<artpdetstr<<endl
		    <<" ArtPRemovalTime:               "<<artpremstr<<endl
		    <<" IdvDetectionTime:              "<<idvdetstr<<endl
		    <<" IdvRemovalTime:                "<<idvremstr<<endl
		    <<" BfsOrderTime:                  "<<bfsorderstr<<endl
		    <<" NumofRemovedEdges:             "<<numof_removed_edges<<endl
		    <<" NumofArtPoints:                "<<numof_art_points<<endl
		    <<" NumofNewlyCreatedVertices:     "<<numof_newly_created_vertices<<endl
		    <<" NumofIdenticalVerticesRemoved: "<<numof_identical_vertices<<endl
		    <<" sizeOfBiggestCcBefore:         "<<biggest_cc_before<<endl
		    <<" sizeOfBiggestCcAfter:          "<<biggest_cc_after<<endl
		    <<" numComp:                       "<<num_comp<<endl
		    <<" partialBCfactor:               "<<partial_bc_factor<<endl
		    <<" MaxBC:"<<bc[maxv]<<" at "<<maxv<<endl
		    <<algo_out<<std::endl;


	free(bc);
	free(adj);
	free(xadj);
	free(tadj);

	return 0;
}
