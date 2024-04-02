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

#include <vector>
#include <map>
#include "ulib.h"
//#include "mmio.h"
#include <string>
#include <sstream>

#define MAXLINE 128*(1024*1024)


using namespace std;
/* reads the MeTiS format graph (Cannot handle multiple vertex weights!)
   file pointed by filename into pxadj, padjncy, padjncyw, ppvw:
   if you don't want edge weights pass NULL as padjncyw param 
   same for vertex weights pass NULL as ppvw 
   note that *pxadj, *padjncy, *padjncyw is allocated by by this function and must be freed by user of 
   this function*/


typedef struct pair {
	long long f;
	long long s;
} Pair;

int pcmp(const void *v1, const void *v2){
	long long diff = (((Pair *)v1)->f - ((Pair *)v2)->f);
	if (diff != 0)
		return diff;
	else
		return (((Pair *)v1)->s - ((Pair *)v2)->s);
}

template <typename VtxType, typename WeightType>
void ReadGraphFromFile(FILE *fpin, VtxType *numofvertex, VtxType **pxadj, VtxType **padjncy, WeightType **padjncyw,
					WeightType **ppvw)
{
  VtxType *xadj, *adjncy,  nvtxs, nedges, fmt, readew, readvw, edge, i, k, ncon;
  WeightType*adjncyw=NULL, *pvw=NULL;
  char *line;
    
    line = (char *)malloc(sizeof(char)*(MAXLINE+1));
    
    do {
        fgets(line, MAXLINE, fpin);
    } while (line[0] == '%' && !feof(fpin));
    
    if (feof(fpin)) 
        errexit("empty graph!!!");
    
    fmt = 0;
    {
      std::string s = line;
      std::stringstream ss (s);
      ss>>nvtxs>>nedges>>fmt>>ncon;
    }
    *numofvertex = nvtxs;
    
    readew = (fmt%10 > 0);
    readvw = ((fmt/10)%10 > 0);
    if (fmt >= 100) 
        errexit("invalid format");
    
    nedges *=2;
    
    xadj = *pxadj = imalloc(nvtxs+2, "ReadGraph: xadj");
    adjncy = *padjncy = imalloc(nedges, "ReadGraph: adjncy");
    if (padjncyw)
        adjncyw = *padjncyw = imalloc(nedges, "ReadGraph: adjncyw");
    if (ppvw)
        pvw = *ppvw = imalloc(nvtxs+1, "ReadGraph: adjncyw");
    
    for (xadj[0]=0, k=0, i=0; i<nvtxs; i++) {
        char *oldstr=line, *newstr=NULL;
        int  ewgt=1, vw=1;
        
        do {
            fgets(line, MAXLINE, fpin);
	} while (line[0] == '%' && !feof(fpin));
        
        if (strlen(line) >= MAXLINE-5) 
            errexit("\nBuffer for fgets not big enough!\n");
        
        if (readvw) {
            vw = (int)strtol(oldstr, &newstr, 10);
            oldstr = newstr;
	}
        
        if (ppvw)
            pvw[i] = vw;	
        
        for (;;) {
            edge = (int)strtol(oldstr, &newstr, 10) -1;
            oldstr = newstr;
            
            if (readew) {
                ewgt = (int)strtol(oldstr, &newstr, 10);
                oldstr = newstr;
	    }
            
            if (edge < 0) {
                break;
            }

            if (edge==i)
                errexit("Self loop in the graph for vertex %d\n", i);
            adjncy[k] = edge;
            if (padjncyw)
                adjncyw[k] = ewgt;
            k++;
	} 
        xadj[i+1] = k;
    }
    
    if (k != nedges) 
        errexit("k(%d)!=nedges(%d) and i:%d", k, nedges,i);
    
    free(line);

    return;
}


/* reads the Matrix Market format graph */

template <typename VtxType, typename WeightType>
void ReadGraphFromMMFile(FILE *matfp, VtxType *numofvertex, VtxType **pxadj, VtxType **padjncy, WeightType **padjncyw,
		WeightType **ppvw)
{

	Pair *coords, *new_coords;
	int m, n, itemp, jtemp;
	long long nnz, tnnz, i, j, onnz;
	int *xadj, *adj, *adjncyw, *pvw;

	int maxLineLength = 1000000;
	int value = 0;
	char line[maxLineLength];
	int num_items_read;

	/* set return null parameter values, in case we exit with errors */
	m = nnz = 0;

	/* now continue scanning until you reach the end-of-comments */
	do {
		if (fgets(line, 1000000, matfp) == NULL)
			value = 1;
	} while (line[0] == '%');

	/* line[] is either blank or has M,N, nz */
	if (sscanf(line, "%d %d %lld", &m, &n, &nnz) == 3) {
		value = 1;
	}
	else {
		do {
			num_items_read = fscanf(matfp, "%d %d %lld", &m, &n, &nnz);
			if (num_items_read == EOF)
				return;
		}
		while (num_items_read != 3);
		value = 0;
	}

//	printf("matrix banner is read %d - %d, %lld nnz\n", m, n, nnz);
	coords = (Pair*) malloc(sizeof(Pair) * 2 * nnz);

	tnnz = 0;
	for(i = 0; i < nnz; i++) {
		fscanf(matfp, "%d %d\n", &itemp, &jtemp);

		if(itemp != jtemp) {
			coords[tnnz].f = itemp;
			coords[tnnz++].s = jtemp;
			coords[tnnz].f = jtemp;
			coords[tnnz++].s = itemp;
		}
	}

	qsort(coords, tnnz, sizeof(Pair), pcmp);

	onnz = 1;
	for(i = 1; i < tnnz; i++) {
		if(coords[i].f != coords[onnz-1].f || coords[i].s != coords[onnz-1].s) {
			coords[onnz].f = coords[i].f;
			coords[onnz++].s = coords[i].s;
		}
	}

	*numofvertex = n;
	xadj = *pxadj = (int*) malloc((n+1) * sizeof(int));
	adj = *padjncy = (int*) malloc(onnz * sizeof(int));
    if (padjncyw)
        adjncyw = *padjncyw = imalloc (nnz, "ReadGraph: adjncyw");
    if (ppvw)
        pvw = *ppvw = imalloc (n+1, "ReadGraph: adjncyw");

    map <long, int> reed;
    map <long, int>::iterator reed_it;
    new_coords = (Pair*) malloc(sizeof(Pair) * 2 * nnz);
    long vno = 0;
    // map the ids
    for(i = 0; i < onnz; i++) {
    	long temp = coords[i].f;
    	reed_it = reed.find(temp);
    	if (reed_it == reed.end()) {
    	    reed.insert (make_pair (temp, vno));
    	    new_coords[i].f = vno++;
    	}
    	else
    		new_coords[i].f = reed_it->second;

    	temp = coords[i].s;
    	reed_it = reed.find(temp);
    	if (reed_it == reed.end()) {
    	    reed.insert (make_pair (temp, vno));
    	    new_coords[i].s = vno++;
    	}
    	else
    		new_coords[i].s = reed_it->second;

    }




    vector<vector<int> > entire_graph;
    entire_graph.resize(n);
    for(i = 0; i < onnz; i++) {
    	entire_graph[new_coords[i].f].push_back(new_coords[i].s);
    }


    xadj[0] = 0;
    j = 0;
    for(i = 1; i < n+1; i++) {
    	xadj[i] = xadj[i-1] + entire_graph[i-1].size();
    	for (unsigned int k = 0; k < entire_graph[i-1].size(); k++) {
    		adj[j++] = entire_graph[i-1][k];
    	}
    }


	free(coords);
	for(i = 0; i < m; i++)
		entire_graph[i].clear();

	entire_graph.clear();

	return;
}



template <typename VtxType, typename WeightType>
void ReadGraph(char *filename, VtxType *numofvertex, VtxType **pxadj, VtxType **padjncy, int** ptadj, WeightType **padjncyw,
		WeightType **ppvw)
{
    FILE *fpin = ufopen(filename, "r", "main: fpin");
    
    char * pch;
    pch = strstr (filename,".graph");

    if (pch != NULL)
    	ReadGraphFromFile (fpin, numofvertex, pxadj, padjncy, padjncyw, ppvw);
    else
    	ReadGraphFromMMFile (fpin, numofvertex, pxadj, padjncy, padjncyw, ppvw);


    int* adj = *padjncy;
    int* xadj = *pxadj;
    int n = *numofvertex;
    *ptadj = (int*) malloc(sizeof(int) * xadj[n]);
    int* tadj = *ptadj;

    int i;
    int* degs = (int*)malloc(sizeof(int) * n);
    int* myedges = (int*)malloc(sizeof(int) * xadj[n]);
    
    memcpy(degs, xadj, sizeof(VtxType) * n);    
    int ptr, j;
    for(i = 0; i < n; i++) {
      for(ptr = xadj[i]; ptr < xadj[i+1]; ptr++) {
	j = adj[ptr];
	myedges[degs[j]++] = i;
      }
    }    

    for(i = 0; i < n; i++) {
    	if(xadj[i+1] != degs[i]) {
    		printf("something is wrong\n");
    	}
    }

    memcpy(adj, myedges, sizeof(VtxType) * xadj[n]);       
    for(i = 0; i < n; i++) {
    	for(ptr = xadj[i]+1; ptr < xadj[i+1]; ptr++) {
    		if(adj[ptr] <= adj[ptr-1]) {
    			printf("is not sorted\n");
			}
		}
    }

    memcpy(degs, xadj, sizeof(VtxType) * n);
    for(i = 0; i < n; i++) {
    	for(ptr = xadj[i]; ptr < xadj[i+1]; ptr++) {
    		j = adj[ptr];
    		if (i < j) {
    			tadj[ptr] = degs[j];
    			tadj[degs[j]++] = ptr;
    		}
    	}
    }
    
    free(degs);
    free(myedges);
    
    for(i = 0; i < n; i++) {
    	for(ptr = xadj[i]; ptr < xadj[i+1]; ptr++) {
    		j = adj[ptr];
    		if((adj[tadj[ptr]] != i) || (tadj[ptr] < xadj[j]) || (tadj[ptr] >= xadj[j+1])) {
    			printf("error i %d j %d ptr %d\n", i, j, ptr);
				printf("error  xadj[j] %d  xadj[j+1] %d\n",  xadj[j], xadj[j+1]);
				printf("error tadj[ptr] %d\n", tadj[ptr]);
				printf("error adj[tadj[ptr]] %d\n", adj[tadj[ptr]]);
				exit(1);
    		}
    	}
    }
      
    ufclose(fpin);
    return;
}


