/*
                BCMatch4Graphs
     Bootstrap AMG based on Compatible weighted Matching for Graphs, version 1.0
    (C) Copyright 2021
                       Pasqua D'Ambra         IAC-CNR, IT
                       Panayot S. Vassilevski Portland State University, OR USA

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions, and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
    3. The name of the BCMatch4Graphs group or the names of its contributors may
       not be used to endorse or promote products derived from this
       software without specific written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BCMATCH4GRAPHS GROUP OR ITS CONTRIBUTORS
  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/
/* Contributor Clara De Santis 22/03/2021 */

#include "bcm_metadata_exp.h"
#include "bcm_util.h"

bcm_CSRMatrix* bcm_MetadataAugmentation
	(
		bcm_CSRMatrix* A, //adjacency matrix of the original graph (N x dim)
 		mge_Metadata* metadata//matrix of the labels (N x nl)
 	)
{    
	
	int nl = mge_MetadataNumProperties(metadata); //number of labels 
	int il,i,j,k;
    int *A_i,*A_j,N,num_cols,nnz  ;
	A_i    = bcm_CSRMatrixI(A);
  	A_j    = bcm_CSRMatrixJ(A);
  	N    = bcm_CSRMatrixNumRows(A);//num_rows/number of supernodes 
  	num_cols    = bcm_CSRMatrixNumCols(A);
  	nnz         = bcm_CSRMatrixNumNonzeros(A);
    double* A_data = bcm_CSRMatrixData(A); //adjacency matrix
	double *l; //vector of the values that the label assunes on all the nodes (1 x N)
	l=(double *) calloc(N, sizeof(double));
	int *nva; //number of values that each label can assume 
	nva=(int*) calloc(nl, sizeof(int));
	int * r ;
	r=(int *) calloc(N, sizeof(int)); //recurrences of each value of the il-th attribute, the possible values  can be N at most   
	double *va; //vector of the values that the label assumes on all the nodes (1 x N)
	va=(double *) calloc(N*nl, sizeof(double));
	int num_new_nodes=0;
	for (il=0;il<nl;il++ )
    	{
		for (i=0;i<N;i++)
    		l[i]=metadata -> metadata_matrix[i*nl+il];
		//let's count the number of values that the il-th label can assume
		order(l,0,N-1);
		nva[il]=0; 
		va[il*N+nva[il]]=l[0]; 
		for (i=1; i<N;i++) //over all the values of l (first excluded )
			{
				if (l[i]==l[i-1])
					r[nva[il]]++; //count the recurrences of each value of the il-th attribute
				 else 
						{
							nva[il]++; //new value
							va[il*N+nva[il]]=l[i]; //save it 
					 	}
			}
		num_new_nodes +=nva[il]+1;
		}
	int new_edge=0; //new edges numbering 
    int new_vertex=N; //new vertex start numbering 
	int*matrix_i;
	matrix_i=calloc(N+num_new_nodes+1,sizeof(int));
	for (i=0;i<N;i++)
		{
			matrix_i[i+1]=A_i[i+1]-A_i[i];
			//printf("original edges of node %d=%d\n",i,matrix_i[i+1]);
		}
	double newedge_weight = 1.0; // weight of the new edge 
	for (il=0;il<nl;il++)//for each label 
    	for (j=0;j<=nva[il];j++) //for each possible value of the k-th attribute 
        {	
			for (i=0;i<N;i++) //for each supernode 
            	if((metadata -> metadata_matrix[i*nl+il])==va[il*N+j])
            		{
						matrix_i[i+1]++;
						matrix_i[new_vertex+1]++;
            			new_edge++;
            			
					}
			new_vertex++;			
		}
	printf("new nodes = %d, new edges= %d\n",num_new_nodes,new_edge );
  	for(j=0;j<new_vertex; j++) //for each row/node 
  		{
  			matrix_i[j+1]+=matrix_i[j];
		  }
    
	bcm_CSRMatrix* A_ext = bcm_CSRMatrixCreate(new_vertex, new_vertex, matrix_i[new_vertex]);
	bcm_CSRMatrixI(A_ext) = matrix_i;
    bcm_CSRMatrixInitialize(A_ext);
	int* matrix_j = bcm_CSRMatrixJ(A_ext); 
	double* matrix_data = bcm_CSRMatrixData(A_ext); 
    int original_edges,considered_edges ;
    int nval;
    int * cons_edges = calloc(num_new_nodes,sizeof(int));
	
    for (i=0;i<N;i++) //for each supernode 
    	{
			considered_edges=0;
			original_edges = A_i[i+1] - A_i[i]; //number of connections/edges  for that supernode in the original graph 
        	for (j= 0; j < original_edges;j++)  //for each edge of that node in the original graph 
        		{//I just copy that part of the adjacency matrix of the original graph 
					matrix_j[matrix_i[i] + j] = A_j[A_i[i] + j]; 
            		matrix_data[matrix_i[i] + j] = A_data[A_i[i] + j]; 
				}
    		new_vertex=N;
    		nval=0;
			for (il=0;il<nl;il++)//for each label 
				for (k=0;k<=nva[il];k++) //for each possible value of the k-th attribute 
      				{	
      				    if((metadata -> metadata_matrix[i*nl+il])==va[il*N+k])
							{
							matrix_j[matrix_i[new_vertex]+cons_edges[nval]]=i;
							matrix_data[matrix_i[new_vertex]+cons_edges[nval]] = newedge_weight; 
							cons_edges[nval]++;
								
								
							matrix_j[matrix_i[i] + original_edges + considered_edges ] = new_vertex; //connection of the supernode  with the label node  
			    			matrix_data[matrix_i[i] + original_edges + considered_edges ] = newedge_weight; //the weight is always 1 / create an edge 
            				considered_edges++; 
							}
						new_vertex++;
						nval++;
					}		
		}
free(va);
free(r);
free(nva);
free(l);
free(cons_edges);
return A_ext;
}
