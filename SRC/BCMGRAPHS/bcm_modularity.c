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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include"bcm_graphs.h"

double modularity
	(
        bcm_CSRMatrix *A, //adjacency matrix of the graph
		int *K,    //vector of the labels
		int NC   //number of clusters
	)
  {
   double Q,W,W1;
   int i,j,ic,ik,jidi,jbp,nzrows_A,jj;
   Q=0.0;
   int num_rows = bcm_CSRMatrixNumRows( A );
   int num_cols = bcm_CSRMatrixNumCols( A );
   int nnz = bcm_CSRMatrixNumNonzeros( A );
   double *A_data=bcm_CSRMatrixData(A);
   int     *A_i = bcm_CSRMatrixI(A);
   int     *A_j = bcm_CSRMatrixJ(A);
   int *idi;
   idi = (int *) calloc(num_rows, sizeof(int));

   for (ic=0;ic<NC;ic++)   //for each cluster
   		{  W=0.0;W1=0.0;
		   jidi=0;
		   for (ik=0;ik<num_rows;ik++) //check in the label vector K
   		    {   
			    if (K[ik]==ic)    //which nodes belong to the ic-th cluster
   		        	{idi[jidi]=ik;  //and save their index in the vector idi
					jidi++; 
				    }
   		    }
   		    for(i=0; i< jidi; i++) // for each node of the ic-th cluster
   		    { 
			 jbp= A_i[idi[i]]-A_i[0]; //compute the number of non-zeros until the first element of the idi-th raw 
			 nzrows_A=A_i[idi[i]+1]-A_i[idi[i]]; // number of nonzeros in the idi-th raw
	   	     for(j=0; j<nzrows_A; j++)  //for each element of the idi-th raw
	         		{
					 W1=W1+A_data[jbp+j]; //sum to W1 the element A(idi,j)
					 for (jj=0;jj<jidi;jj++) //for each element of idi
			 		 	{if (A_j[jbp+j]==idi[jj])  //if the element belongs to the idi[jj]-th column                                 
				     	 W=W+A_data[jbp+j];  //sum to W the element A(idi,idi[jj])
   		             	}
					}
   		    }
   		    W=W/nnz; //average	
			W1=W1/nnz; //average
			Q=Q+(W-W1*W1); //modularity
   	    }
   free(idi);
   return Q; //OUTPUT = modularity 
  }
