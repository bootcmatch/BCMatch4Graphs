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
#include "bcm_util.h"
#include"bcm_graphs.h"


double entropy
	(
		mge_Metadata* metadata, //matrix of the attributes
		int *K,    //vector of the labels
		int NC   //number of clusters
	)
  {
    double e;	
	double pac,coeff;			
	double *ent; //internal sum for the entropy 
	ent=(double *) calloc(NC, sizeof(double)); 
	int ic,il,i,j;
	int nl;
	nl= metadata -> num_properties;//number of labels
	int N;
	N=metadata -> num_verts;  //number of supernodes (original graph nodes)
	int nc;
	int  nva;
	double *v; //internal sum for the entropy 
	v=(double *) calloc(N, sizeof(double)); 
	int * r;
	
	e=0.0;
    for (ic=0;ic<NC;ic++) //for each cluster
	{
		for (il=0;il<nl;il++) //for each label
		{
			nc=0;  
			for (j=0; j<N;j++) //for each supernode		    
				if (K[j]==ic) //that belongs to the cluster 
					{	
						v[nc]=(metadata -> metadata_matrix)[j*nl+il]; //take the value of the il-th label on that node 
						nc++;  //count the elements of the cluster	
					}
			if (nc>0)
			{
			order( v,0,nc-1); //and order the values of the attribute on the cluster elements 	
			r=(int *) calloc(N, sizeof(int)); //recurrences of each value of the il-th attribute
			nva=0;
			for (i=1; i<nc;i++) //over all the cluster elements  
				{
					if (v[i]==v[i-1])
						r[nva]++; //count the recurrences of each value of the il-th attribute
						else 
							nva++; //next value 
				}
			for (i=0;i<=nva;i++)	 //for each value that label can assume
				{
				pac=r[i]+1;//number of recurrences of the i-th value that the ikl-th label can assume
				// it's one more because the C numbering strats from 0
				pac/=nc; // in percentage of the elements of the ic-th cluster 
				pac=pac*log(pac)/log(2); //<0 because it's a log of a precentage (always<1)
				ent[ic]+=pac;
				
				}
			free(r);
			}
			
		}
		e-=ent[ic]*nc/N; //change sign because entropy is always >0
	}
   free(ent);
   free(v);
   return e; //OUTPUT = entropy 
  }
