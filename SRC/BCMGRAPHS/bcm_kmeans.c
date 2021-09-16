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

double distance(
				double * XI, //data of the i-th node of dim coordinates
				double * C,  //cluster centroid (1 x dim)
				int dim, //dimensionality
				int N, //number of nodes/rows of X
				int NC //number of clusters
				) 
	{ 	
		int id;
	  	double dist=0.0; 
	  	for (id=0;id<dim;id++)  //for each dimension
	  	  dist+=pow((XI[id*N]-C[id*NC]),2);  
	  	dist=sqrt(dist); //distance of the node from the cluster center 
		return dist;
	}

void getclosestCentroids(
	double *X, //data
	double *centroids, //vector of the clusters centroids (NC x dim)
	int N,   //number of nodes/rows of X
	int dim,  //dimensionality
	int NC, //number of clusters
	//OUTPUT
	int * indices //indices is the vector of the labels of each node/row, its dimension is N (number of nodes/raws)
						)
	{
	
		int i,idx,j;
		double min_dist ,dist ;
		for (i=0;i<N;i++)  //  for each node/row of X
			{
			 	idx=0; //assign the node to the first  cluster  
				min_dist =distance(X+i,centroids,dim,N,NC); //set the minimum distance as the distance from the first cluster 
				for (j=1;j<NC;j++) //check the other clusters 
					{
						dist =distance(X+i,centroids+j,dim,N,NC); //compute the distance of the i-th node from the j-th cluster centroid 
						 if(dist < min_dist)  //if it's closer 
						 	{
						 		min_dist = dist; //update distance 
          						idx = j; //update label--> assign the i-th node to the j-th cluster 
						 	}
					}
				indices[i] = idx; //store the labels 
				
			}
			
	}


void computeCentroids
	(
		double *X,   //data
		int N, //number of nodes/rows of X
		int dim, // dimensionality
		int * indices,  //indices is the vector of the labels of each node/row of X
		int NC,  //number of clusters  
		//OUTPUT
		double * centroids	//centroids 		
	)
	{ 
		int i,ic,j,cc,icc;

		double * xi;
		
		for (ic=0; ic<NC;ic++) //for each cluster 
		{   
			cc=0; //initialize cluster count to zero 
			for (i=0;i<N;i++) //for each node 
				if (indices[i]==ic) //if the i-th node belongs to the ic-th cluster 
					cc++; //count cluster members
			xi  = (double *) calloc(cc*dim, sizeof(double));
			icc=0;
			for (i=0;i<N;i++) //for each node 
				if (indices[i]==ic) //if the i-th node belongs to the ic-th cluster 
					{
						for (j=0;j<dim;j++) //save in x_i the coordinates of the node 	
							xi[icc+j*cc]=X[i+j*N];  //which is a row of X
						icc++;
					}
			//now xi is a matrix containing all the nodes (rows of X) of the ic-th cluster
			//its dimension is cc x dim , ordered by columns 
			for (j=0;j<dim;j++) //for each coordinate of the node
				{
					centroids[ic+j*NC]=0.0; //initialize
					for (i=0;i<cc;i++)  //for all the nodes  of the ic-th cluster 
						centroids[ic+j*NC]+=xi[i+j*cc]; 
					centroids[ic+j*NC]=centroids[ic+j*NC]/cc; 
					//compute the centroid coordinates as the mean value in each column 
				}
			free(xi);	
		}
		
	}

void kmeans(//input parameters
	double *X,  //data
	int NC, //number of clusters
	int max_iterations, //max iterations for stopping criteria of the algorithm
	int dim, //dimensionality ( nhrc)
	int N, //number of original graph nodes/rows of X
	//OUTPUT parameters
	int *indices, //vector of the labels ( output)
	double *centroids, //cluster centroids( output)
	int *iter // number of iterations ( output)
				) 
{   
	int ii,jj,stop,i;
   	int *indicesold; 
    indicesold = (int *) calloc(N, sizeof(int));
	stop=1; //initialize stop so that the cycle doesn't stop
	for (ii=0;ii<N;ii++)  //initialize indices 
			indices[ii]=1;		     
    i=0;  
    while(stop && i <  max_iterations) // reiterate until either the indices don't change anymore
	                                  // or you reach the max number of iterations 
  		{   i++; //increment i 
		    getclosestCentroids(X,centroids,N,dim,NC,indices);//assign each node to the closest centroid 
  		    computeCentroids(X,N,dim,indices, NC,centroids); //compute the new cluster centroids 
  		    ii=0;
			while ( ii<N && indices[ii]==indicesold[ii]) //if clusters haven't changed
			     ii++;     // see the next difference 
			if (ii==N)  stop=0;  //stop if no component of indices has changed 
			for (ii=0;ii<N;ii++) 
				indicesold[ii]=indices[ii];  //prepare for the next iteration
  		  	
		}
    *iter=i;		
    free(indicesold);
}
