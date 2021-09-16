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

Contribution by Clara De Santis, July 2021 */

#include "bcm_metadata_exp.h"
#include "bcm.h"
#include "ioutil.h"
#include <float.h>
// #define DUMP_HIER
void dump_on_file(const char *prefix, bcm_BootAMG *boot_amg)
{
  int j,k;
  bcm_AMGHierarchy **Harray;
  Harray=bcm_BootAMGHarray(boot_amg);
  bcm_CSRMatrix    **A_array;
  bcm_CSRMatrix    **P_array;
  bcm_CSRMatrix    **L_array;
  bcm_CSRMatrix    **U_array;
  bcm_Vector       **D_array;

  for (k=0; k<1; k++) {
    A_array           = bcm_AMGHierarchyAArray(Harray[k]);
    P_array           = bcm_AMGHierarchyPArray(Harray[k]);
    L_array           = bcm_AMGHierarchyLArray(Harray[k]);
    U_array           = bcm_AMGHierarchyUArray(Harray[k]);
    D_array           = bcm_AMGHierarchyDArray(Harray[k]);
    int num_levels    = bcm_AMGHierarchyNumLevels(Harray[k]);
    char filename[81];
    for (j=0; j<num_levels; j++) {
      sprintf(filename,"%s-P-l%3.3d.mtx",prefix,j);
      if (P_array[j]!=NULL) bcm_CSRMatrixPrintMM(P_array[j],filename);
      sprintf(filename,"%s-AC-l%3.3d.mtx",prefix,j);
      bcm_CSRMatrixPrintMM(A_array[j],filename);
    }
  }      
}

typedef struct {
  int matrixformat;
  char *matrixfile;
  int solver_type;
  int max_hrc;
  double conv_ratio;
  int solver_it;
  int matchtype;
  int aggrsweeps;
  int aggrtype;
  int max_levels;
  int cycle_type;
  int coarse_solver;
  int relax_type;
  int prerelax_sweeps;
  int postrelax_sweeps;
  int NC;
  char *labelfile;
} parms_t;

void free_inparms(parms_t *inparms)
{
  if (inparms->matrixfile!=NULL) free(inparms->matrixfile);
  if (inparms->labelfile!=NULL) free(inparms->labelfile);
}

int get_inp_data(FILE *fp, parms_t *inparms)
{
  if (fp != NULL) {
    inparms->matrixformat      = int_get_fp(fp);
    inparms->matrixfile        = string_get_fp(fp);
    inparms->solver_type       = int_get_fp(fp);
    inparms->max_hrc           = int_get_fp(fp);
    inparms->conv_ratio        = double_get_fp(fp);
    inparms->solver_it         = int_get_fp(fp);
    inparms->matchtype         = int_get_fp(fp);
    inparms->aggrsweeps        = int_get_fp(fp);
    inparms->aggrtype          = int_get_fp(fp);
    inparms->max_levels        = int_get_fp(fp);
    inparms->cycle_type        = int_get_fp(fp);
    inparms->coarse_solver     = int_get_fp(fp);
    inparms->relax_type        = int_get_fp(fp);
    inparms->prerelax_sweeps   = int_get_fp(fp);
    inparms->postrelax_sweeps  = int_get_fp(fp);    
    inparms->NC                = int_get_fp(fp);    
    inparms->labelfile         = string_get_fp(fp);
  } else {
    return(-1);
  }
  return(0);
}

   /*This program applies the clustering of attributed graphs described in:
 
   P. D'Ambra, C. De Santis, P.S. Vassilevski, and L. Cutillo.
   Network clustering by embedding of attribute-augmented graphs, 2021.

   After the graph augmentation, this program builds a Bootstrap AMG with a 
   desired convergence rate for the Laplacian of a general graph and uses the estimated
   smooth vectors obtained at each bootstrap step as an estimate of its eigenvalues.
   Before computing and applying the Bootstrap AMG, the program reads the adjacency matrix of
   the graph and then compute the Laplacian matrix, finally
   apply a rank-1 update to have a s.p.d. matrix.
   The program also applies the modified k-means algorithm described in the above paper to obtain 
   final clustering of the original graph with the largest modularity.*/


int  main(int argc, char *argv[])
{
   bcm_BootAMGBuildData *bootamg_data;
   bcm_AMGApplyData *amg_cycle;
   bcm_CSRMatrix *A,*A_original;
   bcm_Vector *w;
   int i, *num_grid_sweeps, j, k;
  
   double ratio, normnew, normold;
   parms_t inparms;
   FILE    *fp=NULL;

   fprintf(stdout,
	   "Welcome to BCMATCH4GRAPHS version %s\n This is the clustering_AttributedGraph program\n",
	   BCM_VERSION_STRING);

   if (argc > 1) {
     fp = fopen(argv[1], "r");
   } else {
     fp = stdin;
   }

   if (get_inp_data(fp,&inparms) != 0) {
     fprintf(stderr,"Error getting input parms \n");
     exit(1);
   }
   /* read sparse matrix from file */
   switch(inparms.matrixformat)
   {
     case 0:
   {
     fprintf(stderr,"Reading in COO\n");
     A_original=bcm_COO2CSRMatrixRead(inparms.matrixfile);
   }
   break;
     case 1:
   {
     fprintf(stderr,"Reading in MM\n");
     A_original=bcm_MM2CSRMatrixRead(inparms.matrixfile);
   }
   break;
     case 2:
   {
     fprintf(stderr,"Reading in CSR\n");
     A_original=bcm_CSRMatrixRead(inparms.matrixfile);
   }
   break;
   }
   int num_rows_original = bcm_CSRMatrixNumRows( A_original );
   
   fclose(fp);
    //AUGMENTATION OF THE ORIGINAL GRAPH A_original with metadata 
    fp  = fopen(inparms.labelfile, "r");

    int  num_verts, num_labels, ik, label, ierr;
    int rowcount = 0;          //rows 
    int first_colcount = 0;   //number of columns in first row
    int colcount = 0;         //subsequent number of columns	
    //char delimiter='\t';      //delimiter for entries in the matrix
    char delimiter=' ';      //delimiter for entries in the matrix
    char c; // To store a character read from file 
      
	//Get first row
	for (c = getc(fp); rowcount != 1; c = getc(fp)){
		if (c == delimiter) first_colcount = first_colcount + 1;
		if (c == '\n') // New row once we reach newline
			rowcount = rowcount + 1;
	}

	//Get Other rows and check that they are compatible with the first row
	for (c = getc(fp); c != EOF; c = getc(fp)){
		if (c == delimiter) colcount = colcount + 1;
		if (c == '\n') { // New row once we reach newline
			rowcount = rowcount + 1;
			if(colcount!=first_colcount){
				printf("Error in row %i\n",rowcount);
				printf("Uneven number of columns %i %i",colcount,first_colcount);
				return 0;
			}
			colcount = 0;
		}
	}
	

   //Add one to the colcount since the last column is never counted
    num_labels = first_colcount + 1; 
    num_verts = rowcount;
    fclose(fp);
    mge_Metadata* metadata = mge_MetadataCreate(num_verts, num_labels);
    
    printf("Labels Read: numverts=%d numlabels=%d\n", num_verts, num_labels);
    fp  = fopen(inparms.labelfile, "r");
    i=0;
    for (j = 0; j < num_verts; j++)
    {
    	for (k = 0; k < num_labels; k++)
    	{
      		fscanf(fp, "%d", &label);
     	 	(metadata -> metadata_matrix)[i] = label; 
			    i++;
    	}
    }
    fclose(fp);
    
    //AUGMENTED GRAPH --> A
   A= bcm_MetadataAugmentation(A_original,metadata); 
   int num_rows = bcm_CSRMatrixNumRows( A );
   int num_cols = bcm_CSRMatrixNumCols( A );
   int nnz = bcm_CSRMatrixNumNonzeros( A );
   printf("numrows: %d\n",num_rows);
   printf("numcols: %d\n",num_cols);
   printf("nnz: %d\n",nnz);

   /* initialize data structure for AMG building. See the setup
      routines for changing default values */
   bootamg_data = bcm_BootAMGBuildDataInitialize();
   bcm_AMGBuildData *amg_data;
   amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);
   amg_cycle= bcm_AMGCycleInitialize();

/* building edge-vertex incidence matrix */
   bcm_CSRMatrix *B, *E, *ET, *LA, *LASPD;
   int L=1;
   B=bcm_CSRMatrixTriU(A,L);

   int  nrows_B  =  bcm_CSRMatrixNumRows(B);
   int  nnz_B  =  bcm_CSRMatrixNumNonzeros(B);
   int  *B_i = bcm_CSRMatrixI(B);
   int  *B_j = bcm_CSRMatrixJ(B);
   printf("numrowsB: %d\n",nrows_B);
   printf("numcolsB: %d\n",nnz_B);

   int nonzerosE=2*nnz_B;
   int *i_edge_dof, *j_edge_dof;
   double *E_data;

   i_edge_dof = (int *) calloc(nnz_B+1, sizeof(int));
   j_edge_dof= (int *) calloc(nonzerosE, sizeof(int));
   E_data= (double *) calloc(nonzerosE, sizeof(double));

   int l, nnz_rowB, kk, rowedge, irowedge;
   kk=0;
   rowedge=0;
   irowedge=0;
   for(ik=0; ik<nrows_B; ++ik)
   {
     nnz_rowB=B_i[ik+1]-B_i[ik];
     for(l=0; l< nnz_rowB; ++l)
     {
       i_edge_dof[irowedge]=2*irowedge;
       irowedge=++irowedge;

       j_edge_dof[kk]=ik;
       E_data[kk]=1;
       kk=kk+1;
       j_edge_dof[kk]=B_j[rowedge+l];
       E_data[kk]=-1;
       kk=kk+1;
     }
     rowedge=rowedge+nnz_rowB;
   }
   int num_edges=rowedge;
   int num_dofs=nrows_B;
   i_edge_dof[num_edges]=nonzerosE;
   printf("numrowsE: %d\n",num_edges);
   printf("numcolsE: %d\n",num_dofs);
   printf("numnnzE: %d\n",nonzerosE);

   E= bcm_CSRMatrixCreate(num_edges, num_dofs, nonzerosE);
   bcm_CSRMatrixI(E)=i_edge_dof;
   bcm_CSRMatrixJ(E)=j_edge_dof;
   bcm_CSRMatrixData(E)=E_data;

   /* transposing edge-vertex incidence matrix */

   printf("Transposing E \n");
   int *i_dof_edge, *j_dof_edge;
   ierr = bcm_CSRMatrixTranspose(E,&ET,1);
   i_dof_edge=bcm_CSRMatrixI(ET);
   j_dof_edge=bcm_CSRMatrixJ(ET);
   int it = bcm_CSRMatrixNumRows(ET);
   int jt = bcm_CSRMatrixNumCols(ET);

   printf("ET numrows %d\n", it);
   printf("ET numcols %d\n", jt);

   printf("computing laplacian of A \n");
   LA=bcm_CSRMatrixMultiply(ET,E);
   int nrow_LA = bcm_CSRMatrixNumRows(LA);

   bcm_Vector *DIAGLA;
   DIAGLA = bcm_CSRMatrixDiag(LA);
   double *diagla_val;
   diagla_val=bcm_VectorData(DIAGLA);

   double maxdeg;
   maxdeg=diagla_val[0];
   for(l=1; l<nrow_LA; l++)
   {
      if(diagla_val[l] > maxdeg) maxdeg=diagla_val[l];
   }
   printf("maxdeg=%e\n",maxdeg);

   double mindeg;
   mindeg=diagla_val[0];
   for(l=1; l<nrow_LA; l++)
   {
      if(diagla_val[l] < mindeg) mindeg=diagla_val[l];
   }
   printf("mindeg=%e\n",mindeg);

   double avgdeg;
   avgdeg=diagla_val[0];
   for(l=1; l<nrow_LA; l++)
   {
      avgdeg=avgdeg+diagla_val[l];
   }
    avgdeg=avgdeg/nrow_LA;
   printf("avgdeg=%e\n",avgdeg);

  /* applying rank-1 update for dealing with semipositiveness */

   printf("computing rank-1 update of laplacian of A \n");
   int *i_LA, *j_LA;

   i_LA=bcm_CSRMatrixI(LA);
   j_LA=bcm_CSRMatrixJ(LA);
   int nnz_row1 = i_LA[1]-i_LA[0];
   //int nnz_row1 = i_LA[2]-i_LA[1];
   if(nnz_row1 > 1) 
   {
     int jidx=j_LA[1];
     //int jidx=j_LA[2];
     bcm_CSRMatrix *spe, *spet;
     int *i_spe, *j_spe;
     double *data_spe;
     i_spe= (int *) calloc(2, sizeof(int));
     i_spe[0]=0;
     i_spe[1]=2;
     j_spe= (int *) calloc(2, sizeof(int));
     j_spe[0]=0;
     j_spe[1]=jidx;
     data_spe= (double *) calloc(2, sizeof(double));
     data_spe[0]=1.0;
     data_spe[1]=1.0;
     spe=bcm_CSRMatrixCreate(1,nrow_LA,2);
     bcm_CSRMatrixI(spe)=i_spe;
     bcm_CSRMatrixJ(spe)=j_spe;
     bcm_CSRMatrixData(spe)=data_spe;

     ierr = bcm_CSRMatrixTranspose(spe,&spet,1);

     bcm_CSRMatrix *SPE;
     SPE=bcm_CSRMatrixMultiply(spet,spe);
     LASPD=bcm_CSRMatrixAdd(LA,SPE);

   bcm_CSRMatrixDestroy(SPE);
   bcm_CSRMatrixDestroy(spe);
   bcm_CSRMatrixDestroy(spet);
   }
   else printf("Error: laplacian of A has only 1 nonzero in the first row.\n");

   bcm_AMGBuildDataCSRMatrix(amg_data)=LASPD;

   /* reading AMG algorithmic parameters */

   /* composition type  and maximum number of
      components for the bootstrap AMG */
   /* Single Hierarchy Type: matchingtype, number of aggregation sweeps
      for aggressive coarsening and aggregation type (smoothed or unsmoothed) */

   bcm_BootAMGBuildSetSolverType(bootamg_data, inparms.solver_type);
   bcm_BootAMGBuildSetMaxHrc(bootamg_data, inparms.max_hrc);
   bcm_BootAMGBuildSetDesiredRatio(bootamg_data, inparms.conv_ratio);
   bcm_BootAMGBuildSetSolverIt(bootamg_data, inparms.solver_it);
   bcm_AMGBuildSetAggMatchType(amg_data, inparms.matchtype);
   bcm_AMGBuildSetSweepNumber(amg_data, inparms.aggrsweeps);
   bcm_AMGBuildSetAggInterpType(amg_data, inparms.aggrtype);
   bcm_AMGBuildSetMaxLevels(amg_data, inparms.max_levels);
   bcm_AMGBuildSetCoarseSolver(amg_data, inparms.coarse_solver);
   bcm_AMGSetCycleType(amg_cycle, inparms.cycle_type);
   bcm_AMGSetRelaxType(amg_cycle, inparms.relax_type);
   bcm_AMGSetPreRelaxSteps(amg_cycle, inparms.prerelax_sweeps);
   bcm_AMGSetPostRelaxSteps(amg_cycle, inparms.postrelax_sweeps);


   printf("Composite Solver Type: %d\n",bcm_BootAMGBuildDataCompType(bootamg_data));
   printf("Max number of components: %d\n",bcm_BootAMGBuildDataMaxHrc(bootamg_data));
   printf("Desired Conv. Ratio (required in case of bootstrap): %e\n",bcm_BootAMGBuildDataDesRatio(bootamg_data));
   printf("Number of iteration to conv. rate test. (required in case of bootstrap): %d\n",bcm_BootAMGBuildDataCompIt(bootamg_data));
   printf("Matching Type: %d\n",bcm_AMGBuildDataAggMatchType(amg_data));
   printf("Aggregation sweeps: %d\n",bcm_AMGBuildDataSweepNumber(amg_data));
   printf("Aggregation Type: %d\n",bcm_AMGBuildDataAggInterpType(amg_data));
   printf("max_levels: %d\n",bcm_AMGBuildDataMaxLevels(amg_data));
   printf("coarse_solver: %d\n",bcm_AMGBuildDataCoarseSolver(amg_data));
   printf("cycle_type: %d\n",bcm_AMGApplyDataCycleType(amg_cycle));
   printf("relax_type: %d\n",bcm_AMGApplyDataRelaxType(amg_cycle));
   printf("prerelax_sweeps: %d\n",bcm_AMGApplyDataPreRelax(amg_cycle));
   printf("postrelax_sweeps: %d\n",bcm_AMGApplyDataPostRelax(amg_cycle));

   //if (argc > 1) fclose(fp);


   /* set maximum coarse size */

   //int maxcoarse=40*pow((double)num_rows,(double)1/3);
   int maxcoarse=10;
   bcm_AMGBuildSetMaxCoarseSize(amg_data, maxcoarse);
   printf("maxcoarsesize: %d\n",bcm_AMGBuildDataMaxCoarseSize(amg_data));

   /* No relaxation of smooth vectors in AMG building */
   int CRit=0;
   bcm_AMGBuildSetCRIterations(amg_data, CRit);
   printf("CRit: %d\n",bcm_AMGBuildDataCRIterations(amg_data));

   /* initialize num_grid_sweeps parameter w.r.t. the number of levels
      and the chosen cycle.
      We have to manage this in a setup routine after hierarchy building */
   
   num_grid_sweeps = (int *) calloc(inparms.max_levels-1, sizeof(int));
   for(i=0; i<inparms.max_levels-1; i++) num_grid_sweeps[i]=1;
   switch(inparms.cycle_type)  {
   case 1: /* H-cycle */
     {
       for(i=0; i<inparms.max_levels-1; i++) {
	 j=i%2; /*step is fixed to 2; it can be also different */
	 if(j==0) num_grid_sweeps[i]=2; /* if num_grid_sweeps is 2, we are using a hybrid V-W cycle */
       }
     }
     break;
   case 2: /* W-cycle */
     {
       for(i=0; i<inparms.max_levels-2; i++) num_grid_sweeps[i]=2;
     }
     break;
   }
   bcm_AMGApplyDataGridSweeps(amg_cycle)=num_grid_sweeps;


   /* set arbitrary initial (smooth) vector: generally
      we use unitary vector or random vector */

   w=bcm_VectorCreate(num_rows);
   bcm_VectorInitialize(w);
   int w_size=bcm_VectorSize(w);
   printf("wsize: %d\n",w_size);
   bcm_VectorSetConstantValues(w,1.0); 

   bcm_AMGBuildDataSmoothVector(amg_data)= w;

   /* start bootstrap process */

   printf("Bootstrap starting \n");

   bcm_BootAMG *boot_amg;
   double time1=time_getWallclockSeconds();
   boot_amg=bcm_Bootstrap(bootamg_data,amg_cycle);
   double time2=time_getWallclockSeconds()-time1;

   printf("Bootstrap ended\n");

   bcm_AMGHierarchy **Harray;

   printf("Number of components:  %d\n", bcm_BootAMGNHrc(boot_amg));
   printf("Estimated convergence %e \n", bcm_BootAMGEstRatio(boot_amg));
   printf("Information on the Components\n");
   Harray=bcm_BootAMGHarray(boot_amg);
   double avgcmpx=0;
   double avgwcmpx=0;
   double avgnumlev=0;
   for(k=0;k<bcm_BootAMGNHrc(boot_amg); k++) {
     printf("Component:  %d\n", k);
     printf("Number of levels:  %d\n", bcm_AMGHierarchyNumLevels(Harray[k]));
     printf("Operator cmplx for V-cycle %e \n", bcm_AMGHierarchyOpCmplx(Harray[k]));
     printf("Operator cmplx for W-cycle %e \n", bcm_AMGHierarchyOpCmplxW(Harray[k]));
     avgcmpx=avgcmpx+ bcm_AMGHierarchyOpCmplx(Harray[k]);
     avgwcmpx=avgwcmpx+ bcm_AMGHierarchyOpCmplxW(Harray[k]);
     avgnumlev=avgnumlev+ bcm_AMGHierarchyNumLevels(Harray[k]);
   }
   printf("Wall Clock Time for Building:  %e\n", time2); 
   printf("Wall Clock Time for building:  %e WT per million of laplacian nnz: %e\n", time2, time2/(nnz/1.0E+06));
   double timeall=time2;

   avgcmpx=avgcmpx/bcm_BootAMGNHrc(boot_amg);
   avgwcmpx=avgwcmpx/bcm_BootAMGNHrc(boot_amg);
   avgnumlev=avgnumlev/bcm_BootAMGNHrc(boot_amg);
   printf("AVG cmpx  %e\n", avgcmpx); 
   printf("AVG wcmpx  %e\n", avgwcmpx); 
   printf("AVG numlev  %e\n", avgnumlev); 
   
   /* generate vector including the nboot coordinates of the graph vertices */
   int n_hrc;
  
   n_hrc= bcm_BootAMGNHrc(boot_amg); //number of components of bootstrap amg/ n_hrc+1=dimensionality
   bcm_Vector **ws;
   ws=bcm_BootAMGSmoothVectors(boot_amg);
   double *ws_data;
   int lk;

   /* char filename[81];
   for(l=0; l<=n_hrc; l++)
   {
    sprintf(filename,"vectorw-%3.3d.mtx",l);
    if (ws[l]!=NULL) bcm_VectorPrint(ws[l],filename);
   } */
  

   double **W_data;
   double *w_svd, **v_svd;
   v_svd= (double **) calloc(n_hrc+1, sizeof(*v_svd));
   for (i=0; i<=n_hrc; i++) v_svd[i] = (double *) calloc(n_hrc+1, sizeof(double));
   w_svd= (double *) calloc(n_hrc+1, sizeof(double));
   W_data= (double **) calloc(num_rows, sizeof(*W_data));
   for (i=0; i<num_rows; i++) W_data[i] = (double *) calloc(n_hrc+1, sizeof(double));
   for(l=0; l<=n_hrc; l++)
   {
    ws_data=bcm_VectorData(ws[l]); 
    for(lk=0; lk<num_rows; lk++) 
    {
      W_data[lk][l]=ws_data[lk]; 
    }
   }
    
  
   /* computing SVD of W_data, left singular vectors will overwrite W_data while w_svd will include singular values */

  dsvd(W_data, num_rows, n_hrc+1, w_svd, v_svd); 
  //preparing for CLUSTERING with k-means-block
   printf("dim=%d\n",(n_hrc+1));
   double *cluster_centroid;
   double *X;
   int *cluster_assignment;
   X  = (double *) calloc((num_labels+1)*(n_hrc+1)*num_rows_original , sizeof(double));
   cluster_centroid  = (double *) calloc((n_hrc+1)*(num_labels+1)*inparms.NC, sizeof(double));
   cluster_assignment  = (int *) calloc(num_rows_original, sizeof(int));
   int id; 
   /*  
	with N= number of supernodes , ml=number of labels, dim= dimensions
	W_data is  (N+ml*N) x dim ,  X is N x (ml+1)*dim
	TOPOLOGICAL INFO
	The first N rows of W_data are the topological info of the supernodes in the id-th dimension,
	that become the (id*ml)th columns of X .
	LABELS
	(N+1)th row of W_data ------>  1st row  2nd column of X
	(N+2)th row of W_data ------>  1st row  3rd column of X
	*/
	
	//TOPOLOGICAL INFO
	for (id=0;id<=n_hrc;id++) //for each dimension
		for (i=0; i<num_rows_original;i++)  //for each row of X/supernode
			{
				X[i+id*num_rows_original*(num_labels+1)]=W_data[i][id]; //copy the topological info
			}		
	//LABELS
	int *A_I,*A_J,*A_original_I;
	A_original_I=bcm_CSRMatrixI(A_original);
	A_I=bcm_CSRMatrixI(A);
	A_J=bcm_CSRMatrixJ(A);
	int original_edges ,label_edges,label_node;
	for (i=0; i<num_rows_original;i++)  //for each row of X/supernode
	{
		k=0; //let's explore all the columns of X
		for (id=0;id<=n_hrc;id++) //for each dimension	
		{
			k++; //jump one	= topological info 
			original_edges = A_original_I[i+1] - A_original_I[i]; //number of connections/edges  for that supernode in the original graph 
        	for (j=original_edges; j < (original_edges+num_labels);j++) //for each new connection of the extended graph 
			{
				l=A_J[A_I[i]+j]; //the supernode is connected with this additional node/label
				X[i+k*num_rows_original]=W_data[l][id]; //take the l-th row as a coordinate of the supernode 
				k++;	// go to the next column
			}
				
		}
   	}
   
   
   printf("starting clustering\n");
   time1=time_getWallclockSeconds();
 
   int *KMAX;
   KMAX= (int *) calloc(num_rows_original, sizeof(int));
   int im;
   double mod;                      
   int *idx;
   idx  = (int *) calloc(inparms.NC, sizeof(int));
   int ii,jj; 
   double MODMAX=-DBL_MAX; 
   double entr;	
   int max_iterations;
   int iter;	
   for (im=0;im<100;im++)
   {
   
   max_iterations=10;
   iter=0; 
   
   /*INITIALIZING CENTROIDS before calling k_means_block
    centroids is NC x nb*dim, ordered for columns = 
	[x_topo_1,..,x_topo_NC,x_label1_node1,..,x_label1_nodeNC,..]
    (with N= number of supernodes , ml=number of labels, dim=dimensions, nb=ml+1) */  
    	//generation of NC random indexes
   		for(ii=0; ii<inparms.NC; ii++)  // for each cluster 
   			{
       			idx[ii]=random();   //compute a random number
       			idx[ii]=idx[ii]%num_rows_original; //idx[ii] must contain a number in [0;N-1] since we have N nodes
   				//    printf("printing random idx %d\n", idx[lk]);
   			}
   		//Assuming NC random rows of X as the cluster centroids
   		for (jj=0;jj<(num_labels+1)*(n_hrc+1);jj++) //for each column of centroid
   			{	  
			   	for(ii=0; ii<inparms.NC; ii++)       //for each cluster/raw of centroids                       
   					cluster_centroid [ii+jj*inparms.NC] = X[idx[ii]+jj*num_rows_original]; // cluster centroid= (idx-th) random raw of X
   			}    	
   			
   //           X,      NC  ,   nb=ml+1  , max_iter ,      dim ,   N               ,   indices    ,      centroids ,     iter
   kmeansblock(X,inparms.NC,num_labels+1,max_iterations,n_hrc+1,num_rows_original,cluster_assignment,cluster_centroid,&iter);
   time2=time_getWallclockSeconds()-time1;
   timeall=timeall+time2;
  
   mod= modularity(A_original,cluster_assignment,inparms.NC);
   
   
   if(mod > MODMAX)
            { 
			for(i=0; i<num_rows_original; i++) 
			{
			 KMAX[i]=cluster_assignment[i];
			}
            MODMAX=mod;
            }
    
   }
   entr=entropy(metadata, KMAX,inparms.NC );
   printf("MAX modularity=%e\n",MODMAX);
   printf(" entropy=%e \n",entr);
   FILE *fp3;
   fp3 = fopen("clustering_labels_AttributedGraph.mtx", "w");
   printf("printing kmeans output \n");
   fprintf(fp3, "NC=%d\t", inparms.NC);
   fprintf(fp3, "MOD=%e\t", MODMAX);
   fprintf(fp3, "ENT=%e\n",entr );
   for(i=0; i<num_rows_original; i++) 
   {
    fprintf(fp3, "%d \t %d\n", i, KMAX[i]);
   }
   fclose(fp3); 

   free(X);
   free(idx);
   free(cluster_centroid);
   free(cluster_assignment);
   free(KMAX);

#ifdef DUMP_HIER
   dump_on_file("BCM-pre-",boot_amg);
#endif

   free(num_grid_sweeps);
   free_inparms(&inparms);
   bcm_AMGApplyDataDestroy(amg_cycle);
   bcm_BootAMGBuildDataDestroy(bootamg_data);
   bcm_BootAMGDestroy(boot_amg);
   bcm_CSRMatrixDestroy(A);
   bcm_CSRMatrixDestroy(A_original);
   bcm_VectorDestroy(w);
   bcm_VectorDestroy(DIAGLA);
   bcm_CSRMatrixDestroy(B);
   bcm_CSRMatrixDestroy(E);
   bcm_CSRMatrixDestroy(ET);
   bcm_CSRMatrixDestroy(LA);
   bcm_CSRMatrixDestroy(LASPD);
   mge_MetadataDestroy(metadata);
   for (i=0; i<=n_hrc; i++) free(v_svd[i]);
   free(v_svd);
   free(w_svd);
   for (i=0; i<num_rows; i++) free(W_data[i]);
   free(W_data);

   return (0);
}

