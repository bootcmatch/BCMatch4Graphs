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
#include "bcm_linsolvers.h"

int 
bcm_bootpcgnoflex(bcm_Vector *x_vector, bcm_Vector *rhs_vector, bcm_BootAMGBuildData *bootamg_data,
            bcm_BootAMG * boot_amg, bcm_AMGApplyData *amg_cycle, int precon, int max_iter,
            double rtol, int *num_iter, double *timetot)
{

  int ierr=0;
  double rhs_norm ,delta0, delta_old, delta, eps = DBL_EPSILON;
  bcm_Vector *v_vector, *w_vector, *d_vector; 
  int num_dofs;
  
  int iter=0;
  double tau, alpha, beta;

  double l2_norm;
  double   zero=0.0, one=1.0;

 /* matrix data */
      bcm_AMGBuildData *amg_data;
      amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);

      bcm_CSRMatrix *A;
      A=bcm_AMGBuildDataCSRMatrix(amg_data);

      num_dofs=bcm_CSRMatrixNumRows(A);

      double time1, time2;
      *timetot=0.0;

/* Create and Initialize Workspaces */
      v_vector=bcm_VectorCreate(num_dofs);
      ierr=bcm_VectorInitialize(v_vector);
      w_vector=bcm_VectorCreate(num_dofs);
      ierr=bcm_VectorInitialize(w_vector);
      d_vector=bcm_VectorCreate(num_dofs);
      ierr=bcm_VectorInitialize(d_vector);

  bcm_VectorCopy(rhs_vector,w_vector);
  /* sparse-matrix vector product: --------------------------*/
  ierr=bcm_CSRMatrixMatvec(one,A,x_vector,zero,v_vector);
  /* compute residual  */
  ierr=bcm_VectorAxpy(-one,v_vector,w_vector);

  delta0 = bcm_VectorNorm(w_vector);
  rhs_norm = bcm_VectorNorm(rhs_vector);
  printf("residual: %e\n", delta0);

  if (delta0 <= eps * rhs_norm)
    {
      *num_iter = 0;
      return ierr;
    }

     if(precon)
     {
 /* apply preconditioner to w */
      time1=time_getWallclockSeconds();
      ierr=bcm_PrecApply(bootamg_data, boot_amg, amg_cycle, w_vector, d_vector);
      time2=time_getWallclockSeconds()-time1;
      *timetot=*timetot+time2;
     }

  delta_old = bcm_VectorInnerProd(w_vector, d_vector);
  if (delta_old < 0.e0)
    {
      printf("\n ERROR: indefinite preconditioner in cg_iter_coarse: %e\n", delta_old);
      return -1;
    }
    
loop:

  /* sparse-matrix vector product: --------------------------*/

  ierr=bcm_CSRMatrixMatvec(one,A,d_vector,zero,v_vector);

  tau = bcm_VectorInnerProd(v_vector, d_vector);
  if (tau <= 0.e0) 
    {
      printf("\n ERROR: indefinite matrix in cg_iter_coarse: %e\n", tau);
      return -1;                               
    }

  alpha = delta_old/tau;
  /* update solution  */
  ierr=bcm_VectorAxpy(alpha,d_vector,x_vector);
  iter++;
  /* update residual  */
  ierr=bcm_VectorAxpy(-alpha,v_vector,w_vector);

  l2_norm = bcm_VectorNorm(w_vector);

  bcm_VectorSetConstantValues(v_vector, zero);
      if(precon)
      {
 /* apply preconditioner to w */
      time1=time_getWallclockSeconds();
      ierr=bcm_PrecApply(bootamg_data, boot_amg, amg_cycle, w_vector, v_vector);
      time2=time_getWallclockSeconds()-time1;
      *timetot=*timetot+time2;
      }

  delta = bcm_VectorInnerProd(v_vector, w_vector);
  if (delta <= 0.e0)
    {
      printf("\n ERROR: indefinite preconditioner in cg_iter_coarse: %le; \n", delta);
      return -1;
    }
    
  beta = delta /delta_old;
  delta_old = delta;

  /* update direction  */
  ierr=bcm_VectorAxpy(beta,d_vector,v_vector);
  bcm_VectorCopy(v_vector,d_vector);

  printf("bootpcgnoflex iteration: %d;  residual: %e, relative residual: %e\n", iter, l2_norm, l2_norm/delta0);
  if (l2_norm > rtol * delta0 && iter < max_iter) goto loop;

  *num_iter = iter;
  return ierr;
}
