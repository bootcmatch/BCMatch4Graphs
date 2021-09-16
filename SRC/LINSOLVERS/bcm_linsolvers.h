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
#ifndef BCM_KSOLVERS_H_
#define BCM_KSOLVERS_H_

#include "bcm_bootamg.h"

/* bcm_boot_solver.c */
int bcm_boot_solver(bcm_BootAMGBuildData *bootamg_data, bcm_BootAMG *boot_amg,
		bcm_AMGApplyData *amg_cycle,
                bcm_Vector *rhs, bcm_Vector *x,  int solver_it, double RTOL);
/* csr_bootpcg.c */
int bcm_bootpcg(bcm_Vector *x_vector, bcm_Vector *rhs_vector, bcm_BootAMGBuildData *bootamg_data,
		bcm_BootAMG * boot_amg, bcm_AMGApplyData *amg_cycle, int precon, int max_iter,
		double rtol, int *num_iter, double *timetot);

/* csr_bootpcgnoflex.c */
int
bcm_bootpcgnoflex(bcm_Vector *x_vector, bcm_Vector *rhs_vector, bcm_BootAMGBuildData *bootamg_data,
                  bcm_BootAMG * boot_amg, bcm_AMGApplyData *amg_cycle, int precon, int max_iter,
                  double rtol, int *num_iter, double *timetot);


#endif
