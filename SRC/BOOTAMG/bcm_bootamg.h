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
#ifndef BCM_BOOTAMG_H_
#define BCM_BOOTAMG_H_
#include "bcm_amg.h"
#include "bcm_util.h"

/*--------------------------------------------------------------------------
 * bcm_BootAMG includes data defining number of hierarchies coming from bootstrap process
 * and estimated convergence rate of a corresponding solver
 *--------------------------------------------------------------------------*/

typedef struct
{
   /* data generated in the setup phase */

   bcm_AMGHierarchy **H_array; /*array of AMG hierarchies */
   int              n_hrc; /* number of hierarchies  */
   double           estimated_ratio; /* estimated convergence ratio obtained after built*/
   bcm_Vector       **SMOOTH_array; /* array of smooth vectors */

} bcm_BootAMG;

/*--------------------------------------------------------------------------
 * Accessor functions for the bcm_BootAMG structure
 *--------------------------------------------------------------------------*/

/* data generated by the setup phase */
#define bcm_BootAMGHarray(boot_amg) ((boot_amg)->H_array)
#define bcm_BootAMGNHrc(boot_amg) ((boot_amg)->n_hrc)
#define bcm_BootAMGEstRatio(boot_amg) ((boot_amg)->estimated_ratio)
#define bcm_BootAMGSmoothVectors(boot_amg) ((boot_amg)->SMOOTH_array)

/*--------------------------------------------------------------------------
 * bcm_BootAMGBuildData for building hierarchies by the bootstrap process
 *--------------------------------------------------------------------------*/

typedef struct
{

   /* setup params for bootstrap process */
   int      max_hrc; /* max number of hierarchies to be built*/
   double   conv_ratio; /* desired convergence ratio */
   int      solver_type; /* type of composition for applying composite solver */
   int      mvsolver_type; /* type of composition for applying composite solver in MVbootstrap*/
  /* solver_type =0 -> multiplicative composition
   * solver_type =1 -> symmetrized multiplicative composition
   * solver_type =2 -> additive composition */
   int      solver_it; /* number of iterations to be applied for conv. ratio estimating */
   int      num_vec; /* number of smooth vectors used for building multi-vector AMG */
   int      mvagg_type; /* type of aggregation used for building multi-vector AMG */

  /* setup params per each AMG component */
   bcm_AMGBuildData *amg_data;

} bcm_BootAMGBuildData;

/*--------------------------------------------------------------------------
 * Accessor functions for the bcm_BootAMGBuildData structure
 *--------------------------------------------------------------------------*/

#define bcm_BootAMGBuildDataMaxHrc(bootamg_data) ((bootamg_data)->max_hrc)
#define bcm_BootAMGBuildDataDesRatio(bootamg_data) ((bootamg_data)->conv_ratio)
#define bcm_BootAMGBuildDataCompType(bootamg_data) ((bootamg_data)->solver_type)
#define bcm_BootAMGBuildDataCompIt(bootamg_data) ((bootamg_data)->solver_it)
#define bcm_BootAMGBuildDataCompData(bootamg_data) ((bootamg_data)->amg_data)
#define bcm_BootAMGBuildDataNumVec(bootamg_data) ((bootamg_data)->num_vec)
#define bcm_BootAMGBuildDataMVAggType(bootamg_data) ((bootamg_data)->mvagg_type)
#define bcm_BootAMGBuildDataMVCompType(bootamg_data) ((bootamg_data)->mvsolver_type)

/* bcm_bootamg.c */
void * bcm_BootAMGBuildDataInitialize();
int bcm_BootAMGBuildSetMaxHrc(void *data, int  max_hrc);
int bcm_BootAMGBuildSetDesiredRatio(void *data, double con_ratio);
int bcm_BootAMGBuildSetSolverType(void  *data, int solver_type);
int bcm_BootAMGBuildSetMVSolverType(void  *data, int mvsolver_type);
int bcm_BootAMGBuildSetSolverIt(void  *data, int solver_it);
int bcm_BootAMGBuildSetNumVec( void *data, int   num_vec);
int bcm_BootAMGBuildSetMVAggType( void *data, int   mvagg_type);
int bcm_BootAMGBuildDataDestroy(void *data);
void * bcm_BootAMGCreate(int max_hrc);
int bcm_BootAMGInitialize();
int bcm_BootAMGHarrayDestroy(int nhrc, int start, bcm_AMGHierarchy **Harray);
int bcm_BootAMGDestroy(bcm_BootAMG *boot_amg);

/* bcm_bootstrap.c */
bcm_BootAMG * bcm_Bootstrap(bcm_BootAMGBuildData *bootamg_data, bcm_AMGApplyData *amg_cycle);
int bcm_InnerIterations(bcm_BootAMGBuildData *bootamg_data, bcm_BootAMG *boot_amg,
			bcm_AMGApplyData *amg_cycle, int solver_type);

/* bcm_gcycle.c */
int bcm_GAMGCycle(int k, bcm_BootAMGBuildData *bootamg_data,
		  bcm_BootAMG *boot_amg, bcm_AMGApplyData *amg_cycle, bcm_Vector **Rhs, bcm_Vector **Xtent, int l);

/* bcm_boot_prec.c */
int bcm_PrecApply(bcm_BootAMGBuildData *bootamg_data, bcm_BootAMG *boot_amg, bcm_AMGApplyData *amg_cycle, bcm_Vector *rhs, bcm_Vector *x);

/* bcm_inneritkcycle */
int bcm_inneritkcycle(int kh, bcm_Vector *x_vector, bcm_Vector *rhs_vector,
                  bcm_BootAMGBuildData *bootamg_data,
                  bcm_BootAMG * boot_amg, bcm_AMGApplyData *amg_cycle, double rtol, int l);

/* bcm_kc_apply.c */
int bcm_KCApply(int kk, bcm_BootAMGBuildData *bootamg_data,
                bcm_BootAMG *boot_amg, bcm_AMGApplyData *amg_cycle,
                int l, bcm_Vector *rhs, bcm_Vector *x);

#endif