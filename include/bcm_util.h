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
/*--------------------------------------------------------------------------
 * Define various functions
 *--------------------------------------------------------------------------*/
#ifndef BCM_UTIL_H_
#define BCM_UTIL_H_

#ifndef bcm_max
#define bcm_max(a,b)  (((a)<(b)) ? (b) : (a))
#endif
#ifndef bcm_min
#define bcm_min(a,b)  (((a)<(b)) ? (a) : (b))
#endif

#ifndef bcm_round
#define bcm_round(x)  ( ((x) < 0.0) ? ((int)(x - 0.5)) : ((int)(x + 0.5)) )
#endif



/* bcm_timer.c */
double time_getWallclockSeconds(void);
double time_getCPUSeconds(void);
double time_get_cpu_seconds_(void);
double time_get_cpu_seconds_(void);

/* bcm_random.c */
double  bcm_Rand();

/* bcm_vsd.c */
int dsvd(double **a, int m, int n, double *w, double **v);

/* bcm_quicksort.c */
void order(double v[],int first,int last);
#endif
