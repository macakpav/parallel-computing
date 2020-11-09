/******************************************************************************
* FILE: omp_bug1.c
* DESCRIPTION:
*   This example attempts to show use of the parallel for construct.  However
*   it will generate errors at compile time.  Try to determine what is causing
*   the error.  See omp_bug1fix.c for a corrected version.
* AUTHOR: Blaise Barney  5/99
* LAST REVISED: 04/06/05
* Modified by Emil Loevbak
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N       50
#define CHUNKSIZE   5 

int main (int argc, char *argv[]) 
{
int i, chunk, tid=-1;
float a[N], b[N], c[N] ;

/* Some initializations */
for (i=0; i < N; i++)
  a[i] = b[i] = i * 1.0 ;
chunk = CHUNKSIZE ;
  {
    #pragma omp parallel for    \
    shared(a,b,c,chunk)         \
    private(i)                  \
    schedule(static,chunk)      \
    firstprivate(tid)
    for (i=0; i < N; i++)
    {
      if (tid==-1) tid = omp_get_thread_num() ;
      c[i] = a[i] + b[i] ;
      printf("Thread id= %d i= %d c[i]= %f, diff = %f\n", tid, i, c[i], 2*i-c[i]) ;
    }
  }  /* end of parallel for construct */

}

