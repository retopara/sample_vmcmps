/********************************************************************************************************
*   FILE NAME:  AdjustTemper.cpp
*
*   FUNCTION:   AdjustTemper
*   PURPOSE:    Adjust beta distribution(i.e. Temperature distribution) of the processes
*  
*   INPUT:      beta[], the temperature values of replicas
*               exch_acpt_count[], exchange temperature acceptance count
*               irec, the present AdjustTemper loop count
*               nrec, AdjustTemper total loop times
*               PARASIZE, the degree of parallelism, i.e, the number of threads
*   OUTPUT:     modified beta[] values
*********************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "remc.h"
#include "toolbox.h"


extern void AdjustTemper(double *beta, int *exch_acpt_count, const int *irec, const int *nrec)
{
    //internal definitions
    //const *irec, *nrec, *PARASIZE;
  	double weigh_scale = 0.7f;
    double lamda, sum_tmp, wk_tmp0;
    double *beta_tmp;

    static double *wb;
    static double wk;

    int i;
    //internal definitions end

    //malloc
    wb = (double *)malloc(sizeof(double)*PARASIZE);
    if(wb == NULL) printf("wb malloc fail!\n");

    beta_tmp = (double *)malloc(sizeof(double)*PARASIZE);
    if(beta_tmp == NULL) printf("beta_tmp malloc fail!\n");
    //malloc end

    //initialzing
    if(*irec == 0)
    {
        wk = 0;
        for(i=0; i<PARASIZE; i++)
        {
            wb[i] = 0;
        }
    }
    //initialzing end

    wk_tmp0 = (double)Minvec(PARASIZE-1, exch_acpt_count);  //caution! is zero?
#ifdef DEBUG
	if(myid == 0)
	{
		printf("wk_tmp0 = %f\n", wk_tmp0);
	}
#endif
		
    wk += (double)wk_tmp0;
    for(i=0; i<PARASIZE; i++)
    {
        wb[i] += wk_tmp0 * beta[i];
    }

	sum_tmp = 0.0f;
	for(i=1; i<PARASIZE; i++)
	{
		sum_tmp += (beta[i] - beta[i-1]) * exch_acpt_count[i-1];
	}
    lamda = (beta[PARASIZE-1] - beta[0]) / sum_tmp;

	beta_tmp[0] = beta[0];
	for(int i=1; i<PARASIZE; i++)
	{
		beta_tmp[i] = beta_tmp[i-1] + (1.0f-weigh_scale + weigh_scale * exch_acpt_count[i-1] * lamda) * (beta[i] - beta[i-1]);
	}

    for(i=0; i<PARASIZE; i++)
	{
		beta[i] = beta_tmp[i];
	}

    //free
    free(beta_tmp);
    //free end

    if(*irec == (*nrec-1))
    {
        for(i=0; i<PARASIZE; i++)
        {
            beta[i] = wb[i]/wk;
        }
        //free
        free(wb);
        //free end
    }
	return;
}
