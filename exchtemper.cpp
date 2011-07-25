/*******************************************************************************************
*   FILE NAME:  ExchTemper.cpp
*   FUNCTION:   ExchTemper
*
*   PURPOSE:    Exchange temperatures of neighbor replicas(same as exchanging comformations)
*
*   INPUT:      nscheme, the number of total exchange schemes
*               exch_count[], exchange times count of each temperature
*               exch_acpt_count[], exchange acceptance count of each temperature
*               exbeta[][], stores the exchange scheme
*               pro_index[], process index sort by beta rank
*               beta_index[], beta rank index sort by process id
*               fvec, the value of target function
*               fvec_recv[], on thread rank 0, for receiving fvec from other 2 thread
*               myid, the process id
*               PARASIZE, the degree of parallelism
*
*   OUTPUT:     The processes after the exchange
********************************************************************************************/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "remc.h"
#include "toolbox.h"

extern void ExchTemper(const int *nscheme, int *exch_count, int *exch_acpt_count, int **exbeta, 
                       int *pro_index, int *beta_index, double *fvec, double *fvec_recv, 
                       double *beta)
{
    MPI_Status status;
    int tag = 99;

    //const input
    //const *nscheme, *myid, *fvec;

    static int ischeme = 0;

    int tmppro, tmpbeta;
    int seed_exch = -5000; //??????????????????
    int ipair = 0;
    int low_beta, high_beta;
    double deltah;
    //internal definitions end

	if(ischeme > (*nscheme-1))
	{
		ischeme = 0;
	}
#ifdef DEBUG2
	if(myid == 0)
	printf("ischeme = %d\n", ischeme);
#endif

    while(exbeta[ischeme][ipair] != -1)
    {
        low_beta = exbeta[ischeme][ipair];
        high_beta = exbeta[ischeme][ipair] + 1;
		if(high_beta >= PARASIZE){
			break;
		}
#ifdef DEBUG
		if(myid == 0)
		{
			putchar(10);
			printf("low_beta = %d, high_beta = %d\n", low_beta, high_beta);
		}
#endif
        exch_count[low_beta]++;
#ifdef DEBUG
		printf("exch_count[%d] = %d\n", low_beta, exch_count[low_beta]);
#endif

        if(myid == pro_index[low_beta]){//PARA_PART
#ifdef DEBUG
			//printf("pro_index[low_beta] = %d\n", pro_index[low_beta]);
			printf("low_beta_myid = %d\n", myid);
			printf("fvec = %f\n", fvec);
#endif
        	MPI_Send(fvec, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
#ifdef DEBUG
        	printf("MPI_Send_low_beta complete!\n");
#endif

       } 

        if(myid == pro_index[high_beta]){//PARA_PART
            MPI_Send(fvec, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
#ifdef DEBUG
			printf("pro_index[high_beta] = %d\n", pro_index[high_beta]);
            printf("MPI_Send_high_beta complete!\n");
#endif
        }
        if(myid == 0){//PARA_PART
            MPI_Recv(&fvec_recv[0], 1, MPI_DOUBLE, pro_index[low_beta], tag, MPI_COMM_WORLD, &status);
#ifdef DEBUG
            printf("MPI_Recv_low_beta complete!\n");
#endif
            MPI_Recv(&fvec_recv[1], 1, MPI_DOUBLE, pro_index[high_beta], tag, MPI_COMM_WORLD, &status);
#ifdef DEBUG
            printf("MPI_Recv_high_beta complete!\n");
#endif
            deltah = (fvec_recv[0] - fvec_recv[1]) * (beta[high_beta] - beta[low_beta]);
#ifdef DEBUG
			printf("fvec_recv[1] = %f, fvec_recv[0] = %f\n", fvec_recv[1], fvec_recv[0]);
			printf("beta_high = %f, beta_low = %f\n", beta[high_beta], beta[low_beta]);
			printf("deltah = %f\n", deltah);
#endif

			if(MetroFlag(deltah, &seed_exch))
            {
                /* exchange beta index */
                tmpbeta = beta_index[pro_index[low_beta]];
                beta_index[pro_index[low_beta]] = beta_index[pro_index[high_beta]];
                beta_index[pro_index[high_beta]] = tmpbeta;
                
                /* exchange process index */
                tmppro = pro_index[high_beta];
                pro_index[high_beta] = pro_index[low_beta];
                pro_index[low_beta] = tmppro;

                /* exchange accept count ++ */
                exch_acpt_count[low_beta]++;                
            }
        }//end if(myid == 0)
        //BCAST TO ALL//PARA_PART
        ipair++;
    }//end while(exbeta[ischeme][ibeta] != -1)
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
    printf("MPI_exch_Barrier complete!\n");
#endif
    MPI_Bcast(beta_index, PARASIZE, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(pro_index, PARASIZE, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(exch_acpt_count, PARASIZE, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef DEBUG
    printf("MPI_exch_Bcast complete!\n");
#endif
    ischeme ++;

	return;
}
