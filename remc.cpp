/***********************************************************************************************
 *	 AUTHOR:		 Zhen Wang(USTC LQCC)
 *							 wangzhen.ustc@gmail.com
 *	 
 *   FILE NAME:  REMC.cpp
 *   FUNCTION:   REMC
 *
 *   PURPOSE:    The entrance of Replica-exchange Monte Carlo program
 *
 *   INPUT:      fvec, the value of target function
 *               r2, the error
 *               betamin, the min beta(max temperature) allowed for the replicas
 *               betamax, the max beta(min temperature) allowed for the replicas
 *               n, the number of variables
 *               nswp, nexch, nrec, loop times of CSA, ExchTemper, AdjustTemper
 *               (*target_func), function that calculate fvec value from x
 *               (*x_maker), function that generate new x[] according to csa acceptance rate
 *   OUTPUT:     yb, and r2
 ***********************************************************************************************/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "remc.h"
#include "toolbox.h"
#include "screenio.h"
#include "fileio.h"
#include "share.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
using namespace std;

extern void REMC(double *fvec_min, double *r2_min, const double *tol, const double *betamin, 
		const double *betamax, const int *n, const int *nswp, const int *nexch, 
		const int *nrec, const int *nequi1, const int *nequi2, double *x_min, 
		double (*target_func)(double *x_new), 
		void (*x_maker)(const double *PI, const int *n, const int *i_x, double *x, double *x_new, 
			double *sigma, int *csa_acpt_count, int *csa_count, int *seed_cau))
{
	//result
	MPI_Status status;
	int i, j;

	//global variables
	int csa_count = 0;
	double r2;
	double fvec_new, fvec;
	int nscheme;
	int seed_met, seedx, seed_cau;
	//global variables end

	//global pointers
	int **exbeta;
	int *seed;
	int *csa_acpt_count;
	int *pro_index, *beta_index;
	double *x, *sigma, *x_new;
	double *beta;
	double *fvec_recv;
	//	double *localmin_in, *localmin_out;
	struct reduce{
		double value;
		int rank;
	}localmin_in, localmin_out;

	int *exch_count, *exch_acpt_count;
	//global pointers end

	//internal definitions
	int irec, iexch, iswp, iequi;
	int ischeme = 0;
	int conv_count = 0;
	int CONVMAX;
	//internal definitions end

	//ReplicaInit
	{
		//internal definitions

		int seed0 = -2340;
		double sigmainit;
		//internal definitions end

		//file input
		//c++ version input
		ifstream ifs("./input/replicainit.input");
		if(!ifs)
		{
			cerr << "\nCan't find ./input/replicainit.input." << endl;
			abort();
		}
		string str;
		ifs >> PI;
		getline(ifs, str);

		ifs >> CONVMAX;
		getline(ifs, str);

		ifs >> nscheme;
		getline(ifs, str);

		ifs >> sigmainit;
		getline(ifs, str);
		ifs.close();

		if(myid == 0)
		{   
			cout << "\nPI=" << PI;
			cout << "\nCONVMAX=" << CONVMAX;
			cout << "\nnscheme=" << nscheme;
			cout << "\nsigmainit=" <<sigmainit << endl;

			//c++ version input
		}
		ifs.close();
		/*
			 MPI_Bcast(&PI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			 MPI_Bcast(&sigmainit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			 MPI_Bcast(&CONVMAX, 1, MPI_INT, 0, MPI_COMM_WORLD);
			 MPI_Bcast(&nscheme, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef DEBUG
printf("nscheme = %d\n", nscheme);
#endif
*/
		//file input end

		//malloc
		fvec_recv = (double *)malloc(sizeof(double)*2);
		if(fvec_recv == NULL) printf("fvec_recv malloc fail!\n");

		exch_count = (int *)malloc(sizeof(int)*PARASIZE);
		if(exch_count == NULL) printf("exch_count malloc fail!\n");

		exch_acpt_count = (int *)malloc(sizeof(int)*PARASIZE);
		if(exch_acpt_count == NULL) printf("exch_acpt_count malloc fail!\n");

		csa_acpt_count = (int *)malloc(sizeof(int)*(*n));
		if(csa_acpt_count == NULL) printf("csa_acpt_count malloc fail!\n");

		exbeta = (int **)malloc(sizeof(int *)*(nscheme));
		if(exbeta == NULL) printf("exbeta[][] malloc fail!\n");
		for(i=0; i<nscheme; i++){
			exbeta[i] = (int *)malloc(sizeof(int)*PARASIZE);
			if(exbeta[i] == NULL) printf("exbeta[%d] malloc fail!\n", i);
		}

		x_new = (double *)malloc(sizeof(double)*(*n));
		if(x_new == NULL) printf("x_new[] malloc fail!\n");

		x = (double *)malloc(sizeof(double)*(*n));
		if(x == NULL) printf("x[] malloc fail!\n");

		x_min = (double *)malloc(sizeof(double)*(*n));
		if(x_min == NULL) printf("x_min[] malloc fail!\n");

		pro_index = (int *)malloc(sizeof(int)*PARASIZE);
		if(pro_index == NULL) printf("pro_index[] malloc fail!\n");

		beta_index = (int *)malloc(sizeof(int)*PARASIZE);
		if(beta_index == NULL) printf("beta_index[] malloc fail!\n");

		sigma = (double *)malloc(sizeof(double)*(*n));
		if(sigma == NULL) printf("sigma[] malloc fail!\n");

		beta = (double *)malloc(sizeof(double)*PARASIZE);
		if(beta == NULL) printf("beta[] malloc fail!\n");

		seed = (int *)malloc(sizeof(int)*PARASIZE);
		if(seed == NULL) printf("seed[] malloc fail!\n");

		//malloc end

		//initializing 
		for(i=0; i<*n; i++){
			csa_acpt_count[i] = 0;
		}
		for(i=0; i<*n; i++){
			sigma[i] = sigmainit;
		}

		//initializing seed
		for(i=0; i<PARASIZE; i++){
			seed[i] = int(-59370 * Ran(&seed0));
		}
		for(i=0; i<PARASIZE; i++){
			for(j=0; j<PARASIZE; j++){
				if((i!=j)&&(seed[i] == seed[j])){
					printf("different processes have same seeds!\n");
				}
			}
		}
		seed_met = seed[myid];
		seedx = seed_met - 1234;
#ifdef DEBUG
		if(myid == 0){
			printf("seedx = %d\n", seedx);
		}
#endif
		seed_cau = seedx;
		//initializing seed end    

		//initialzing x[]
#ifdef DEBUG
		printf("seedx = %d\n", seedx);
#endif
		//    xinput(*n, x);
		for(i=0; i<*n; i++){
			x[i] = Ran(&seedx);
		}
		//initializing x[] end

		//initializing beta[]
		if(PARASIZE == 1)
		{
			beta[0] = *betamin;
		}
		else
		{
			for(i=0; i<PARASIZE; i++)
			{
				beta[i] = *betamin + i * (*betamax - *betamin)/(PARASIZE-1);
			}
		}

		for(i=0; i<PARASIZE; i++)
		{
			pro_index[i] = i;
			beta_index[i] = i;
#ifdef DEBUG
			printf("pro_index[%d] = %d\t", i, pro_index[i]);
#endif
		}
		//initializing beta[] end

		//initializing exbeta[]
		for(i=0; i<nscheme; i++){
			for(j=0; j<PARASIZE; j++){
				exbeta[i][j] = -1;
			}
		}
		for(i=0; i<nscheme; i++){
			for(j=0; j<PARASIZE; j++){
				int k = i + nscheme*j;
				if(k+1 <= PARASIZE) exbeta[i][j] = k;
#ifdef DEBUG
				printf("exbeta[%d][%d] = %d\t", i, j, exbeta[i][j]);
#endif
			}
		}
		//initializing exbeta[] end
	}//ReplicaInit end

	/*
		 ReplicaInit(n, betamin, betamax, 
		 x, x_new, x_min, exbeta, 
		 seed, &seed_met, &seedx, &seed_cau, 
		 beta, pro_index, beta_index, &nscheme, 
		 &PI, fvec_recv, localmin_in, localmin_out, 
		 csa_acpt_count, sigma, exch_count,
		 exch_acpt_count, &CONVMAX);
		 */

	//initialize the value of fvec
#ifdef DEBUG
	for(i=0; i<*n; i++){
		printf("x[%d] = %f\n", i, x[i]);
	}
#endif
	fvec = (*target_func)(x);

	//*fvec_min = fvec;//set initial value for fvec_min
	*fvec_min = 10000.0;//set initial value for fvec_min, set a huge value, this way is safer for fewer n_x

	if(myid == 0)
	{
		divider("The First Time");
		printf("\nthe initial beta is:\n");
		for(i=0; i<PARASIZE; i++)
		{
			printf("i = %d, beta = %f\n", i, beta[i]);
		}
		printf("\ninitial fvec = %f, fvec_min = %f\n", fvec, *fvec_min);
	}

	for (iequi=0; iequi<*nequi1; iequi++)
	{
		CSA(n, &PI, x_new, x, 
				x_min, sigma, csa_acpt_count, &fvec_new, 
				&fvec, fvec_min, 
				&csa_count, &beta[beta_index[myid]], &seed_met, &seed_cau, 
				target_func, x_maker);
	}

	for(irec=0; irec<*nrec; irec++)
	{
		if(irec != 0)
		{
			recflag_send(fvec_min, 10000.0);
			if(myid == 0) cout << "====================================\nFOR NEXT RECURSION, fvec_min = " << *fvec_min << endl;
		}
		for(i=0; i<PARASIZE; i++)
		{
			exch_acpt_count[i] = 1;
			exch_count[i] = 1;
		}
#ifdef DEBUG
		printf("run to exch loop!\n");
#endif
		for(iexch=0; iexch<*nexch; iexch++)
		{
			for(iswp=0; iswp<*nswp; iswp++)
			{
#ifdef DEBUG
				printf("run all the way down to BEFORE CSA!\n");
#endif
				CSA(n, &PI, x_new, x, 
						x_min, sigma, csa_acpt_count, &fvec_new, 
						&fvec, fvec_min, 
						&csa_count, &beta[beta_index[myid]], &seed_met, &seed_cau, 
						target_func, x_maker);
#ifdef DEBUG
				printf("myid = %d, fvec_min = %f\n", myid, *fvec_min);
#endif
			}
			ExchTemper(&nscheme, exch_count, exch_acpt_count, exbeta, 
					pro_index, beta_index, &fvec, fvec_recv, 
					beta);
		}

		if(myid == 0)
		{
			divider("Current Recursion");
			printf("\nirec = %d\n", irec);
#ifdef DEBUG
			for(i=0; i<PARASIZE; i++)
			{
				printf("i = %d, beta = %f\n", i, beta[i]);
			}
			putchar(10);
#endif
			printf("exch acpt rate is:\n");
			for(i=0; i<PARASIZE; i++)
			{
				printf("i = %d %8d %8d\n", i, exch_acpt_count[i], exch_count[i]);
			}
			printf("\npro_index and beta_index are:\n");
			for(i=0; i<PARASIZE; i++)
			{
				printf("i = %d %8d %8d\n", i, pro_index[i], beta_index[i]);
			}
		}
		AdjustTemper(beta, exch_acpt_count, &irec, nrec);

		for (iequi=0; iequi<*nequi2; iequi++)
		{
			CSA(n, &PI, x_new, x, 
					x_min, sigma, csa_acpt_count, &fvec_new, 
					&fvec, fvec_min, 
					&csa_count, &beta[beta_index[myid]], &seed_met, &seed_cau, 
					target_func, x_maker);
		}

		// find min //
		localmin_in.value = *fvec_min;
		localmin_in.rank = myid;

		MPI_Reduce(&localmin_in, &localmin_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
		MPI_Bcast(&localmin_out.rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//send x_min
		if(myid == localmin_out.rank){
			MPI_Send(x_min, *n, MPI_DOUBLE, 0, 80, MPI_COMM_WORLD);
		}
		if(myid == 0){
			MPI_Recv(x_min, *n, MPI_DOUBLE, localmin_out.rank, 80, MPI_COMM_WORLD, &status);
		}

		//MPI_Bcast(x_min, *n, MPI_DOUBLE, localmin_out.rank, MPI_COMM_WORLD);

		/*      //alternate way to create r2  
						fvec_now = localmin_out[1];
		 *r2 = fvec_last - fvec_now;
		 if(*r2 < 0) *r2 = -(*r2);
		 fvec_last = fvec_now;
		 */

		if(myid == 0)
		{
			cout << "\nirec = " << irec;
			cout << "\nmin_rank = " << localmin_out.rank;
			cout << "\nfvec_min = " << scientific << setprecision(9) << localmin_out.value << endl;

			if(irec == *nrec-1)
			{
				divider("x_min[] is:");
				for(i=0; i<*n; i++)
				{
					printf("x_min[%d] = %e\n", i, x_min[i]);
					//printf("\np[%d] = %f", i, x_min[i]*x_min[i]);
				}
				//        Output(x_min, localmin_out.rank, localmin_out.value);
			}
		}
		/*
		//judge convergence
		if(localmin_out.r2 < *tol){
		conv_count++;
		}
		else{
		conv_count = 0;
		}
		if(conv_count == CONVMAX){
//free
for(i=0; i<nscheme; i++){
free(exbeta[i]);
}
free(exbeta);
free(csa_acpt_count);
free(x_new);
free(pro_index);
free(beta_index);
free(x);
free(x_min);
free(sigma);
free(beta);
free(exch_count);
free(exch_acpt_count);
free(fvec_recv);
//free end

return;
}
*/
}//end for(irec=0; irec<*nrec; irec++)


//free
for(i=0; i<nscheme; i++){
	free(exbeta[i]);
}
free(exbeta);

free(csa_acpt_count);
free(x_new);
free(pro_index);
free(beta_index);
free(x);
free(x_min);
free(sigma);
free(beta);
free(exch_count);
free(exch_acpt_count);
free(fvec_recv);
//free end

return;
}
