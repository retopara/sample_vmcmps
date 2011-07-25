/***********************************************************************************************
 *  FILE NAME:  CSA.cpp
 *  FUNCTION:   CSA
 *
 *  PURPOSE:    Implement Constrained simulated annealing(CSA)
 *
 *  INPUT:      n, the number of variables
 *              PI, the value of PI
 *              x_new[], the new variables generated by X_Maker
 *              x[], the variables from last CSA procedure
 *              x_min[], the relevant variables of the min fvec
 *              sigma[], the step length between x_new[] and x[]
 *              csa_acpt_count[], csa acceptance count
 *              fvec_new, new fvec value, calculated from x_new[]
 *              fvec, calculated from x[]
 *              fvec_min, min fvec value, on the thread, from program start to now
 *              r2, represent the accuracy of fvec, the error
 *              r2_min, the relevant r2 of fvec_min
 *              csa_count, csa run time count
 *              beta_myid, beta value of current thread
 *              seed_met, seed for metropolis
 *              int *seed_cau, seed for generating x_new[]
 *              (*target_func), function that calculate fvec value from x
 *              (*x_maker), function that generate new x[] according to csa acceptance rate
 *
 *  OUTPUT:     Updated value after one sweep of CSA
 ***********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "remc.h"
#include "toolbox.h"
#include "fileio.h"

#include <iostream>
#include <iomanip>
using namespace std;

extern void CSA(const int *n, const double *PI, double *x_new, double *x,
		double *x_min, double *sigma, int *csa_acpt_count, double *fvec_new, 
		double *fvec, double *fvec_min,  
		int *csa_count, const double *beta_myid, int *seed_met, int *seed_cau, 
		double (*target_func)(double *x_new), 
		void (*x_maker)(const double *PI, const int *n, const int *i_x, double *x, double *x_new, 
			double *sigma, int *csa_acpt_count, int *csa_count, int *seed_cau))
{
	//internal definitions    
	int i, j;

	//internal definitions end

	(*csa_count)++; //csa run time +1//!!!!!!!!!!!!!must put ahead to implement same as fortran

	for(i=0; i<*n; i++){

		/*restore x_new[] from x[]*/
		for(j=0; j<*n; j++)
		{
			x_new[j] = x[j];
		}
#ifdef DEBUG
		if(myid == 0)
		{
			printf("\nx_new[%d] = %16.15f\n", i, x_new[i]);
			printf("seed_cau = %d\n", *seed_cau);
		}
#endif
//    (*x_maker)(PI, n, &i, x, x_new, 
//        sigma, csa_acpt_count, csa_count, seed_cau);
#ifdef DEBUG
		if(myid == 0)
		{
			printf("x_new[%d] = %16.15f\n", i, x_new[i]);
			printf("seed_cau = %d\n", *seed_cau);
		}
#endif
		*fvec_new = (*target_func)(x_new);
#ifdef DEBUG
		if(myid == 0)
		{
			printf("x_new[%d] = %16.15f\n", i, x_new[i]);
		}
#endif

		/* using Metropolis to judge whether to accept */
		if(MetroFlag((*fvec_new - *fvec) * (*beta_myid), seed_met))
		{
			csa_acpt_count[i]++;
			for(j=0; j<*n; j++){
				x[j] = x_new[j];
			}
			*fvec = *fvec_new;

			if((*fvec) < (*fvec_min))   //if fvec_new is smaller than fvec_min, update all the min values
			{
#ifdef DEBUG
				printf("\nrun up to here. myid = %d", myid);
				printf("fvec_min should be changed\nfvec = %f, fvec_new = %f, fvec_min = %f\n", *fvec, *fvec_new, *fvec_min);
#endif
				*fvec_min = *fvec;
				for(j=0; j<*n; j++){
					x_min[j] = x[j];
#ifdef DEBUG
					if(myid == 0) printf("\nin csa x_min[%d] = %e", i, x_min[i]);
#endif
				}
			}   //if(*fvec < *fvec_min) 
		}   //if(MetroFlag)

		//DISPLAY
		etotoutput(*fvec);
		//DISPLAY

	}   //for(i=1; i<*n; i++)
}