/************************************************************************
 *   FILE NAME:  XGenerator.cpp
 *   FUNCTION:   XGenerator
 *
 *   PURPOSE:    Generate x_new[] from x[], adjust accroding
 *               to CSA acceptance rate.
 *   INPUT:      PI, the value of PI
 *               n, the number of variables, i.e, the elements in x[]
 *               i_x, indicates which element in x[] need to be generated
 *               x[], the variables
 *               x_new[], the new x[] generated from x[]
 *               sigma[], the step length between x_new[] and x[]
 *               csa_acpt_count[], csa accpetance count
 *               csa_count, csa run time count
 *               seed_cau, the seed used to generate x_new[]
 *   OUTPUT:     x_new[]
 ************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <iostream>
#include "target.h"
#include "toolbox.h"

using namespace std;

extern void XGenerator(const double *PI, const int *n, const int *i_x, double *x, double *x_new, 
		double *sigma, int *csa_acpt_count, int *csa_count, int *seed_cau)
{
	//internal definitions
	const static double ACPT_MAX = 0.3f;
	const static double ACPT_MIN = 0.2f;
	const static double CONST1 = 2.0f;
	const static double CONST2 = 2.0f;
	const static int COUNT_MAX = 100;

	double ran_num, r0;
	double csa_rate;

	int i;
	//internal definitions end

	/*generate x_new[i]*/
	ran_num = Ran(seed_cau);
	r0 = tan(*PI * (ran_num - 0.5f));
	x_new[*i_x] = x[*i_x] + r0*sigma[*i_x];

	/* adjust sigma[] according to csa acceptance rate */
	if(*csa_count == COUNT_MAX)
	{
		for(i=0; i<*n; i++)
		{
			//			cout << "\n i = " << i;
			csa_rate = double(csa_acpt_count[i])/double(COUNT_MAX);
			if(csa_rate > ACPT_MAX){
				sigma[i] *= (1 + CONST1 * (csa_rate-ACPT_MAX) / (1 - ACPT_MAX));
				if(sigma[i] > pow(10.0, 3)) sigma[i] = pow(10.0, 3);
			}
			else if(csa_rate < ACPT_MIN){
				sigma[i] /= (1 + CONST2 * (ACPT_MIN-csa_rate) / ACPT_MIN);
				if(sigma[i] > pow(10.0, 3)) sigma[i] = pow(10.0, 3);
			}
			else{}
		}
		*csa_count = 0;
		for(i=0; i<*n; i++){
			csa_acpt_count[i] = 0;
		}
	}

	return;
}

