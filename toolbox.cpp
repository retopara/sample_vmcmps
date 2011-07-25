/*******************************************************************************************
*   FILE NAME:  Toolbox.cpp
*   PURPOSE:    Small and useful fuctions
*
*   FUNCTION:   MetroFlag(), Implement metropolis rule, return exchange flag.
*               Max(), return the max number of 2 double.
*               MaxInt(), return the max number of 2 int.
*               Min(), return the min number of 2 double.
*               CalDenmat(), calculate density matrix from a vector and return
*               Ptrace3(), this function is written specifically to partially trace a 3-partite
*               by anyone particle, turning it into a 2-partite 2-level system.
*               Ptrace4(), this function is written specifically to partially trace a 4-partite
*               by anyone particle, turning it into a 3-partite 2-level system.
********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "toolbox.h"

extern int MetroFlag(double delta, int *seed)
{
    int acpt = 0;
    if(delta <= 0.0){   //fortran .lt. seems equal to <=
        acpt = 1;
    }
    else{
        double tmp = Ran(seed);
        if(exp(-delta) > tmp) acpt = 1;
    }
    return acpt;
}

extern int Minvec(int size, int *vec)
{
    int i;
	int minval = vec[0];
	for(i=0; i<size; i++)
	{
		if( vec[i] < minval )
		{
			minval = vec[i];
		}
	}
	return minval;
}
/*
extern void DecConvert(const int NPAR, const int NLEV, const int num, int *veck)
{
	//internal definitions
	int dec, i;
	double num_tmp = num;

	for(i=0; i<NPAR; ++i)
	{
		dec = pow(double(NLEV), NPAR-i-1);
		veck[i] = int(num_tmp/dec);
		num_tmp -= veck[i]*dec;
	}
}

extern void NewIntVec(const int row, const int col, int** vec)
{
	vec = new int*[row];
	for(int i=0; i<row; i++)
	{
		vec[i] = new int[col];
	}
}
*/
