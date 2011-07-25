/*********************************************************************************
*   FILE NAME:  Normalize.cpp
*   FUNCTION:   Normalize
*
*   PURPOSE:    Normalize x[]
*
*   INPUT:      q[], cr[][], ci[][], the new coefficients copied from x[]
*               NC, NP, numbers of coefficients c and p
*
*   OUTPUT:     Normalized q[], cr[][], ci[][]
********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include "remc.h"
#include "target.h"
#include "screenio.h"

//void Normalize(const int i_x, const int nx, double *x)
void Normalize(const int nx, double *x)
{
	double max = x[0];
	int maxi = 0;
	for(int i=0; i<nx; i++)
	{
		if(abs(x[i]) > max)
		{
			max = abs(x[i]);
			maxi = i;
		}
	}

	if(max > 1.0)
	{
//    divider("Normalize Do!");
		for(int i=0; i<nx; i++)
		{
			x[i] = x[i] / max;
		}
	}
}
