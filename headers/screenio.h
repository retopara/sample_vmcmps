#ifndef SCREEN_H
#define SCREEN_H

#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

void divider(string name);
void enter(int count, int interval);
void tensoroutput(mat *tensor, const int dim);
void tensoroutputfirst(mat *tensor, const int dim);
void screenesiteerr(
		mat *tensor,
		const int dim, 
		double *esitesave, 
		int **spinsave, 
		const int Nsample, 
		const int L, 
		const double etot, 
		const double thr);
#endif
