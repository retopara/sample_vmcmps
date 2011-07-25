#ifndef FILEIO_H
#define FILEIO_H 

#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

double nxinput();
void xinput(const int nx, double *x);
void tpoutput();
void etotoutput(const double etot);
void esiteerr(
		mat *tensor,
		const int dim, 
		double *esitesave, 
		int **spinsave, 
		const int Nsample, 
		const int L, 
		const double etot, 
		const double thr);
#endif
