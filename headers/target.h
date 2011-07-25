#ifndef TARGET_H
#define TARGET_H

#include "armadillo"
using namespace arma;

int flip(int spin);
double Energy(double *x_new);
double esite(mat *left, mat* right, mat *tensor, int *spin, const int L, const double h);
void Normalize(const int nx, double *x);
extern void XGenerator(const double *PI, const int *n, const int *i_x, double *x, double *x_new, 
                       double *sigma, int *csa_acpt_count, int *csa_count, int *seed_cau);

#endif
