#ifndef FUNC
#define FUNC

#define DEBUG3

//function declarations
extern void REMC(double *fvec_min, double *r2_min, const double *tol, const double *betamin, 
                 const double *betamax, const int *n, const int *nswp, const int *nexch, 
                 const int *nrec, const int *nequi1, const int *nequi2, double *x_min, 
                 double (*target_func)(double *x_new), 
                 void (*x_maker)(const double *PI, const int *n, const int *i_x, double *x, double *x_new, 
                 double *sigma, int *csa_acpt_count, int *csa_count, int *seed_cau));

extern void CSA(const int *n, const double *PI, double *x_new, double *x,
                 double *x_min, double *sigma, int *csa_acpt_count, double *fvec_new, 
                 double *fvec, double *fvec_min, 
                 int *csa_count, const double *beta_myid, int *seed_met, int *seed_cau, 
                 double (*target_func)(double *x_new), 
                 void (*x_maker)(const double *PI, const int *n, const int *i_x, double *x, double *x_new, 
                                 double *sigma, int *csa_acpt_count, int *csa_count, int *seed_cau));

extern void AdjustTemper(double *beta, int *exch_acpt_count, const int *irec, const int *nrec);

extern void ExchTemper(const int *nscheme, int *exch_count, int *exch_acpt_count, int **exbeta, 
                       int *pro_index, int *beta_index, double *fvec, double *fvec_recv, 
                       double *beta);
//function declarations end

//variables
extern int PARASIZE;
extern int myid;
extern double PI;
//variables end

#endif
