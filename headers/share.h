#ifndef SHARE_H 
#define SHARE_H 

//variables
extern int RECURSION_FLAG;
//variables end

void recflag_send(double *fvec_min, const double value);
void recflag_recv(int *Nsample, const int Nsample_inc);

#endif
