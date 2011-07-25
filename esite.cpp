#include <iostream>
#include "armadillo"
#include "target.h"

using namespace arma;
using namespace std;

double esite(mat *left, mat* right, mat *tensor, int *spin, const int L, const double h)
{
	//calculating energy of current configuration using right and left 
	//calculate the normalized sigmaZ*sigmaZ energy
	double delta = 0.5, scale = 4.0;
	double ezz = scale * (spin[0] - delta) * (spin[L-1] - delta); //scale normalize spin to 1 & -1
	for(int i=0; i<L-1; i++)
	{
		ezz = ezz + scale * (spin[i] - delta) * (spin[i+1] - delta);
	}
	//      cout << "ezz=" << ezz << endl;

	//calculate sigmaX energy
	double prob = trace(right[0]);
	double probflip;
	double ex = 0;
	for(int Jsite=0; Jsite<L; Jsite++)
	{
		int f = flip(spin[Jsite]);
		if(Jsite == 0)
		{
			probflip = trace(tensor[f] * right[Jsite+1]);
		}
		else if(Jsite == L-1)
		{
			probflip = trace(left[Jsite-1] * tensor[f]);
		}
		else
		{
			probflip = trace(left[Jsite-1] * tensor[f] * right[Jsite+1]);
		}
		ex = ex + probflip;
	}
	ex = ex/prob;
	//      cout << "ex=" << ex << endl;
	double esite = -(ezz + h*ex)/L;

	//DISPLAY
	//      cout << fixed << setprecision (6) << "esite[" << Isite << "] = " << esite << "\t";
	//      enter(Isite, 4);
	//DISPLAY

	return esite;
}
