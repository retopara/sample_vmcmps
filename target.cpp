#include <iostream>
#include <iomanip>

#include "armadillo"
#include "target.h" 
#include "screenio.h"
#include "fileio.h"
#include "share.h"
#include "remc.h"

using namespace arma;
using namespace std;

double Energy(double *x_new)
{
	//define global constant
	static int count = 0;
	count ++;
	static int	
		dim, //spin only have up or down, spin 1/2 system, dim = 2
		schnum, //schmidt number of decomposition is 8
		N, //the length of spin chain is 2^N
		Nrelax, //the number of relaxation sweeps
		Nsample, //the number of sweeps when calculating energy.
		L, //the length of spin chain
		Nsample_inc, //the increment of Nsample per recursion.
		nx; //the number of variables
	static double h; //the strength of transverse field

	bool flipflag; //flipflag is true when a flip is accepted
	double prob, prob0, prob1, //to calculate metropolis acceptance rate
				 etot; //total energy

	if(count == 1)
	{
		//read in input file
		string str;

		ifstream ifs("./input/target.input");
		if(!ifs)
		{
			cerr << "\nCannot find the file: ./input/target.input" << endl;
			abort();
		}

		ifs >> dim;
		getline(ifs, str);

		ifs >> schnum;
		getline(ifs, str);

		ifs >> N;
		getline(ifs, str);

		ifs >> Nrelax;
		getline(ifs, str);

		ifs >> Nsample;
		getline(ifs, str);

		ifs >> h;
		getline(ifs, str);

		ifs >> Nsample_inc;
		getline(ifs, str);
		ifs.close();
		//end read in

		//    L = int(pow(2.0, N)); //length of spin chain
		L = N; //length of spin chain
		nx = dim * schnum * schnum; //the number of variables
	}

	if(RECURSION_FLAG == 1)
	{
		recflag_recv(&Nsample, Nsample_inc);
		if(myid == 0) cout << "\nNsample= " << Nsample << endl;
	}

	//Normalize x_new[];
	Normalize(nx, x_new);

	//initialize tensors from x_new[]
	mat *tensor = new mat[dim];
	for(int i=0; i<dim; i++)
	{
		tensor[i] = mat(schnum, schnum);
		tensor[i].zeros();
	}

	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<schnum; j++)
			for(int k=0; k<schnum; k++)
			{
				int l = i*(schnum*schnum) + j*schnum + k;
				tensor[i](j, k) = x_new[l];
			}
	}

	//DISPLAY
//  tensoroutput(tensor, dim);

	//new spinsave and esitesave
	double *esitesave = new double[Nsample*L];
	int **spinsave = new int*[Nsample*L];
	for(int i=0; i<Nsample*L; i++)
	{
		spinsave[i] = new int[L];
	}

	//initialize the spin
	int *spin = new int[L];
	for(int i=0; i<L; i++)
	{
		double r = rand()/(RAND_MAX*1.0);
		if(r >= 0.5)
			spin[i] = 1;
		else
			spin[i] = 0;
	}

	//create left
	mat *left = new mat[L];
	for(int i=0; i<L; i++)
	{
		left[i] = mat(schnum, schnum);
		left[i].zeros();
	}

	//create right
	mat *right = new mat[L];
	for(int i=0; i<L; i++)
	{
		right[i] = mat(schnum, schnum);
		right[i].zeros();
	}

	//calculate all the right
	right[L-1] = tensor[spin[L-1]];
	for(int i=L-2; i>=0; i--)
	{
		right[i] = tensor[spin[i]] * right[i+1];
	}

	//relaxing the spin chain, doing metropolis for Nrelax times
	for(int i=0; i<Nrelax; i++)
	{
		//sweeping from left to right
		for(int Isite=0; Isite<L; Isite++)
		{
			flipflag = false;
			int f = flip(spin[Isite]); //flip current site for a trial

			if(Isite == 0)
			{
				prob0 = trace(right[Isite]);
				prob1 = trace(right[Isite+1] * tensor[f]);
			}
			else if(Isite == L-1)
			{
				prob0 = trace(left[Isite-1] * right[Isite]);
				prob1 = trace(left[Isite-1] * tensor[f]);
			}
			else
			{
				prob0 = trace(left[Isite-1] * right[Isite]);
				prob1 = trace(left[Isite-1] * tensor[f] * right[Isite+1]);
			}

			//judge if the flip should be accepted, prob0 is old energy, prob1 is new energy
			double prob = pow((prob1/prob0), 2.0);
			double r = rand()/(RAND_MAX*1.0);
			if(r < min(prob, 1.0))
			{
				spin[Isite] = f;
				flipflag = true;
			}

			//update left while using right
			if(Isite == 0)
			{
				left[Isite] = tensor[spin[Isite]];
			}
			else
			{
				left[Isite] = left[Isite-1] * tensor[spin[Isite]];
			}
		}//end sweeping from left to right

		//sweeping from right to left
		for(int Isite=L-1; Isite>=0; Isite--)
		{
			flipflag = false;
			int f = flip(spin[Isite]);
			if(Isite ==  L-1)
			{
				prob0 = trace(left[Isite]);
				prob1 = trace(left[Isite-1] * tensor[f]);
			}
			else if(Isite == 0)
			{
				prob0 = trace(left[Isite] * right[Isite+1]);
				prob1 = trace(tensor[f] * right[Isite+1]);
			}
			else
			{
				prob0 = trace(left[Isite] * right[Isite+1]);
				prob1 = trace(left[Isite-1] * tensor[f] * right[Isite+1]);
			}

			//judge if the flip should be accepted, prob0 is old energy, prob1 is new energy
			double r = rand()/(RAND_MAX*1.0);
			prob = pow((prob1/prob0), 2.0);
			if(r < min(prob, 1.0))
			{
				spin[Isite] = f;
				flipflag = true;
			}

			//updating right while using left
			if(Isite == L-1)
			{
				right[Isite] = tensor[spin[Isite]];
			}
			else
			{
				right[Isite] = tensor[spin[Isite]] * right[Isite+1];
			}
		}
	}//end relaxing spin chain

	//calculation of energy
	etot = 0;
	for(int Isample=0; Isample<Nsample; Isample++)
	{
		for(int Isite=0; Isite<L; Isite++)
		{
			flipflag = false;
			int f = flip(spin[Isite]); //flip current site for a trial
			if(Isite == 0)
			{
				prob0 = trace(right[Isite]);
				prob1 = trace(tensor[f] * right[Isite+1]);
			}
			else if(Isite == L-1)
			{
				prob0 = trace(left[Isite-1] * right[Isite]);
				prob1 = trace(left[Isite-1] * tensor[f]);
			}
			else
			{
				prob0 = trace(left[Isite-1] * right[Isite]);
				prob1 = trace(left[Isite-1] * tensor[f] * right[Isite+1]);
			}

			//judge if the flip should be accepted, prob0 is old energy, prob1 is new energy
			prob = pow((prob1/prob0), 2.0);
			double r = rand()/(RAND_MAX*1.0);
			if(r < min(prob, 1.0))
			{
				spin[Isite] = f;
				flipflag = true;
			}
			//calculate all the Left and Right for sampling energy
			if(Isite == 0)
			{
				left[Isite] = tensor[spin[Isite]];
				right[Isite] = tensor[spin[Isite]] * right[Isite+1];
				for(int Jsite=1; Jsite<L; Jsite++)
				{
					left[Jsite] = left[Jsite-1] * tensor[spin[Jsite]];
				}
			}
			else if(Isite == L-1)
			{
				left[Isite] = left[Isite-1] * tensor[spin[Isite]] ;
				right[Isite] = tensor[spin[Isite]];
				for(int Jsite=L-2; Jsite>=0; Jsite--)
				{
					right[Jsite] = tensor[spin[Jsite]] * right[Jsite+1];
				}
			}
			else
			{
				for(int Jsite=Isite; Jsite<L; Jsite++)
				{
					left[Jsite] = left[Jsite-1] * tensor[spin[Jsite]];
				}
				for(int Ksite=Isite; Ksite>=0; Ksite--)
				{
					right[Ksite] = tensor[spin[Ksite]] * right[Ksite+1];
				}
			}

			double es = esite(left, right, tensor, spin, L, h);
			etot = etot + es;

			//save current spin and esite
			esitesave[Isample*L+Isite] = es;
			for(int i=0; i<L; i++)
			{
				spinsave[Isample*L+Isite][i] = spin[i];
			}

			//      divider("esitesave and spinsave");
			//      cout << "esitesave[" << Isample*L+Isite << "]=" << esitesave[Isample*L+Isite] << endl;
			//      for(int i=0; i<L; i++)
			//      {
			//        cout << spinsave[Isample*L+Isite][i] << " ";
			//      }
			//      cout << endl;
		}//end of one sweep
	}//end of Nsample sweeps
	etot = etot/(Nsample * L * 1.0);

	//  double thr = -2.0;
	//  if(etot < thr)
	//  {
	//    screenesiteerr(tensor, dim, esitesave, spinsave, Nsample, L, etot, thr);
	//  }

	//delete! very important!
	delete [] tensor;
	delete [] esitesave;

	for(int i=0; i<Nsample*L; i++)
	{
		delete [] spinsave[i];
	}
	delete [] spinsave;

	delete [] spin;
	delete [] left;
	delete [] right;

	return etot;
}
