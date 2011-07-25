/***********************************************************************************************
 *   FILE NAME:  fileio.cpp
 *		FUNCTION:		double nxinput();
 *								void xinput(const int nx, double *x);
 *								void tpoutput();
 *
 *		PURPOSE:		double nxinput(); 
 *								read target.input and calculate nx, the number of variables
 *
 *								void xinput(const int nx, double *x);
 *								read in given data to create x[]
 *
 *								void tpoutput();
 *								read in target.input to output parameters for target function.
 *
 ************************************************************************************************/
#include "target.h"
#include "toolbox.h"
#include "remc.h"
#include "screenio.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "armadillo"

using namespace arma;
using namespace std;

//input nx from file
double nxinput()
{
	int	dim, //spin only have up or down, spin 1/2 system, dim = 2
			schnum, //schmidt number of decomposition is 8
			nx; //the number of variables for target function

	ifstream ifs;
	string str;

	ifs.open("./input/target.input");
	if(!ifs)
	{
		cerr << "\nCannot find the file: ./input/target.input" << endl;
		abort();
	}

	ifs >> dim;
	getline(ifs, str);

	ifs >> schnum;
	getline(ifs, str);
	ifs.close();

	nx = dim * schnum * schnum;

	return nx;
}

//input x from file
void xinput(const int nx, double *x)
{
	ifstream ifs;
	string str;

	ifs.open("./input/x.input");
	if(!ifs)
	{
		cerr << "\nCannot find the file: ./input/x.input" << endl;
		abort();
	}
	for(int i=0; i<nx; i++)
	{
		ifs >> x[i];
		getline(ifs, str);
	}
	ifs.close();
}

//output to screen parameters of target function 
void tpoutput()
{
	int	dim, //spin only have up or down, spin 1/2 system, dim = 2
			schnum, //schmidt number of decomposition is 8
			N, //the length of spin chain is 2^N
			Nrelax, //the number of relaxation sweeps
			Nsample, //the number of sweeps when calculating energy.
			L, //the length of spin chain
			nx; //the number of variables
	double h; //the strength of transverse field

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
	ifs.close();
	//end read in

//  L = int(pow(2.0, N)); //length of spin chain
	L = N; //length of spin chain
	nx = dim * schnum * schnum; //the number of variables

	cout << "Parameters for Target Function:" << endl;
	cout << "dim=" << dim << endl;
	cout << "schnum=" << schnum << endl;
	cout << "N=" << N << endl;
	cout << "Nrelax=" << Nrelax << endl;
	cout << "Nsample=" << Nsample << endl;
	cout << "L=" << L << endl;
	cout << "nx=" << nx << endl;
}

//output sampling etot values
void etotoutput(const double etot)
{
	static int count = 0;
	count ++;
	ofstream ofs;
	if(count == 1)
	{
		ofs.open("./data/etot.data");
		if(!ofs)
		{
			cerr << "\ncan't open ./data/etot.data!" << endl;
			abort();
		}
		ofs << "etot values" << endl << endl;
	}
	else
	{
		ofs.open("./data/etot.data", ios_base::app);
		if(!ofs)
		{
			cerr << "\ncan't open ./data/etot.data!" << endl;
			abort();
		}
	}
	ofs << fixed << setprecision(8) << count << ", " << etot << endl;
	ofs.close();
}

//output error sampling information to different files(according to myid) when etot is too low
void esiteerr(
		mat *tensor,
		const int dim, 
		double *esitesave, 
		int **spinsave, 
		const int Nsample, 
		const int L, 
		const double etot, 
		const double thr)
{
	static int count = 0;
	count ++;

	ostringstream oss;
	oss << "./data/esiteerr" << myid << ".data";
	cout << oss.str() << endl;
	
	ofstream ofs;
	if(count == 1)
	{
		ofs.open(oss.str().c_str());
		if(!ofs)
		{
			cerr << "\ncan't open esiteerr data file!" << endl;
			abort();
		}
		ofs << "myid = " << myid << endl << endl;
	}
	else
	{
		ofs.open(oss.str().c_str(), ios_base::app);
		if(!ofs)
		{
			cerr << "\ncan't open esiteerr data file!" << endl;
			abort();
		}
	}

	ofs << "\nCOUNT = " << count << " ===============================================================" << endl;
	ofs << "etot < " << thr << ". ERROR! Total energy is :" << endl;
	ofs << "etot = " << etot << endl;

	ofs << "\ntensors are:" << endl;
	
	for(int i=0; i<dim; i++)
	{
		ofs << tensor[i] << endl;
	}
	
	ofs << "\nesitesave and spinsave : " << endl;
	for(int Isample=0; Isample<Nsample; Isample++)
	{
		for(int Isite=0; Isite<L; Isite++)
		{
			ofs << "\nIsample=" << Isample;
			ofs << "\tIsite=" << Isite << endl;

			ofs << "esite=" << esitesave[Isample*L+Isite] << "\t";
			for(int i=0; i<L; i++)
			{
				ofs << spinsave[Isample*L+Isite][i] << " ";
			}
			ofs << endl;
		}
	}
}
