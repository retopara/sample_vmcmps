#include <iostream>
#include <string>
#include "armadillo"
#include "screenio.h"
#include "remc.h"

using namespace arma;
using namespace std;

//printer a divider for output
void divider(string name)
{
	int MAX = 100;
	cout << endl;
	for(int i=0; i<MAX; i++)
	{
		cout << "=";
	}
	cout << endl << name << " : " << endl;
	for(int i=0; i<MAX; i++)
	{
		cout << "=";
	}
	cout << endl << endl;
}

//print data with a tabbed format
void tab(double data, int count, int interval)
{
}

//print return "\n" for a certain interval of count
void enter(int count, int interval)
{
	count++;
	if(count%interval == 0)
	{
		cout << endl;
	}
}

//output tensors
void tensoroutputfirst(mat *tensor, const int dim)
{
	static int count = 0;
	count ++;
	if(count == 1)
	{
		divider("tensors are");
		for(int i=0; i<dim; i++)
		{
			cout << tensor[i] << endl;
		}
	}
}

void tensoroutput(mat *tensor, const int dim)
{
	divider("tensors are");
	for(int i=0; i<dim; i++)
	{
		cout << tensor[i] << endl;
	}
}

//readin input file of certain format
void readin()
{}


void screenesiteerr(
		mat *tensor,
		const int dim, 
		double *esitesave, 
		int **spinsave, 
		const int Nsample, 
		const int L, 
		const double etot, 
		const double thr)
{
	if(myid == 0)
	{
		static int count = 0;
		count ++;

		cout << "\nCOUNT = " << count << " ===============================================================" << endl;
		cout << "myid = " << myid << endl;
		cout << "etot < " << thr << ". ERROR! Total energy is :" << endl;
		cout << "etot = " << etot << endl;

		cout << "\ntensors are:" << endl;
		tensoroutput(tensor, dim);

		cout << "\nesitesave and spinsave : " << endl;
		for(int Isample=0; Isample<Nsample; Isample++)
		{
			for(int Isite=0; Isite<L; Isite++)
			{
				cout << "\nIsample=" << Isample;
				cout << "\tIsite=" << Isite << endl;

				cout << "esite=" << esitesave[Isample*L+Isite] << "\t";
				for(int i=0; i<L; i++)
				{
					cout << spinsave[Isample*L+Isite][i] << " ";
				}
				cout << endl;
			}
		}
	}
}
