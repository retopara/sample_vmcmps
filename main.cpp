/*************************************************************
 *	 AUTHOR:		 Zhen Wang(USTC LQCC)
 *							 wangzhen.ustc@gmail.com
 *
 *   FILE NAME:  Main.cpp
 *   FUNCTION:   Main
 *
 *   PURPOSE:    The main function
 *************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include "target.h"
#include "remc.h"
#include "fileio.h"
#include "screenio.h"

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

//global variables
int myid;
int PARASIZE;
double PI;

int main(int argc, char **argv)	//the int main() must follow ANSI C standard, in order to use MPI.int argc;
{
	time_t t1, t2;
	t1 = time(NULL);
	//Parallel part initializing

	MPI_Init(&argc, &argv);
#ifdef DEBUG
	printf("MPI_Init complete!\n");
#endif

	MPI_Comm_size(MPI_COMM_WORLD, &PARASIZE);

#ifdef DEBUG
	printf("MPI_Comm_size complete!\n");
	printf("SIZE = %d\n", PARASIZE);
#endif

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#ifdef DEBUG
	printf("MPI_Comm_rank complete!\n");
	printf("myid = %d\n", myid);
#endif
	if(myid == 0)
	{
		divider("Program Start");
	}

	//file input
	double bmin, bmax, yb, r2, tol, *x_min;
	int nx, nswp, nexch, nrec, nequi1, nequi2;
	int NP, NC, NPAR, NONE, NLEV;

	string str;
	ifstream ifs("./input/main.input");
	if(!ifs)
	{
		cerr << "\nCannot find the file: ./input/main.input" << endl;
		abort();
	}

	ifs >> nrec;
	getline(ifs, str);

	ifs >> nexch;
	getline(ifs, str);

	ifs >> nswp;
	getline(ifs, str);

	ifs >> nequi1;
	getline(ifs, str);

	ifs >> nequi2;
	getline(ifs, str);

	ifs >> bmin;
	getline(ifs, str);

	ifs >> bmax;
	getline(ifs, str);

	ifs >> tol;
	getline(ifs, str);
	ifs.close();

	nx = nxinput(); //calculate the number of variables in x[]

	if(myid == 0)
	{
		cout << "\nThe number of cores: "; 
		cout << "\nPARASIZE=" << PARASIZE << endl << endl;
		tpoutput();
		cout << "\nParameters for REMC: "; 
		cout << "\nnx=" << nx;
		cout << "\nnrec=" <<nrec;
		cout << "\nnexch=" <<nexch;
		cout << "\nnswp=" <<nswp;
		cout << "\nnequi1=" <<nequi1;
		cout << "\nnequi2=" <<nequi2;
		cout << "\nbmin=" <<bmin;
		cout << "\nbmax=" <<bmax;
		cout << "\ntol=" <<tol << endl;
	}

	//start REMC
	REMC(&yb, &r2, &tol, &bmin, 
			&bmax, &nx, &nswp, &nexch, 
			&nrec, &nequi1, &nequi2, x_min, 
			Energy, XGenerator);

	//output running time
	t2 = time(NULL);

	if(myid == 0){
		divider("Program Running Time");
		cout << "myid = 0, running time is:";
		printf("\nprogram run time(hour) = %f", (t2-t1)/3600.0);
		printf("\nprogram run time(minute) = %f", (t2-t1)/60.0);
		printf("\nprogram run time(second) = %d\n", (t2-t1));
	}

	//finalize
	MPI_Finalize();
}
