/***********************************************************************************************
 *	 AUTHOR:		 Zhen Wang(USTC LQCC)
 *							 wangzhen.ustc@gmail.com
 *	 
 *   FILE NAME:  share.cpp
 *   FUNCTION:   void recflag_send()
 *							 void recflag_recv();
 *
 *   PURPOSE:    Though built separated, the REMC part and the target function part of the code 
 *							 sometimes needs to communicate with each other. In order to do this, we set 
 *							 flags here to carry these info, by setting these flags in REMC part and receving
 *							 in target funcion part, the communication is done.
 *
 ***********************************************************************************************/

int RECURSION_FLAG(0);

#include <iostream>

using namespace std;

//IN REMC, set the flag when there is need to send info
void recflag_send(double *fvec_min, const double value)
{
	RECURSION_FLAG = 1;
	*fvec_min = value;
//  cout << "\nRECURSION_FLAG = " << RECURSION_FLAG << endl;
}

//IN TARGET FUNCTION, receive the flag and do certain actions, use when flag == 1
void recflag_recv(int *Nsample, const int Nsample_inc)
{
	*Nsample += Nsample_inc;
	RECURSION_FLAG = 0;
//  cout << "\nRECURSION_FLAG = " << RECURSION_FLAG << endl;
}

