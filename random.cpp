/***************************************************************************************
*   FILE NAME:  Ran.cpp
*
*   PURPOSE:    "Minimal" random number generator of Park and Miller with
*               Bays-Durham shuffle and added safeguards. Returns a uniform
*               random deviate between 0.0 and 1.0 (exclusive of the end points).
*               Call with idum a negative integer to initialize; thereafter do
*               not alter idum between successive deviates in a sequence. RNMX
*               should approximate the largest doubleing point value that is 
*               less than 1.
*
*   INPUT:      the address of idum(the seed).
*   OUTPUT:     random number between 0.0 and 1.0, and change idum value in the process
****************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "toolbox.h"

extern double Ran(int *idum)
{
    //internal definitions
    static int IA, IM, IQ, IR, NTAB;
    static double EPS;

    static int i, j, ndiv;
    static double am, rnmx;

    static int ran_count = 0;
    static int iy = 0;
    static int *iv;
    //internal definitions end

    IA=16807;   IM=2147483647;
    IQ=127773;  IR=2836;    NTAB=32;
    EPS = 1.2f * pow(10.0f, -7.0f);

    am = 1.0f/IM;
    rnmx = 1.0f - EPS;
    ndiv = 1 + (IM-1)/NTAB;

    if(ran_count == 0)
    {
        //malloc
        iv = (int *)malloc(sizeof(int)*NTAB +1);
        if(iv == NULL) printf("iv malloc fail!\n"); //how to free????????????
        //malloc end

        //initializing
        for(i=0; i<NTAB; i++){
            iv[i] = 0;
        }
    }
    ran_count++;
    //initializing end

    //generate random number, and change the seed
    if((*idum <= 0) || (iy == 0))
    {
        if(-*idum < 1){
            *idum = 1;
        }
        else{
            *idum = -*idum;
        }
        for(i=NTAB+7; i>=0; i--)
        {
            j = *idum / IQ;
            *idum = IA * (*idum - j*IQ) - IR*j;
            if(*idum < 0){
                *idum += IM;
            }
            if(i < NTAB){
                iv[i] = *idum;
            }
        }
        iy = iv[0];
    }

    j = *idum / IQ;
    *idum = IA * (*idum - j*IQ) - IR*j;
    if(*idum < 0){
        *idum += IM;
    }
    i = iy/ndiv;
    iy = iv[i];
    iv[i] = *idum;
    if(am*iy < rnmx){
        return am*iy;
    }
    else{
        return rnmx;
    }
}

