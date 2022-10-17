/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#ifdef _MSC_VER
#pragma warning( disable : 4514 )  
#endif

#include "commfunc.h"

using namespace std;

void Reverse2(void *ptr, int num)
{
	char *cptr = (char *) ptr;
	for (int i=0; i<num; i++) 
	{
		char tmp = *cptr;
		*cptr = *(cptr + 1);
		*(cptr + 1) = tmp;
		cptr += 2;
	}
}

void Reverse4(void *ptr, int num)
{
	char *cptr = (char *) ptr;
	for (int i=0; i<num; i++) 
	{
		char tmp1 = *cptr;
		char tmp2 = *(cptr + 1);
		*cptr = *(cptr + 3);
		*(cptr + 1)= *(cptr + 2);
		*(cptr + 2)= tmp2;
		*(cptr + 3)= tmp1;
		cptr += 4;
	}
}

void Reverse8(void *ptr, int num)
{
	char *cptr = (char *) ptr;
	for (int i=0; i<num; i++) 
	{
		char tmp1 = *cptr;
		char tmp2 = *(cptr + 1);
		char tmp3 = *(cptr + 2);
		char tmp4 = *(cptr + 3);
		*cptr = *(cptr + 7);
		*(cptr + 1)= *(cptr + 6);
		*(cptr + 2)= *(cptr + 5);
		*(cptr + 3)= *(cptr + 4);
		*(cptr + 4)= tmp4;
		*(cptr + 5)= tmp3;
		*(cptr + 6)= tmp2;
		*(cptr + 7)= tmp1;
		cptr += 8;
	}
}
