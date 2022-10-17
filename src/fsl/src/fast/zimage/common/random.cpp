/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "random.h"

/****************************************************************
Minimal" random number generator of Park and Miller with Bays-Durham shu.e and added
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
successive deviates in a sequence. RNMX should approximate the largest oating value that is
less than 1. 
****************************************************************/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7f
#define RNMX (1.0f-EPS)

float rand1(long& idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (idum <= 0 || !iy) 
	{
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) 
		{
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=float(AM*iy)) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/***********************************************************************
Long period (> 2*10^18 ) random number generator of L'Ecuyer with Bays-Durham shu.e
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest oating
value that is less than 1.
**************************************************************************************/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7f
#define RNMX (1.0f-EPS)

float rand2(long& idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (idum <= 0) 
	{
		if (-idum < 1) idum=1;
		else idum = -idum;
		idum2=idum;
		for (j=NTAB+7;j>=0;j--) 
		{
			k=idum/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	if ((temp=float(AM*iy)) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

static long normal_random_seed = -10;
void set_noraml_rand_seed(long seed) { normal_random_seed = seed; }

// normal random variate generator
// return only kernel; needed to be multiplied by sigma and then added by the mean
// m + k * s ~ N(m, s)
double normal_random_kernel()
{
	double x1, x2, w, y1;

	static double y2;	
	static int use_last = 0;
	static long ini = normal_random_seed;

	if (use_last)		        // use value from previous call
	{		
		y1 = y2;
		use_last = 0;	
	}	
	else	
	{		
		do 
		{			
			x1 = 2.0 * rand1(ini) - 1.0;
			x2 = 2.0 * rand1(ini) - 1.0;			
			w = x1 * x1 + x2 * x2;		
		} 
		while ( w >= 1.0 );
		w = sqrt( (-2.0 * log( w ) ) / w );		
		y1 = x1 * w;		
		y2 = x2 * w;		
		use_last = 1;
	}	
	return y1;
}
