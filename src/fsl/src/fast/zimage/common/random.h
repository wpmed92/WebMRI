/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for generate random numbers.

#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <cstdlib>
#include <cmath>

inline float rrand(float mi, float ma)
{
	return float(rand()) / (RAND_MAX+1) * (ma - mi) + mi;
}

float rand1(long & idum);
float rand2(long & idum);

// To initialize normal random generator
// Call with idum a negative integer to initialize
void set_noraml_rand_seed(long seed);

// normal random variate generator
// return only kernel; needed to be multiplied by sigma and then added by the mean
// m + k * s ~ N(m, s)
double normal_random_kernel();

// normal random variate generator with mean and sigma
inline double normal_random(double mean, double sigma) { return mean + normal_random_kernel() * sigma; }


inline double rayleigh_random(double sigma)
{
	double x1 = normal_random(0, sigma);
	double x2 = normal_random(0, sigma);
	return sqrt(x1*x1 + x2*x2) ;
}

#endif
