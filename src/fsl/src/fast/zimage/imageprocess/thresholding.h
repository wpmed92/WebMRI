/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for graylevel thresholding.

#ifndef __THRESHOLDING__H_
#define __THRESHOLDING__H_

#include <vector>

#include "imagecore.h"
#include "zmath.h"

void CreateHistogramImage(ZGrayByteImage& image, const float* hist, float ratio = 1);
void CreateHistogramImage(ZGrayByteImage& image, const int* hist, float ratio = 1);
void CreateHistogramImage(ZGrayByteImage& image, const unsigned int* hist, float ratio = 1);
void CreateHistogramImage(ZGrayByteImage& image, const unsigned long* hist, float ratio = 1);

int  VarianceBilevelThresholding(const UINT* hist);
bool VarianceThreelevelThresholding(const UINT* hist, int& thres1, int& thres2);
bool VarianceMultilevelThresholding(const UINT* hist, int nblevel, int* thres);

void Threshold(ZImageBase& image, float thres1, float thres2=-1, float thres3=-1);

int FindHistogramMode(const UINT* hist, bool& max);

template <class TPixel>
void Threshold(ZImage<TPixel>& image, float thres1, float thres2=-1, float thres3=-1)
{
	typename ZImage<TPixel>::iterator pixel(image);
	for(; pixel.more(); pixel++)
	{
		if(*pixel < thres1) *pixel = 0;
		else if(thres2 > thres1)
		{
			if(*pixel < thres2) *pixel = 1;
			else if(thres3 < thres2) 
				if(*pixel < thres3) *pixel = 2;
				else *pixel = 3;
		}
		else *pixel = 1;
	}
}

template <class T> 
void SmoothHistogram(T* hist, int nIter=1)
{
	float kernel[9], sum=0;

	for(int i=-4; i<=4; i++)
	{
		kernel[i+4] = float(std::pow(float(2.71), -float(i*i)/32));
		sum += kernel[i+4];
	}

	T hh[256];
	for(int iter=0; iter<nIter; iter++)
	{
		memcpy(hh, hist, 256*sizeof(T));

		for(int i=10; i<=245; i++)
		{
			float value = 0;
			for(int j=0; j<9; j++) value += kernel[j] * hh[i-4+j];
			
			hist[i] = T(value / sum);
		}
	}
}

/*
error: cannot convert 'float*' to 'const UINT*' for argument '1' to 
'int VarianceBilevelThresholding(const UINT*)'template <class TPixel>

float VarianceBilevelThresholding(ZImage<TPixel>& image)
{
	float		hist[256];
	memset(hist, 0, 256 * sizeof(float));

	float min=0, max=0; image.MinMax(min, max);
	double	nor = 255.0 / (max - min);

	typename ZImage<TPixel>::iterator pixel(image);
	for(; pixel.more(); pixel.Next()) hist[int((*pixel - min) * nor)] ++;

	for(int i=0; i<256; i++) hist[i] /= image.Size();

	float thres = float(VarianceBilevelThresholding(hist));

	thres = float(thres / nor + min);

	return thres;
}
*/

/*
error: cannot convert 'float*' to 'const UINT*' for argument '1' 
to 'bool VarianceThreelevelThresholding(const UINT*, int&, int&)'

template <class TPixel>
void VarianceThreelevelThresholding(ZImage<TPixel>& image, float& thres1, float& thres2)
{
	float		hist[256];
	memset(hist, 0, 256 * sizeof(float));

	float min=0, max=0; image.MinMax(min, max);
	double	nor = 255.0 / (max - min);

	typename ZImage<TPixel>::iterator pixel(image);
	for(; pixel.more(); pixel.Next()) hist[int((*pixel - min) * nor)] ++;

	for(int i=0; i<256; i++) hist[i] /= image.Size();

	thres1 = 80, thres2 = 170;

	VarianceThreelevelThresholding(hist, thres1, thres2);

	thres1 = float(thres1 / nor + min);
	thres2 = float(thres2 / nor + min);
}
*/
#endif
