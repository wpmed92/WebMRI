/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <math.h>
#include "thresholding.h"

using namespace std;

void CreateHistogramImage(ZGrayByteImage& image, const float* hist, float ratio)
{
	image.Create(256, 100);
	image.FillPixels(200);
	
	PBYTE	histimg = image.PPixel(0,99);

	int i;
	float max = 0;
	for(i=0; i<256; i++) if(max < hist[i]) max = hist[i];
	
	for(i=0; i<256; i++, histimg++)
	{
		PBYTE   ptr = histimg;
		int height = int(hist[i] * ratio / max * 99);

		if(height > 100) height = 100;
		for(int j=99; j>=100-height; j--, ptr-=256) *ptr = 0;
	}
}

void CreateHistogramImage(ZGrayByteImage& image, const int* hist, float ratio)
{
	float fhist[256];

	for(int i=0; i<256; i++) fhist[i] = float(hist[i]);
	CreateHistogramImage(image, fhist, ratio);
}

void CreateHistogramImage(ZGrayByteImage& image, const unsigned int* hist, float ratio)
{
	float fhist[256];

	for(int i=0; i<256; i++) fhist[i] = float(hist[i]);
	CreateHistogramImage(image, fhist, ratio);
}

void CreateHistogramImage(ZGrayByteImage& image, const unsigned long* hist, float ratio)
{
	float fhist[256];

	for(int i=0; i<256; i++) fhist[i] = float(hist[i]);
	CreateHistogramImage(image, fhist, ratio);
}

int VarianceBilevelThresholding(const UINT* hist)
{
	if(!hist)
	{
		cout << "hist array is empty! 128 is returned!" << endl;
		return 128;
	}

	float fhist[256]; 
	int i, max, min;

	float nor = 0;
	for(i=0; i<256; i++) nor += hist[i];
	for(i=0; i<256; i++) fhist[i] = float(hist[i]) / nor;

	for(i=255; fhist[i]==0 && i>0; i--);
	max = i;

	for(i=0; fhist[i]==0 && i<max; i++);
	min = i;

	int l;
	float wk = 0, uk = 0, ut = 0;
	for(l=min; l<max; l++) ut += l * fhist[l];

	float thresvar=0;
	int thres = 128;
	int decrease = 0;
	for(l=min; l<max; l++)
	{
		wk += fhist[l];
		uk += l * fhist[l];

		if(wk * (1 - wk) < VERYTINY) continue;
		float var = ut * wk - uk;
		var = var * var / (wk * (1 - wk));
		if(thresvar < var)
		{
			thresvar = var;
			thres = l;
			decrease = 0;
		}
		else if(thresvar > var + 0.01) decrease ++;
		if(decrease == 10) break;
	}

	return thres;
}

bool VarianceThreelevelThresholding(const UINT* hist, int& thres1, int& thres2)
{
	if(!hist)
	{
		cout << "hist array is empty! 128 is returned!" << endl;
		return false;
	}

	float fhist[256], nor = 0; int i;
	for(i=0; i<256; i++) nor += hist[i];
	for(i=0; i<256; i++) fhist[i] = float(hist[i]) / nor;

	float ut = 0;
	for(int l=0; l<256; l++) ut += l * fhist[l];

	thres1 = 80, thres2 = 170;
	float thresvar=0, w2 = 0, u2 = 0;
	int l1;
	for(l1=10; l1<256; l1++)
	{
		w2 += fhist[l1];
		u2 += l1 * fhist[l1] ;
	}
	
	for(l1=10; l1<256; l1++)
	{
		if(w2 < VERYTINY) break;

		float w0 = 0, w1 = 0, u0 = 0, u1 = 0;
		int l2;
		for(l2=0; l2<l1; l2++)
		{
			w1 += fhist[l2];
			u1 += l2 * fhist[l2];
		}

		for(l2=0; l2<l1; l2++)
		{
			w0 += fhist[l2];
			u0 += l2 * fhist[l2];

			w1 -= hist[l2];
			u1 -= l2 * hist[l2];
	
			if(w0 * w1 < VERYTINY) continue;
			float uu0 = u0 * u0 / w0, uu1 = u1 * u1 / w1, uu2 = u2 * u2 / w2;
			float var = uu0 + uu1 + uu2 - ut * ut;
			if(thresvar < var)
			{
				thresvar = var;
				thres1 = l2;
				thres2 = l1;
			}
		}

		w2 -= fhist[l1];
		u2 -= l1 * fhist[l1];
	}

	return true;
}

bool VarianceMultilevelThresholding(const UINT* hist, int nblevel, int* thres)
{
	if(nblevel > 5) return false;
	if(!hist)
	{
		cout << "hist array is empty! 128 is returned!" << endl;
		return false;
	}

	int k[10]; 

	int i, max, min;
	float fhist[256], nor = 0;
	for(i=0; i<256; i++) nor += hist[i];
	for(i=0; i<256; i++) fhist[i] = float(hist[i]) / nor;

	for(i=255; fhist[i]==0 && i>0; i--);
	max = i;

	for(i=0; fhist[i]==0 && i<max; i++);
	min = i;

	for(i=0; i<=nblevel+1; i++) k[i] = i * (max-min) / (nblevel+1) + min;

	float e[10], u[10];

	bool flag = true;
	for(int iters=0; iters<1000 && flag; iters++)
	{
		int i;
		for(i=0; i<=nblevel; i++)
		{
			float nor = u[i] = 0;
			for(int j=k[i]; j<k[i+1]; j++)
			{
				u[i] += fhist[j] * j;
				nor += fhist[j];
			}
			u[i] /= nor;
		}

		flag = false;
		for(i=1; i<=nblevel; i++)
		{
			e[i] = (u[i-1] + u[i]) / 2 - k[i];
			e[i] += e[i] > 0 ? 0.5f : -0.5f;
			if(int(e[i]) != 0) flag = true;
		}

		for(i=1; i<=nblevel; i++)
			k[i] += int(e[i]);
	}
	if(flag) return false;

	for(i=0; i<nblevel; i++) thres[i] = k[i+1];
	
	return true;
}

void Threshold(ZImageBase& image, float thres1, float thres2, float thres3)
{
	switch(image.PixFmt())
	{
	case epixfmtGrayByte:
		Threshold((ZGrayByteImage&)image, thres1, thres2, thres3);
		break;
	case epixfmtGrayShort:
		Threshold((ZGrayShortImage&)image, thres1, thres2, thres3);
		break;
	case epixfmtGrayInt:
		Threshold((ZGrayIntImage&)image, thres1, thres2, thres3);
		break;
	case epixfmtGrayFloat:
		Threshold((ZGrayFloatImage&)image, thres1, thres2, thres3);
		break;
	case epixfmtRGBByte:
		Threshold((ZRGBByteImage&)image, thres1, thres2, thres3);
		break;
	case epixfmtRGBShort:
		Threshold((ZRGBShortImage&)image, thres1, thres2, thres3);
		break;
	case epixfmtRGBFloat:
		Threshold((ZRGBFloatImage&)image, thres1, thres2, thres3);
		break;
	default:
		return;
	}

	return;
}

int FindHistogramMode(const UINT* hist, bool& max)
{
	UINT hhh[256], mm=0;
	memcpy(hhh, hist, 256*sizeof(int));
	SmoothHistogram(hhh, 7);

	int lev = 0, i;
	for(i=10; i<246; i++)
		if(hhh[i] >= mm) mm = hhh[i], lev = i;

	if(max) return lev;		//return lev of maximum	if mode is 0

	max = true;
	UINT hist1[256];
	memcpy(hist1, hhh, 256*sizeof(int));

	hist1[lev] = mm * 2 / 3 + 1;

	for(i=lev+1; i<256; i++) 
	{
		if(hist1[i] > mm * 2 / 3 && hhh[i] >= hhh[i+1]) hist1[i] = 0;
		else break;
	}
	for(i=lev-1; i>20; i--)
	{
		if(hist1[i] > mm * 2 / 3 && hhh[i] >= hhh[i-1]) hist1[i] = 0;
		else break;
	}
	
	mm = 0;
	int ll = lev;
	for(i=50; i<246; i++)
		if(hist1[i] >= mm) mm = hist1[i], ll = i;
	
	if(ll > 120 && lev != ll)
	{
		int min = Min(ll, lev), max = Max(ll, lev);
		mm = 100000; 
		for(i=min; i<=max; i++)
			if(hhh[i] <= mm) mm = hhh[i], lev = i;

		max = false;
	}
	return lev;
}
