/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions of general image processing.

#ifndef __IMALGORITHM__H_
#define __IMALGORITHM__H_

#include <vector>
#include <algorithm>

#include "imagecore.h"
#include "zmath.h"
#include "filter.h"

void AutoRange(const ZImageBase& image, int& left, int& right, int& top, int& bottom, int& front, int& back);
bool AutoCrop(ZImageBase& image);
bool AutoRange(const std::vector<ZImageBase*>& images, int& left, int& right, int& top, int& bottom, int& front, int& back);
bool AutoCrop(std::vector<ZImageBase*>& images);

void FillPattern(ZGrayByteImage& image, UINT size = 20);

template <class T>
void RemoveOutlierT(ZImage<T>& image, float low, float high, bool keep=false)
{
	for(typename ZImage<T>::iterator p(image); p.more(); p++)
	{
		if(*p < low) *p = 0;
		if(*p > high) *p = keep ? T(high) : T(0);
	}
}

void RemoveOutlier(ZImageBase& image, float low, float high, bool keep=false);
void RemovePercent(ZImageBase& image, float low, float high, bool keep=false);

#endif
