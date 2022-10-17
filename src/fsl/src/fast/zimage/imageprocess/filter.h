/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions of image filtering.

#ifndef __IMFILTER__H_
#define __IMFILTER__H_

#include "imagecore.h"
#include "zmath.h"

void GaussianBlur(ZImageBase& image, float sigma, bool b2d=false);
void GaussianBlur(ZGrayFloatImage& image, float x_sigma, float y_sigma, float z_sigma, bool b2d=false);

ZImageBase*	Gradient(const ZImageBase& image);
ZImageBase*	XGradient(const ZImageBase& image);
ZImageBase*	YGradient(const ZImageBase& image);
ZImageBase*	Gradient3D(const ZImageBase& image);

#endif
