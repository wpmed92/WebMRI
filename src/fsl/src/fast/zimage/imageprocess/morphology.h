/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for morphology operations.

#ifndef __MORPHOLOGY__H_
#define __MORPHOLOGY__H_

#include <vector>
#include <algorithm>

#include "imagecore.h"
#include "zmath.h"

void CreateSquareKernel(std::vector< ZPoint<int> >& kernel, int width, int height, int depth=1);
void CreateRoundKernel(std::vector< ZPoint<int> >& kernel, int radius, bool b3D=false);
void Dilation(ZImageBase& image, const std::vector< ZPoint<int> >& kernel);
void Erosion(ZImageBase& image, const std::vector< ZPoint<int> >& kernel);
void Opening(ZImageBase& image, const std::vector< ZPoint<int> >& kernel);
void Closing(ZImageBase& image, const std::vector< ZPoint<int> >& kernel);

#endif
