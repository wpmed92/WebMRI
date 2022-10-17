/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef __MATRIXD_H_
#define __MATRIXD_H_

#include "mathcommon.h"

class ZMatrixD
{
public:
	typedef float* iterator;
	typedef const float* const_iterator;

	virtual iterator begin() = 0;
	virtual const_iterator begin() const = 0;
	virtual iterator end() = 0;
	virtual const_iterator end() const = 0;
	
	virtual float& operator () (int r, int c) = 0;
	virtual const float& operator() (int r, int c) const  = 0;

	virtual float operator~ () const = 0;
	ZMatrix2D			Inverse (void) const
	{
		float val = ~(*this);
		if ( ! val ) return ZMatrix2D();

		ZMatrix2D temp;
		temp.d11 = d22 / val;
		temp.d22 = d11 / val;
		temp.d12 = - d12 / val;
		temp.d21 = - d21 / val;
		
		return temp;
	}

	ZMatrix2D			Transpose (void) const
	{
		ZMatrix2D		mat;
		mat.d11 = d11; mat.d12 = d21; mat.d21 = d12; mat.d22 = d22;
		return mat;
	}
};

#endif
