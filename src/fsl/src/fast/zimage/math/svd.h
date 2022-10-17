/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for SVD decomposition and matrix inverse.

#ifndef __SVD_H_  // prevent multiple includes
#define __SVD_H_

bool _inv(ZMatrix<double>& U);
bool _SVD(ZMatrix<double>& U, ZVector<double>& W, ZMatrix<double>& V);

template<class T> 
bool SVD(const ZMatrix<T>& m, ZMatrix<double>& U, ZVector<double>& W, ZMatrix<double>& V)
{
	if(m.NElements() == 0) 
	{
		ZError("SVD Decomposition", "Matrix is empty!");
		return false;
	}
	return _SVD(U=m, W, V);
}

template<class T> ZMatrix<double> inv(const ZMatrix<T>& m)
{
	if(m.NElements() == 0) 
	{
		ZError("Matrix Inverse", "Matrix is empty!");
		return false;
	}
	ZMatrix<double> U = m;
        if (_inv(U) == false) {
	  ZError("Matrix Inverse", "Matrix non invertible!");
	}
	return U;
}

#endif





