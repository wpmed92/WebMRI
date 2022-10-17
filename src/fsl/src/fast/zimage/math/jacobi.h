/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for computing jacobi.

#ifndef __JACOBI_H_  // prevent multiple includes
#define __JACOBI_H_

bool _Jacobi(ZMatrix<float>& a, ZVector<float>& d, ZMatrix<float>& v);

template <class T> 
bool Jacobi(const ZMatrix<T>& m, ZVector<float>& d, ZMatrix<float>& v)
{
	ZMatrix<float> a = m;
	return _Jacobi(a, d, v);
}


#endif
