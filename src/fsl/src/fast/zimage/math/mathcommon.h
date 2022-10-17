/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares comon functions for math classes

#ifndef __MATHCOMMON_H_  // prevent multiple includes
#define __MATHCOMMON_H_

#include "matrix.h" 

std::ostream& operator << (std::ostream& s, const ZMatrix<int>& m);
std::ostream& operator << (std::ostream& s, const ZMatrix<float>& m);
std::ostream& operator << (std::ostream& s, const ZMatrix<double>& m);
std::ostream& operator << (std::ostream& s, const ZMatrix2D& m);
std::ostream& operator << (std::ostream& s, const ZMatrix3D& m);
std::ostream& operator << (std::ostream& s, const ZMatrix4D& m);

ZVector<float> RandVector(UINT r);
ZMatrix<float> RandMatrix(UINT r, UINT c);

bool svdcmp(ZMatrix<double>& U, ZVector<double>& W, ZMatrix<double>& V);
bool LUdcmp(ZMatrix<double>& A);
bool LUdcmp(ZMatrix<double>& A, ZVector<int>& indx);
bool LUdcmp(ZMatrix<float>& A);
bool LUdcmp(ZMatrix<float>& A, ZVector<int>& indx);
void LUbksb (const ZMatrix<double>& A, const ZVector<int>& indx, ZVector<double>& B );
void LUbksb (const ZMatrix<float>& A, const ZVector<int>& indx, ZVector<float>& B );

void laplacian(double* image,double* lap, int nx, int ny);
void horn_derives(double* image1, double* image2, double* dx, double* dy, double* dt, int nx,int ny);

template <class T>
void Threshold(ZMatrix<T>& Data, float percent)
{
	if(percent < 0 || percent > 1) percent = 1;
	int size = Data.NRows() * Data.NCols();
	ZVector<T> vec(size);
	typename ZMatrix<T>::uiterator p_d(Data);
	typename ZVector<T>::iterator p_v = vec.begin();
	for(; p_d.FMore(); p_v++, p_d++) *p_v = Abs<T>(*p_d);
	sort(vec.begin(), vec.end());
	T thres = vec[(size-1) * (1 - percent)];
	for(p_d.Reset(); p_d.FMore(); p_d++) if(Abs(*p_d)<thres) *p_d = 0;
}

#endif // __MATHCOMMON_H_
