/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares 4 x 4 Matrix class for faster calculations than the
// general ZMatrix class.

#ifndef __MATRIX4D_H_
#define __MATRIX4D_H_

#include "vectorD.h"

class ZMatrix4D
{
public:
	float				d11, d12, d13, d14;
	float				d21, d22, d23, d24;
	float				d31, d32, d33, d34;
	float				d41, d42, d43, d44;

	typedef float* iterator;
	typedef const float* const_iterator;

	iterator begin() {return &d11; }
	const_iterator begin() const { return ((const_iterator)(&d11)); }
	iterator end() {return (&d44)+1; }
	const_iterator end() const { return (const_iterator)((&d44)+1); }
	
	ZMatrix4D () {for(iterator p=begin(); p!=end(); p++) *p = 0; }
	ZMatrix4D (float d1, float d2, float d3, float d4, float d5, float d6, float d7, float d8, 
			   float d9, float d10, float d11, float d12, float d13, float d14, float d15, float d16) 
		: d11(d1), d12(d2), d13(d3), d14(d4), d21(d5), d22(d6), d23(d7), d24(d8),
		  d31(d9), d32(d10), d33(d11), d34(d12), d41(d13), d42(d14), d43(d15), d44(d16) {}
	ZMatrix4D (float *dat) { *this = *(ZMatrix4D*)dat; }
	ZMatrix4D (float dat[][4]) { *this = *(ZMatrix4D*)dat; }

	ZMatrix4D& operator = (float v) 
	{ 
		d11 = d12 = d13 = d14 = v;  d21 = d22 = d23 = d24 = v; 
		d31 = d32 = d33 = d34 = v;  d41 = d42 = d43 = d44 = v; 
		return *this;
	}

	ZVector4D&	operator[](int d) { return *((ZVector4D*)this+d); }
	const ZVector4D&	operator[](int d) const { return *((ZVector4D*)this+d); }
	
	float&	operator () (int r, int c) { return *(&d11 + 4*r + c); }
	const float& operator() (int r, int c) const { return *(&d11 + 4*r + c); }

	ZMatrix4D&	normalize(int dir=2);	//column
	ZMatrix4D& MakeDiag(float v) { d11=d22=d33=d44=v; return *this; }

	ZMatrix4D& operator += (const ZMatrix4D& m);
	ZMatrix4D& operator -= (const ZMatrix4D& m);
	ZMatrix4D& operator *= (const ZMatrix4D& m);
	ZMatrix4D& operator /= (const ZMatrix4D& m);

	ZMatrix4D& operator += (float f);
	ZMatrix4D& operator -= (float f);
	ZMatrix4D& operator *= (float f);
	ZMatrix4D& operator /= (float f);
};

ZMatrix4D operator + (const ZMatrix4D& m);
ZMatrix4D operator - (const ZMatrix4D& m);
ZMatrix4D operator + (const ZMatrix4D& l, const ZMatrix4D& r);
ZMatrix4D operator - (const ZMatrix4D& l, const ZMatrix4D& r);
ZMatrix4D operator / (const ZMatrix4D& l, const ZMatrix4D& r);
	
ZMatrix4D operator + (const ZMatrix4D& l, float r);
ZMatrix4D operator - (const ZMatrix4D& l, float r);
ZMatrix4D operator * (const ZMatrix4D& l, float r);
ZMatrix4D operator / (const ZMatrix4D& l, float r);

ZMatrix4D operator* (const ZMatrix4D& l, const ZMatrix4D& r);
ZVector4D operator* (const ZMatrix4D& l, const ZVector4D& r);
ZVector4D operator* (const ZVector4D& l, const ZMatrix4D& r);

float det(const ZMatrix4D& v); 
float sum(const ZMatrix4D& v);
float Min(const ZMatrix4D& v);
float Max(const ZMatrix4D& v);
float mean(const ZMatrix4D& v);
ZVector4D mean_row(const ZMatrix4D& v); 
ZVector4D mean_col(const ZMatrix4D& v); 
ZMatrix4D Abs(const ZMatrix4D& v);
ZMatrix4D trans(const ZMatrix4D& v); 
inline float tracs(const ZMatrix4D& v) { return v.d11 + v.d22 + v.d33 + v.d44; }
ZMatrix4D cov(const ZMatrix4D& v);
ZMatrix4D inv(const ZMatrix4D& m);

ZMatrix4D cos(const ZMatrix4D& v);
ZMatrix4D sin(const ZMatrix4D& v);
ZMatrix4D tan(const ZMatrix4D& v);
ZMatrix4D cosh(const ZMatrix4D& v);
ZMatrix4D sinh(const ZMatrix4D& v);
ZMatrix4D tanh(const ZMatrix4D& v);
ZMatrix4D acos(const ZMatrix4D& v);
ZMatrix4D asin(const ZMatrix4D& v);
ZMatrix4D atan(const ZMatrix4D& v);
ZMatrix4D sqrt(const ZMatrix4D& v);
ZMatrix4D exp(const ZMatrix4D& v);
ZMatrix4D log(const ZMatrix4D& v);
ZMatrix4D log10(const ZMatrix4D& v);

#endif
