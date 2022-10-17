/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares 2 x 2 Matrix class for faster calculations than the
// general ZMatrix class.

#ifndef __MATRIX2D_H_
#define __MATRIX2D_H_

#include <algorithm>
#include "vectorD.h"

#undef _VALOP
#define _VALOP(OP) return ZMatrix2D(OP(v.d11), OP(v.d12), OP(v.d21), OP(v.d22)) 

class ZMatrix2D
{
public:
	float				d11, d12;
	float				d21, d22;

	typedef float* iterator;
	typedef const float* const_iterator;

	iterator begin() {return &d11; }
	const_iterator begin() const { return ((const_iterator)(&d11)); }
	iterator end() {return (&d22)+1; }
	const_iterator end() const { return (const_iterator)((&d22)+1); }

	ZMatrix2D () : d11(0), d12(0), d21(0), d22(0) {}
	ZMatrix2D (float d1, float d2, float d3, float d4) : d11(d1), d12(d2), d21(d3), d22(d4) {}
	ZMatrix2D (float *dat) { *this = *(ZMatrix2D*)dat; }
	ZMatrix2D (float dat[][2]) { *this = *(ZMatrix2D*)dat; }

	ZMatrix2D& operator = (float v) { d11 = d12 = v; d21 = d22 = v; return *this; }

	float&	operator () (int r, int c) { return *(&d11 + 2*r + c); }
	const float& operator() (int r, int c) const { return *(&d11 + 2*r + c); }

	ZMatrix2D&	normalize(int dir=2);	//column
	ZMatrix2D& MakeDiag(float v) { d11=d22=v; return *this; }

	ZMatrix2D& operator += (const ZMatrix2D& m) { d11 += m.d11; d12 += m.d12; d21 += m.d21; d22 += m.d22; return *this; }
	ZMatrix2D& operator -= (const ZMatrix2D& m) { d11 -= m.d11; d12 -= m.d12; d21 -= m.d21; d22 -= m.d22; return *this; }
	ZMatrix2D& operator *= (const ZMatrix2D& m) { d11 *= m.d11; d12 *= m.d12; d21 *= m.d21; d22 *= m.d22; return *this; }
	ZMatrix2D& operator /= (const ZMatrix2D& m) { d11 /= m.d11; d12 /= m.d12; d21 /= m.d21; d22 /= m.d22; return *this; }

	ZMatrix2D& operator += (float f) { d11 += f; d12 += f; d21 += f; d22 += f; return *this; }
	ZMatrix2D& operator -= (float f) { d11 -= f; d12 -= f; d21 -= f; d22 -= f; return *this; }
	ZMatrix2D& operator *= (float f) { d11 *= f; d12 *= f; d21 *= f; d22 *= f; return *this; }
	ZMatrix2D& operator /= (float f) { d11 /= f; d12 /= f; d21 /= f; d22 /= f; return *this; }
};

inline ZMatrix2D operator + (const ZMatrix2D& v) { return v; }
inline ZMatrix2D operator - (const ZMatrix2D& v) { _VALOP(-); }
inline ZMatrix2D operator + (const ZMatrix2D& l, const ZMatrix2D& r) { ZMatrix2D v(l); return v+=r; }
inline ZMatrix2D operator - (const ZMatrix2D& l, const ZMatrix2D& r) { ZMatrix2D v(l); return v-=r; }
inline ZMatrix2D operator / (const ZMatrix2D& l, const ZMatrix2D& r) { ZMatrix2D v(l); return v/=r; }
	
inline ZMatrix2D operator + (const ZMatrix2D& l, float r) { ZMatrix2D v(l); return v+=r; }
inline ZMatrix2D operator - (const ZMatrix2D& l, float r) { ZMatrix2D v(l); return v-=r; }
inline ZMatrix2D operator * (const ZMatrix2D& l, float r) { ZMatrix2D v(l); return v*=r; }
inline ZMatrix2D operator / (const ZMatrix2D& l, float r) { ZMatrix2D v(l); return v/=r; }

inline ZMatrix2D operator* (const ZMatrix2D& l, const ZMatrix2D& r)
{		
	ZMatrix2D	t;
	t.d11 = l.d11 * r.d11 + l.d12 * r.d21; 	t.d12 = l.d11 * r.d12 + l.d12 * r.d22;
	t.d21 = l.d21 * r.d11 + l.d22 * r.d21; 	t.d22 = l.d21 * r.d12 + l.d22 * r.d22;
	return t;
}
inline ZVector2D operator* (const ZMatrix2D& l, const ZVector2D& r)
{
	return ZVector2D(l.d11*r.x + l.d12*r.y, l.d21*r.x + l.d22*r.y);
}
inline ZVector2D operator* (const ZVector2D& l, const ZMatrix2D& r)
{
	return ZVector2D(l.x*r.d11 + l.y*r.d21, l.x*r.d12 + l.y*r.d22);
}

inline float det(const ZMatrix2D& v) { return v.d11*v.d22 - v.d12*v.d21; }
inline float sum(const ZMatrix2D& v) { return v.d11+v.d12+v.d21+v.d22; }
inline float Min(const ZMatrix2D& v) { return *(std::min_element(v.begin(), v.end())); }
inline float Max(const ZMatrix2D& v) { return *(std::max_element(v.begin(), v.end())); }
inline float mean(const ZMatrix2D& v) { return float(sum(v))/4; }
inline ZVector2D mean_row(const ZMatrix2D& v) { return ZVector2D((v.d11+v.d21)/2, (v.d12+v.d22)/2); }
inline ZVector2D mean_col(const ZMatrix2D& v) { return ZVector2D((v.d11+v.d12)/2, (v.d21+v.d22)/2); }
inline ZMatrix2D Abs(const ZMatrix2D& v) { _VALOP(Abs); }
inline ZMatrix2D trans(const ZMatrix2D& v) { return ZMatrix2D(v.d11, v.d21, v.d12, v.d22); }
inline float tracs(const ZMatrix2D& v) { return v.d11 + v.d22; }
inline ZMatrix2D cov(const ZMatrix2D& v) 
{
	float m1 = (v.d11+v.d21)/2, m2 = (v.d12+v.d22)/2;
	ZMatrix2D mm(v.d11-m1, v.d12-m2, v.d21-m2, v.d22-m2);
	ZMatrix2D cov = trans(mm) * mm;
	return cov /= 2;
}
inline ZMatrix2D inv(const ZMatrix2D& m)
{
	float v = det(m); if(v < TINY) v = 0;
	if (!v) { ZWarning("Matrix2D inv", "Not invertable!"); 	return ZMatrix2D(); }
	return ZMatrix2D(m.d22/v, -m.d12/v, -m.d21/v, m.d11/v); 
}

inline ZMatrix2D cos(const ZMatrix2D& v) { _VALOP(cos); }
inline ZMatrix2D sin(const ZMatrix2D& v) { _VALOP(sin); }
inline ZMatrix2D tan(const ZMatrix2D& v) { _VALOP(tan); }
inline ZMatrix2D cosh(const ZMatrix2D& v) { _VALOP(cosh); }
inline ZMatrix2D sinh(const ZMatrix2D& v) { _VALOP(sinh); }
inline ZMatrix2D tanh(const ZMatrix2D& v) { _VALOP(tanh); }
inline ZMatrix2D acos(const ZMatrix2D& v) { _VALOP(acos); }
inline ZMatrix2D asin(const ZMatrix2D& v) { _VALOP(asin); }
inline ZMatrix2D atan(const ZMatrix2D& v) { _VALOP(atan); }
inline ZMatrix2D sqrt(const ZMatrix2D& v) { _VALOP(sqrt); }
inline ZMatrix2D exp(const ZMatrix2D& v) { _VALOP(exp); }
inline ZMatrix2D log(const ZMatrix2D& v) { _VALOP(log); }
inline ZMatrix2D log10(const ZMatrix2D& v) { _VALOP(log10); }

#endif
