/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares 3 x 3 Matrix class for faster calculations than the
// general ZMatrix class.

#ifndef __MATRIX3D_H_
#define __MATRIX3D_H_

#include <algorithm>
#include "vectorD.h"

#undef _VALOP
#define _VALOP(OP) \
	return ZMatrix3D(OP(v.d11), OP(v.d12), OP(v.d13), \
					 OP(v.d21), OP(v.d22), OP(v.d23), \
					 OP(v.d31), OP(v.d32), OP(v.d33))

#undef _VALOPF
#define _VALOPF(RHS) \
	d11 RHS; d12 RHS; d13 RHS; \
	d21 RHS; d22 RHS; d23 RHS; \
	d31 RHS; d32 RHS; d33 RHS; \
	return *this
#undef _VALOPM
#define _VALOPM(OP, m) \
	d11 OP m.d11; d12 OP m.d12; d13 OP m.d13; \
	d21 OP m.d21; d22 OP m.d22; d23 OP m.d23; \
	d31 OP m.d31; d32 OP m.d32; d33 OP m.d33; \
	return *this

class ZMatrix3D
{
public:
	float				d11, d12, d13;
	float				d21, d22, d23;
	float				d31, d32, d33;

	typedef float* iterator;
	typedef const float* const_iterator;

	iterator begin() {return &d11; }
	const_iterator begin() const { return ((const_iterator)(&d11)); }
	iterator end() {return (&d33)+1; }
	const_iterator end() const { return (const_iterator)((&d33)+1); }
	
	ZMatrix3D () {for(iterator p=begin(); p!=end(); p++) *p = 0; }
	ZMatrix3D (float d1, float d2, float d3, float d4, float d5, float d6, float d7, float d8, float d9) 
		: d11(d1), d12(d2), d13(d3), d21(d4), d22(d5), d23(d6), d31(d7), d32(d8), d33(d9) {}
	ZMatrix3D (float *dat) { *this = *(ZMatrix3D*)dat; }
	ZMatrix3D (float dat[][3]) { *this = *(ZMatrix3D*)dat; }

	ZMatrix3D& operator = (float v) { d11 = d12 = d13 = v; 	d21 = d22 = d23 = v;  d31 = d32 = d33 = v; 	return *this; }

	ZVector3D&	operator[](int d) { return *((ZVector3D*)this+d); }
	const ZVector3D&	operator[](int d) const { return *((ZVector3D*)this+d); }
	
	float&	operator () (int r, int c) { return *(&d11 + 3*r + c); }
	const float& operator() (int r, int c) const { return *(&d11 + 3*r + c); }

	ZMatrix3D&	normalize(int dir=2);	//column
	ZMatrix3D& MakeDiag(float v) { d11=d22=d33=v; return *this; }

	ZMatrix3D& operator += (const ZMatrix3D& m) { _VALOPM(+=, m); }
	ZMatrix3D& operator -= (const ZMatrix3D& m) { _VALOPM(-=, m); }
	ZMatrix3D& operator *= (const ZMatrix3D& m) { _VALOPM(*=, m); }
	ZMatrix3D& operator /= (const ZMatrix3D& m) { _VALOPM(/=, m); }

	ZMatrix3D& operator += (float f) { _VALOPF(+= f); }
	ZMatrix3D& operator -= (float f) { _VALOPF(-= f); }
	ZMatrix3D& operator *= (float f) { _VALOPF(*= f); }
	ZMatrix3D& operator /= (float f) { _VALOPF(/= f); }
};

inline ZMatrix3D operator + (const ZMatrix3D& v) { return v; }
inline ZMatrix3D operator - (const ZMatrix3D& v) { _VALOP(-); }
inline ZMatrix3D operator + (const ZMatrix3D& l, const ZMatrix3D& r) { ZMatrix3D v(l); return v+=r; }
inline ZMatrix3D operator - (const ZMatrix3D& l, const ZMatrix3D& r) { ZMatrix3D v(l); return v-=r; }
inline ZMatrix3D operator / (const ZMatrix3D& l, const ZMatrix3D& r) { ZMatrix3D v(l); return v/=r; }
	
inline ZMatrix3D operator + (const ZMatrix3D& l, float r) { ZMatrix3D v(l); return v+=r; }
inline ZMatrix3D operator - (const ZMatrix3D& l, float r) { ZMatrix3D v(l); return v-=r; }
inline ZMatrix3D operator * (const ZMatrix3D& l, float r) { ZMatrix3D v(l); return v*=r; }
inline ZMatrix3D operator / (const ZMatrix3D& l, float r) { ZMatrix3D v(l); return v/=r; }

ZMatrix3D operator* (const ZMatrix3D& l, const ZMatrix3D& r);

inline ZVector3D operator* (const ZMatrix3D& l, const ZVector3D& r)
{
	return ZVector3D(l.d11*r.x+l.d12*r.y+l.d13*r.z, 
					 l.d21*r.x+l.d22*r.y+l.d23*r.z,
					 l.d31*r.x+l.d32*r.y+l.d33*r.z);
}
inline ZVector3D operator* (const ZVector3D& l, const ZMatrix3D& r)
{
	return ZVector3D(l.x*r.d11+l.y*r.d21+l.z*r.d31, 
					 l.x*r.d12+l.y*r.d22+l.z*r.d32,
					 l.x*r.d13+l.y*r.d23+l.z*r.d33);
}

inline float det(const ZMatrix3D& v) 
{  return v.d11*v.d22*v.d33+v.d12*v.d23*v.d31+v.d13*v.d21*v.d32-v.d13*v.d22*v.d31-v.d11*v.d23*v.d32-v.d12*v.d21*v.d33; }
inline float sum(const ZMatrix3D& v) { return v.d11+v.d12+v.d13+v.d21+v.d22+v.d23+v.d31+v.d32+v.d33; }
inline float Min(const ZMatrix3D& v) { return *(std::min_element(v.begin(), v.end())); }
inline float Max(const ZMatrix3D& v) { return *(std::max_element(v.begin(), v.end())); }
inline float mean(const ZMatrix3D& v) { return float(sum(v))/9; }
inline ZVector3D mean_row(const ZMatrix3D& v) 
{ return ZVector3D((v.d11+v.d21+v.d31)/3, (v.d12+v.d22+v.d32)/3, (v.d13+v.d23+v.d33)/3); }
inline ZVector3D mean_col(const ZMatrix3D& v) 
{ return ZVector3D((v.d11+v.d12+v.d13)/3, (v.d21+v.d22+v.d23)/3, (v.d31+v.d32+v.d33)/3); }
inline ZMatrix3D Abs(const ZMatrix3D& v) { _VALOP(Abs); }
inline ZMatrix3D trans(const ZMatrix3D& v) 
{ return ZMatrix3D(v.d11, v.d21, v.d31, v.d12, v.d22, v.d32, v.d13, v.d23, v.d33); }
inline float tracs(const ZMatrix3D& v) { return v.d11 + v.d22 + v.d33; }
ZMatrix3D cov(const ZMatrix3D& v);
ZMatrix3D inv(const ZMatrix3D& m);

inline ZMatrix3D cos(const ZMatrix3D& v) { _VALOP(cos); }
inline ZMatrix3D sin(const ZMatrix3D& v) { _VALOP(sin); }
inline ZMatrix3D tan(const ZMatrix3D& v) { _VALOP(tan); }
inline ZMatrix3D cosh(const ZMatrix3D& v) { _VALOP(cosh); }
inline ZMatrix3D sinh(const ZMatrix3D& v) { _VALOP(sinh); }
inline ZMatrix3D tanh(const ZMatrix3D& v) { _VALOP(tanh); }
inline ZMatrix3D acos(const ZMatrix3D& v) { _VALOP(acos); }
inline ZMatrix3D asin(const ZMatrix3D& v) { _VALOP(asin); }
inline ZMatrix3D atan(const ZMatrix3D& v) { _VALOP(atan); }
inline ZMatrix3D sqrt(const ZMatrix3D& v) { _VALOP(sqrt); }
inline ZMatrix3D exp(const ZMatrix3D& v) { _VALOP(exp); }
inline ZMatrix3D log(const ZMatrix3D& v) { _VALOP(log); }
inline ZMatrix3D log10(const ZMatrix3D& v) { _VALOP(log10); }

#endif
