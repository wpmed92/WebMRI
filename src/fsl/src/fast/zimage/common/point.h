/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares/defines the ZPoint class.

#ifndef __ZPOINT_H__
#define __ZPOINT_H__

#include <iostream> // C++ I/O
#include <cmath>

#include "mydefine.h"

template <class T>
class ZPoint
{
public:
	T	x,y,z;

	ZPoint(T xx = T(0), T yy = T(0), T zz = T(0)) : x(xx) , y(yy), z(zz) {}

	template <class T1>
	ZPoint(const ZPoint<T1>& pt) : x(T(pt.x)), y(T(pt.y)), z(T(pt.z)) {}
	
	template <class T1>
	ZPoint&		operator = (const ZPoint<T1> & pt) { x = T(pt.x); y = T(pt.y); z = T(pt.z); return *this; }

	T& operator[](int d) { return *(&x+d); }
	const T& operator[](int d) const { return *(&x+d); }

	T* begin() { return &x;} 
	const T* begin() const { return &x;}
	
	operator T* () { return &x; }
	operator const T* () const { return &x; }

	void		set(T xx, T yy, T zz=0) { x=xx; y=yy; z=zz; }
	void		clear() { x = y = z = 0; }
	double		Angle() const { return atan2(y, x); }
	double		AngleXY() const { return atan2(y, x); }
	double		AngleYZ() const { return atan2(z, y); }
	double		AngleXZ() const {return atan2(z, x); }
	
	ZPoint&		normalize() { *this /= len(*this); return *this; }
	ZPoint&		operator += (const ZPoint & v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
	ZPoint&		operator -= (const ZPoint & v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
	template <class T2>
	ZPoint&		operator *= (const ZPoint<T2> & v) { x*=v.x; y*=v.y; z*=v.z; return *this; }
	template <class T2>
	ZPoint&		operator /= (const ZPoint<T2> & v) { x/=v.x; y/=v.y; z/=v.z; return *this; }
	ZPoint&		operator += (T f) { x+=f; y+=f; z+=f; return *this; }
	ZPoint&		operator -= (T f) { x-=f; y-=f; z-=f; return *this; }
	ZPoint&		operator *= (float f) { x=T(x*f); y=T(y*f); z=T(z*f); return *this; }
	ZPoint&		operator /= (float f) { x=T(x/f); y=T(y/f); z=T(z/f); return *this; }
};	//class define over

template <class T> inline float sqrlen(const ZPoint<T>& v) { return dot(v, v); }
template <class T> inline float len(const ZPoint<T>& v) { return sqrt(sqrlen(v)); }
template <class T> inline ZPoint<T> norm(const ZPoint<T>& v) { return v / len(v); }
template <class T> inline T sum(const ZPoint<T>& v) { return v.x+v.y+v.z; }
template <class T> inline T Min(const ZPoint<T>& v) { return (v.x<v.y) ? Min(v.x, v.z) : Min(v.y, v.z); }
template <class T> inline T Max(const ZPoint<T>& v) { return (v.x>v.y) ? Max(v.x, v.z) : Max(v.y, v.z); }
template<class T> inline float mean(const ZPoint<T>& v) { return float(v.x+v.y+v.z)/3; }
template <class T> inline ZPoint<T> Abs(const ZPoint<T>& v) { return ZPoint<T>(Abs(v.x), Abs(v.y), Abs(v.z)); }
template <class T> inline float sqrdis(const ZPoint<T>& l, const ZPoint<T>& r) 
{ return (l.x-r.x)*(l.x-r.x) + (l.y-r.y)*(l.y-r.y) + (l.z-r.z)*(l.z-r.z); } 
template <class T> inline float dis(const ZPoint<T>& l, const ZPoint<T>& r) { return sqrt(sqrdis(l, r)); }
template <class T> inline float dot(const ZPoint<T>& l, const ZPoint<T>& r) { return l.x*r.x+l.y*r.y+l.z*r.z; }

template <class T> inline ZPoint<T> operator + (const ZPoint<T>& v) { return v; }
template <class T> inline ZPoint<T> operator - (const ZPoint<T>& v) { return ZPoint<T> (-v.x, -v.y, -v.z); }
template <class T> 
inline ZPoint<T> operator + (const ZPoint<T>& l, const ZPoint<T>& r) { ZPoint<T> v(l);  return v+=r; }
template <class T>
inline ZPoint<T> operator - (const ZPoint<T>& l, const ZPoint<T>& r) { ZPoint<T> v(l);  return v-=r; }
template <class T1, class T2>
inline ZPoint<T1> operator * (const ZPoint<T1>& l, const ZPoint<T2>& r) { ZPoint<T1> v(l);  return v*=r; }
template <class T1, class T2>
inline ZPoint<T1> operator / (const ZPoint<T1>& l, const ZPoint<T2>& r) { ZPoint<T1> v(l);  return v/=r; }
template <class T>
inline ZPoint<T> operator + (const ZPoint<T>& l, T f) { ZPoint<T> v(l); return v+=f; }
template <class T>
inline ZPoint<T> operator - (const ZPoint<T>& l, T f) { ZPoint<T> v(l); return v-=f; }
template <class T>
inline ZPoint<T> operator * (const ZPoint<T>& l, float f) { ZPoint<T> v(l); return v*=f; }
template <class T>
inline ZPoint<T> operator / (const ZPoint<T>& l, float f) { ZPoint<T> v(l); return v/=f; }
template <class T> inline bool	operator == (const ZPoint<T>& l, const ZPoint<T>& r)
{ return (l.x==r.x) && (l.y==r.y) && (l.z==r.z); }
template <class T> inline bool operator != (const ZPoint<T>& l, const ZPoint<T>& r) {return !(l==r); }

template <class T> inline ZPoint<double> acos(const ZPoint<T>& v) { return ZPoint<double>(acos(v.x), acos(v.y), acos(v.z)); }
template <class T> inline ZPoint<double> asin(const ZPoint<T>& v) { return ZPoint<double>(asin(v.x), asin(v.y), asin(v.z)); }
template <class T> inline ZPoint<double> atan(const ZPoint<T>& v) { return ZPoint<double>(atan(v.x), atan(v.y), atan(v.z)); }
template <class T> inline ZPoint<double> cos(const ZPoint<T>& v) { return ZPoint<double>(cos(v.x), cos(v.y), cos(v.z)); }
template <class T> inline ZPoint<double> sin(const ZPoint<T>& v) { return ZPoint<double>(sin(v.x), sin(v.y), sin(v.z)); }
template <class T> inline ZPoint<double> cosh(const ZPoint<T>& v) { return ZPoint<double>(cosh(v.x), cosh(v.y), cosh(v.z)); }
template <class T> inline ZPoint<double> sinh(const ZPoint<T>& v) { return ZPoint<double>(sinh(v.x), sinh(v.y), sinh(v.z)); }
template <class T> inline ZPoint<double> tan(const ZPoint<T>& v) { return ZPoint<double>(tan(v.x), tan(v.y), tan(v.z)); }
template <class T> inline ZPoint<double> tanh(const ZPoint<T>& v) { return ZPoint<double>(tanh(v.x), tanh(v.y), tanh(v.z)); }
template <class T> inline ZPoint<double> exp(const ZPoint<T>& v) { return ZPoint<double>(exp(v.x), exp(v.y), exp(v.z)); }
template <class T> inline ZPoint<double> log(const ZPoint<T>& v) { return ZPoint<double>(log(v.x), log(v.y), log(v.z)); }
template <class T> inline ZPoint<double> log10(const ZPoint<T>& v) { return ZPoint<double>(log10(v.x), log10(v.y), log10(v.z)); }
template <class T> inline ZPoint<double> sqrt(const ZPoint<T>& v) { return ZPoint<double>(sqrt(v.x), sqrt(v.y), sqrt(v.z)); }
template <class T> inline ZPoint<double> pow(const ZPoint<T>& l, const ZPoint<T>& r) 
{ return ZPoint<double>(pow(l.x, r.x), pow(l.y, r.y), pow(l.z, r.z)); }
template <class T> inline ZPoint<double> pow(const ZPoint<T>& v, const T& f) 
{ return ZPoint<double>(pow(v.x, f), pow(v.y, f), pow(v.z, f)); }
template <class T> inline ZPoint<double> pow(const T& f, const ZPoint<T>& v) 
{ return ZPoint<double>(pow(f, v.x), pow(f, v.y), pow(f, v.z)); }

template <class T>
inline std::ostream& operator << (std::ostream& s, const ZPoint<T> & p)
{
	s << p.x << ' ' << p.y << ' ' << p.z;;
	return s ;
}

typedef ZPoint<float> ZPoint3D;

#endif
