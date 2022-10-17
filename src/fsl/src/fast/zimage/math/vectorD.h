/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares 2D, 3D, and 4D Vector class for faster calculations 
// than the general ZVector class.

#ifndef __VECTORD_H_
#define __VECTORD_H_

#include "common.h"

class ZVector2D : public ZPoint<float>
{
public:
	ZVector2D()	: ZPoint<float>() {}
	ZVector2D(float xx, float yy) : ZPoint<float>(xx, yy) {}
	ZVector2D(const ZPoint<float>& pt) : ZPoint<float>(pt) {}

	ZVector2D&	normalize() { float l=len(*this); if(l>TINY) *this /= len(*this); return *this; }
	ZVector2D&	operator += (const ZVector2D & v) { x+=v.x; y+=v.y; return *this; }
	ZVector2D&	operator -= (const ZVector2D & v) { x-=v.x; y-=v.y; return *this; }
	ZVector2D&	operator *= (const ZVector2D & v) { x*=v.x; y*=v.y; return *this; }
	ZVector2D&	operator /= (const ZVector2D & v) { x/=v.x; y/=v.y; return *this; }
	ZVector2D&	operator += (float f) { x+=f; y+=f; return *this; }
	ZVector2D&	operator -= (float f) { x-=f; y-=f; return *this; }
	ZVector2D&	operator *= (float f) { x*=f; y*=f; return *this; }
	ZVector2D&	operator /= (float f) { x/=f; y/=f; return *this; }

	ZVector2D&	operator = (const ZPoint<float>& pt) { x=pt.x; y=pt.y; return *this; }
};

inline ZVector2D operator + (const ZVector2D& v) { return v; }
inline ZVector2D operator - (const ZVector2D& v) { return ZVector2D (-v.x, -v.y); }
inline ZVector2D operator + (const ZVector2D& l, const ZVector2D& r) { ZVector2D v(l); return v += r; }
inline ZVector2D operator - (const ZVector2D& l, const ZVector2D& r) { ZVector2D v(l); return v -= r; }
inline ZVector2D operator * (const ZVector2D& l, const ZVector2D& r) { ZVector2D v(l); return v *= r; }
inline ZVector2D operator / (const ZVector2D& l, const ZVector2D& r) { ZVector2D v(l); return v /= r; }
inline ZVector2D operator + (const ZVector2D& l, float f) { ZVector2D v(l); return v += f; }
inline ZVector2D operator - (const ZVector2D& l, float f) { ZVector2D v(l); return v -= f; }
inline ZVector2D operator * (const ZVector2D& l, float f) { ZVector2D v(l); return v *= f; }
inline ZVector2D operator / (const ZVector2D& l, float f) { ZVector2D v(l); return v /= f; }
inline bool	operator == (const ZVector2D& l, const ZVector2D& r) { return (l.x==r.x) && (l.y==r.y); }
inline bool operator != (const ZVector2D& l, const ZVector2D& r) {return !(l==r); }

inline ZVector2D norm(const ZVector2D& v) { return v / len(v); }
inline float sum(const ZVector2D& v) { return v.x+v.y; }
inline float Min(const ZVector2D& v) { return Min(v.x, v.y); }
inline float Max(const ZVector2D& v) { return Max(v.x, v.y); }
inline float mean(const ZVector2D& v) { return (v.x+v.y)/2; }
inline ZVector2D Abs(const ZVector2D& v) { return ZVector2D(Abs(v.x), Abs(v.y)); }
inline float sqrdis(const ZVector2D& l, const ZVector2D& r) { return (l.x-r.x)*(l.x-r.x) + (l.y-r.y)*(l.y-r.y); }
inline float dis(const ZVector2D& l, const ZVector2D& r) { return sqrt(sqrdis(l, r)); }
inline float dot(const ZVector2D& l, const ZVector2D& r) { return l.x*r.x+l.y*r.y; }
inline float sqrlen(const ZVector2D& v) { return dot(v, v); }
inline float len(const ZVector2D& v) { return sqrt(sqrlen(v)); }
inline ZVector2D rotate(const ZVector2D& v, float angle) 
{ double c=cos(angle), s=sin(angle); return ZVector2D(v.x*c-v.y*s, v.x*s+v.y*c); }
inline ZVector2D mirror(const ZVector2D& v, const ZVector2D& a) 
{ 
	float l=len(a); if(l<TINY) return v; ZVector2D t=a/l; l = t.x*t.x-t.y*t.y;
	return ZVector2D(2*t.x*t.y*v.y+l*v.x, 2*t.x*t.y*v.x-l*v.y); 
}

inline ZVector2D acos(const ZVector2D& v) { return ZVector2D(acos(v.x), acos(v.y)); }
inline ZVector2D asin(const ZVector2D& v) { return ZVector2D(asin(v.x), asin(v.y)); }
inline ZVector2D atan(const ZVector2D& v) { return ZVector2D(atan(v.x), atan(v.y)); }
inline ZVector2D cos(const ZVector2D& v) { return ZVector2D(cos(v.x), cos(v.y)); }
inline ZVector2D sin(const ZVector2D& v) { return ZVector2D(sin(v.x), sin(v.y)); }
inline ZVector2D cosh(const ZVector2D& v) { return ZVector2D(cosh(v.x), cosh(v.y)); }
inline ZVector2D sinh(const ZVector2D& v) { return ZVector2D(sinh(v.x), sinh(v.y)); }
inline ZVector2D tan(const ZVector2D& v) { return ZVector2D(tan(v.x), tan(v.y)); }
inline ZVector2D tanh(const ZVector2D& v) { return ZVector2D(tanh(v.x), tanh(v.y)); }
inline ZVector2D exp(const ZVector2D& v) { return ZVector2D(exp(v.x), exp(v.y)); }
inline ZVector2D log(const ZVector2D& v) { return ZVector2D(log(v.x), log(v.y)); }
inline ZVector2D log10(const ZVector2D& v) { return ZVector2D(log10(v.x), log10(v.y)); }
inline ZVector2D sqrt(const ZVector2D& v) { return ZVector2D(sqrt(v.x), sqrt(v.y)); }
inline ZVector2D pow(const ZVector2D& l, const ZVector2D& r) { return ZVector2D(pow(l.x, r.x), pow(l.y, r.y)); }
inline ZVector2D pow(const ZVector2D& v, float f) { return ZVector2D(pow(v.x, f), pow(v.y, f)); }
inline ZVector2D pow(float f, const ZVector2D& v) { return ZVector2D(pow(f, v.x), pow(f, v.y)); }

class ZVector3D : public ZPoint<float>
{
public:
	ZVector3D()	: ZPoint<float>() {}
	ZVector3D(float xx, float yy, float zz) : ZPoint<float>(xx, yy, zz) {}
	ZVector3D(const ZPoint3D& pt) : ZPoint<float>(pt) {}

	ZVector3D&	normalize() { float l=len(*this); if(l>TINY) *this /= len(*this); return *this; }
	ZVector3D&	operator += (const ZVector3D & v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
	ZVector3D&	operator -= (const ZVector3D & v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
	ZVector3D&	operator *= (const ZVector3D & v) { x*=v.x; y*=v.y; z*=v.z; return *this; }
	ZVector3D&	operator /= (const ZVector3D & v) { x/=v.x; y/=v.y; z/=v.z; return *this; }
	ZVector3D&	operator += (float f) { x+=f; y+=f; z+=f; return *this; }
	ZVector3D&	operator -= (float f) { x-=f; y-=f; z-=f; return *this; }
	ZVector3D&	operator *= (float f) { x*=f; y*=f; z*=f; return *this; }
	ZVector3D&	operator /= (float f) { x/=f; y/=f; z/=f; return *this; }

	ZVector3D&	operator = (const ZPoint3D& pt) { x=pt.x; y=pt.y; z=pt.z; return *this; }
};

inline ZVector3D operator + (const ZVector3D& v) { return v; }
inline ZVector3D operator - (const ZVector3D& v) { return ZVector3D (-v.x, -v.y, -v.z); }
inline ZVector3D operator + (const ZVector3D& l, const ZVector3D& r) { ZVector3D v(l); return v += r; }
inline ZVector3D operator - (const ZVector3D& l, const ZVector3D& r) { ZVector3D v(l); return v -= r; }
inline ZVector3D operator * (const ZVector3D& l, const ZVector3D& r) { ZVector3D v(l); return v *= r; }
inline ZVector3D operator / (const ZVector3D& l, const ZVector3D& r) { ZVector3D v(l); return v /= r; }
inline ZVector3D operator + (const ZVector3D& l, float f) { ZVector3D v(l); return v += f; }
inline ZVector3D operator - (const ZVector3D& l, float f) { ZVector3D v(l); return v -= f; }
inline ZVector3D operator * (const ZVector3D& l, float f) { ZVector3D v(l); return v *= f; }
inline ZVector3D operator / (const ZVector3D& l, float f) { ZVector3D v(l); return v /= f; }
inline bool	operator == (const ZVector3D& l, const ZVector3D& r) { return (l.x==r.x) && (l.y==r.y) && (l.z==r.z); }
inline bool operator != (const ZVector3D& l, const ZVector3D& r) {return !(l==r); }

inline ZVector3D norm(const ZVector3D& v) { return v / len(v); }
inline float sum(const ZVector3D& v) { return v.x+v.y+v.z; }
inline float Min(const ZVector3D& v) { return (v.x<v.y) ? Min(v.x, v.z) : Min(v.y, v.z); }
inline float Max(const ZVector3D& v) { return (v.x>v.y) ? Max(v.x, v.z) : Max(v.y, v.z); }
inline float mean(const ZVector3D& v) { return (v.x+v.y+v.z)/3; }
inline ZVector3D Abs(const ZVector3D& v) { return ZVector3D(Abs(v.x), Abs(v.y), Abs(v.z)); }
inline float sqrdis(const ZVector3D& l, const ZVector3D& r) { return (l.x-r.x)*(l.x-r.x) + (l.y-r.y)*(l.y-r.y) + (l.z-r.z)*(l.z-r.z); }
inline float dis(const ZVector3D& l, const ZVector3D& r) { return sqrt(sqrdis(l, r)); }
inline float dot(const ZVector3D& l, const ZVector3D& r) { return l.x*r.x+l.y*r.y+l.z*r.z; }
inline ZVector3D cross(const ZVector3D& l, const ZVector3D& r) { return ZVector3D(l.y*r.z - l.z*r.y, l.z*r.x - l.x*r.z, l.x*r.y - l.y*r.x); }
inline float sqrlen(const ZVector3D& v) { return dot(v, v); }
inline float len(const ZVector3D& v) { return sqrt(sqrlen(v)); }

inline ZVector3D acos(const ZVector3D& v) { return ZVector3D(acos(v.x), acos(v.y), acos(v.z)); }
inline ZVector3D asin(const ZVector3D& v) { return ZVector3D(asin(v.x), asin(v.y), asin(v.z)); }
inline ZVector3D atan(const ZVector3D& v) { return ZVector3D(atan(v.x), atan(v.y), atan(v.z)); }
inline ZVector3D cos(const ZVector3D& v) { return ZVector3D(cos(v.x), cos(v.y), cos(v.z)); }
inline ZVector3D sin(const ZVector3D& v) { return ZVector3D(sin(v.x), sin(v.y), sin(v.z)); }
inline ZVector3D cosh(const ZVector3D& v) { return ZVector3D(cosh(v.x), cosh(v.y), cosh(v.z)); }
inline ZVector3D sinh(const ZVector3D& v) { return ZVector3D(sinh(v.x), sinh(v.y), sinh(v.z)); }
inline ZVector3D tan(const ZVector3D& v) { return ZVector3D(tan(v.x), tan(v.y), tan(v.z)); }
inline ZVector3D tanh(const ZVector3D& v) { return ZVector3D(tanh(v.x), tanh(v.y), tanh(v.z)); }
inline ZVector3D exp(const ZVector3D& v) { return ZVector3D(exp(v.x), exp(v.y), exp(v.z)); }
inline ZVector3D log(const ZVector3D& v) { return ZVector3D(log(v.x), log(v.y), log(v.z)); }
inline ZVector3D log10(const ZVector3D& v) { return ZVector3D(log10(v.x), log10(v.y), log10(v.z)); }
inline ZVector3D sqrt(const ZVector3D& v) { return ZVector3D(sqrt(v.x), sqrt(v.y), sqrt(v.z)); }
inline ZVector3D pow(const ZVector3D& l, const ZVector3D& r) { return ZVector3D(pow(l.x, r.x), pow(l.y, r.y), pow(l.z, r.z)); }
inline ZVector3D pow(const ZVector3D& v, float f) { return ZVector3D(pow(v.x, f), pow(v.y, f), pow(v.z, f)); }
inline ZVector3D pow(float f, const ZVector3D& v) { return ZVector3D(pow(f, v.x), pow(f, v.y), pow(f, v.z)); }

class ZVector4D : public ZPoint<float>
{
public:
	float w;

	ZVector4D()	: ZPoint<float>(), w(0) {}
	ZVector4D (float xx, float yy, float zz, float ww) : ZPoint<float>(xx, yy, zz), w(ww) {}
	ZVector4D (const ZPoint3D& pt) : ZPoint<float>(pt), w(1) {}

	void		set(float xx, float yy, float zz, float ww) { x=xx; y=yy; z=zz; w=ww; }
	ZVector4D&	normalize() { float l=len(*this); if(l>TINY) *this /= len(*this); return *this; }
	
	ZVector4D&	operator += (const ZVector4D & v) { x+=v.x; y+=v.y; z+=v.z; w+=v.w; return *this; }
	ZVector4D&	operator -= (const ZVector4D & v) { x-=v.x; y-=v.y; z-=v.z; w-=v.w; return *this; }
	ZVector4D&	operator *= (const ZVector4D & v) { x*=v.x; y*=v.y; z*=v.z; w*=v.w; return *this; }
	ZVector4D&	operator /= (const ZVector4D & v) { x/=v.x; y/=v.y; z/=v.z; w/=v.w; return *this; }
	ZVector4D&	operator += (float f) { x+=f; y+=f; z+=f; w+=f; return *this; }
	ZVector4D&	operator -= (float f) { x-=f; y-=f; z-=f; w-=f; return *this; }
	ZVector4D&	operator *= (float f) { x*=f; y*=f; z*=f; w*=f; return *this; }
	ZVector4D&	operator /= (float f) { x/=f; y/=f; z/=f; w/=f; return *this; }

	ZVector4D&	operator = (const ZPoint3D& pt) { x=pt.x; y=pt.y; z=pt.z; w=1; return *this; }
};

inline ZVector4D operator + (const ZVector4D& v) { return v; }
inline ZVector4D operator - (const ZVector4D& v) { return ZVector4D (-v.x, -v.y, -v.z, -v.w); }
inline ZVector4D operator + (const ZVector4D& l, const ZVector4D& r) { ZVector4D v(l); return v += r; }
inline ZVector4D operator - (const ZVector4D& l, const ZVector4D& r) { ZVector4D v(l); return v -= r; }
inline ZVector4D operator * (const ZVector4D& l, const ZVector4D& r) { ZVector4D v(l); return v *= r; }
inline ZVector4D operator / (const ZVector4D& l, const ZVector4D& r) { ZVector4D v(l); return v /= r; }
inline ZVector4D operator + (const ZVector4D& l, float f) { ZVector4D v(l); return v += f; }
inline ZVector4D operator - (const ZVector4D& l, float f) { ZVector4D v(l); return v -= f; }
inline ZVector4D operator * (const ZVector4D& l, float f) { ZVector4D v(l); return v *= f; }
inline ZVector4D operator / (const ZVector4D& l, float f) { ZVector4D v(l); return v /= f; }
inline bool	operator == (const ZVector4D& l, const ZVector4D& r) { return (l.x==r.x)&&(l.y==r.y)&&(l.z==r.z)&&(l.w==r.w); }
inline bool operator != (const ZVector4D& l, const ZVector4D& r) {return !(l==r); }

inline ZVector4D norm(const ZVector4D& v) { return v / len(v); }
inline float sum(const ZVector4D& v) { return v.x+v.y+v.z+v.w; }
inline float Min(const ZVector4D& v) { return Min(Min(v.x, v.y), Min(v.z, v.w)); }
inline float Max(const ZVector4D& v) { return Min(Max(v.x, v.y), Max(v.z, v.w)); }
inline float mean(const ZVector4D& v) { return (v.x+v.y+v.z+v.w)/4; }
inline ZVector4D Abs(const ZVector4D& v) { return ZVector4D(Abs(v.x), Abs(v.y), Abs(v.z), Abs(v.w)); }
inline float sqrdis(const ZVector4D& l, const ZVector4D& r) { return (l.x-r.x)*(l.x-r.x)+(l.y-r.y)*(l.y-r.y)+(l.z-r.z)*(l.z-r.z)+(l.w-r.w)*(l.w-r.w); }
inline float dis(const ZVector4D& l, const ZVector4D& r) { return sqrt(sqrdis(l, r)); }
inline float dot(const ZVector4D& l, const ZVector4D& r) { return l.x*r.x+l.y*r.y+l.z*r.z+l.w*r.w; }
inline float sqrlen(const ZVector4D& v) { return dot(v, v); }
inline float len(const ZVector4D& v) { return sqrt(sqrlen(v)); }

inline ZVector4D acos(const ZVector4D& v) { return ZVector4D(acos(v.x), acos(v.y), acos(v.z), acos(v.w)); }
inline ZVector4D asin(const ZVector4D& v) { return ZVector4D(asin(v.x), asin(v.y), asin(v.z), asin(v.w)); }
inline ZVector4D atan(const ZVector4D& v) { return ZVector4D(atan(v.x), atan(v.y), atan(v.z), atan(v.w)); }
inline ZVector4D cos(const ZVector4D& v) { return ZVector4D(cos(v.x), cos(v.y), cos(v.z), cos(v.w)); }
inline ZVector4D sin(const ZVector4D& v) { return ZVector4D(sin(v.x), sin(v.y), sin(v.z), sin(v.w)); }
inline ZVector4D cosh(const ZVector4D& v) { return ZVector4D(cosh(v.x), cosh(v.y), cosh(v.z), cosh(v.w)); }
inline ZVector4D sinh(const ZVector4D& v) { return ZVector4D(sinh(v.x), sinh(v.y), sinh(v.z), sinh(v.w)); }
inline ZVector4D tan(const ZVector4D& v) { return ZVector4D(tan(v.x), tan(v.y), tan(v.z), tan(v.w)); }
inline ZVector4D tanh(const ZVector4D& v) { return ZVector4D(tanh(v.x), tanh(v.y), tanh(v.z), tanh(v.w)); }
inline ZVector4D exp(const ZVector4D& v) { return ZVector4D(exp(v.x), exp(v.y), exp(v.z), exp(v.w)); }
inline ZVector4D log(const ZVector4D& v) { return ZVector4D(log(v.x), log(v.y), log(v.z), log(v.w)); }
inline ZVector4D log10(const ZVector4D& v) { return ZVector4D(log10(v.x), log10(v.y), log10(v.z), log10(v.w)); }
inline ZVector4D sqrt(const ZVector4D& v) { return ZVector4D(sqrt(v.x), sqrt(v.y), sqrt(v.z), sqrt(v.w)); }
inline ZVector4D pow(const ZVector4D& l, const ZVector4D& r) { return ZVector4D(pow(l.x, r.x), pow(l.y, r.y), pow(l.z, r.z), pow(l.w, r.w)); }
inline ZVector4D pow(const ZVector4D& v, float f) { return ZVector4D(pow(v.x, f), pow(v.y, f), pow(v.z, f), pow(v.w, f)); }
inline ZVector4D pow(float f, const ZVector4D& v) { return ZVector4D(pow(f, v.x), pow(f, v.y), pow(f, v.z), pow(f, v.w)); }

#endif
