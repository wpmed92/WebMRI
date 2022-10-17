/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file defines the class for RGB pixels. The class can be used as
// a special pixel type.

#ifndef ___RGBPIXELS_H__
#define ___RGBPIXELS_H__

#include <iostream>
#include <cmath>

#include "common.h"

// Forward declaration

template <class T> struct ZRGB;

// @type ZRGBCharPixel | <c ZRGB>\<c signed char>
typedef ZRGB<signed char>		ZRGBChar;

// @type ZRGBShortPixel | <c ZRGB>\<short>
typedef ZRGB<short>			ZRGBShort;

// @type ZRGBIntPixel | <c ZRGB>\<int>
typedef ZRGB<int>				ZRGBInt;

// @type ZRGBLongPixel | <c ZRGB>\<long>
typedef ZRGB<long>			ZRGBLong;


// @type ZRGBBytePixel | <c ZRGB>\<unsigned char>
typedef ZRGB<unsigned char>	ZRGBByte;
typedef ZRGB<unsigned char>	ZRGBUChar;

// @type ZRGBUShortPixel | <c ZRGB>\<unsigned short>
typedef ZRGB<unsigned short>	ZRGBUShort;

// @type ZRGBUIntPixel | <c ZRGB>\<unsigned int>
typedef ZRGB<unsigned int>	ZRGBUInt;

// @type ZRGBULongPixel | <c ZRGB>\<unsigned long>
typedef ZRGB<unsigned long>	ZRGBULong;


// @type ZRGBFloatPixel | <c ZRGB>\<float>
typedef ZRGB<float>			ZRGBFloat;

// @type ZRGBDoublePixel | <c ZRGB>\<float>
typedef ZRGB<float>			ZRGBDouble;

template <class T> 
struct ZRGB	
{
	T	r,g,b;

	ZRGB() : r(0), g(0), b(0) {}
	ZRGB(float n) : r(T(n)), g(T(n)), b(T(n)) {}
	ZRGB(float rr, float gg, float bb) : r(T(rr)), g(T(gg)), b(T(bb)) {}

	template <class T2>
	ZRGB(const ZRGB<T2>& n) : r(T(n.r)), g(T(n.g)), b(T(n.b)) {}
	
	template <class T2>
	ZRGB&	operator = (const ZRGB<T2>& n) { r = T(n.r), g = T(n.g), b = T(n.b); return *this; }
	ZRGB&	operator = (T n) { r = n, g = n, b = n; return *this; }

	void	set(T rr, T gg, T bb) { r=rr; g=gg; b=bb; }
	
	T		GrayValue() const { return T((r*77 + g*151 + b*28)/256); }
	operator float () const { return GrayValue(); }

	T&	operator () (UINT channel)
	{
		if(channel == 2) return g;
		if(channel == 3) return b;
		return r;
	}

	T   operator () (UINT channel=0) const
	{
		if(channel == 1) return r;
		if(channel == 2) return g;
		if(channel == 3) return b;
		return GrayValue();
	}

	ZRGB&	normalize() { *this /= len(*this); return *this; }

	template <class T2>
	ZRGB&	operator += (const ZRGB<T2>& v) { r+=v.r; g+=v.g; b+=v.b; return *this; }
	template <class T2>
	ZRGB&	operator -= (const ZRGB<T2>& v) { r-=v.r; g-=v.g; b-=v.b; return *this; }
	template <class T2>
	ZRGB&	operator *= (const ZRGB<T2>& v) { r*=v.r; g*=v.g; b*=v.b; return *this; }
	template <class T2>
	ZRGB&	operator /= (const ZRGB<T2>& v) { r/=v.r; g/=v.g; b/=v.b; return *this; }
	ZRGB&	operator += (float f) { r+=f; g+=f; b+=f; return *this; }
	ZRGB&	operator -= (float f) { r-=f; g-=f; b-=f; return *this; }
	ZRGB&	operator *= (float f) { r=T(r*f); g=T(g*f); b=T(b*f); return *this; }
	ZRGB&	operator /= (float f) { r=T(r/f); g=T(g/f); b=T(b/f); return *this; }

	ZRGB&	operator ++ (void) { r++; g++, b++; return *this; }
	ZRGB&	operator ++ (int) { r++; g++, b++; return *this; }
	ZRGB&	operator -- (void) { r--; g--, b--; return *this; }
	ZRGB&	operator -- (int) { r--; g--, b--; return *this; }
};

template <class T> inline float sqrlen(const ZRGB<T>& v) { return dot(v, v); }
template <class T> inline float len(const ZRGB<T>& v) { return sqrt(sqrlen(v)); }
template <class T> inline ZRGB<T> norm(const ZRGB<T>& v) { return v / len(v); }
template <class T> inline T sum(const ZRGB<T>& v) { return v.r+v.g+v.b; }
template <class T> inline T Min(const ZRGB<T>& v) { return (v.r<v.g) ? Min(v.r, v.b) : Min(v.g, v.b); }
template <class T> inline T Max(const ZRGB<T>& v) { return (v.r>v.g) ? Max(v.r, v.b) : Max(v.g, v.b); }
template<class T> inline float mean(const ZRGB<T>& v) { return float(v.r+v.g+v.b)/3; }
template <class T> inline ZRGB<T> Abs(const ZRGB<T>& v) { return ZRGB<T>(Abs(v.r), Abs(v.g), Abs(v.b)); }
template <class T> inline float sqrdis(const ZRGB<T>& x, const ZRGB<T>& y) 
{ return (x.r-y.r)*(x.r-y.r) + (x.g-y.g)*(x.g-y.g) + (x.b-y.b)*(x.b-y.b); } 
template <class T> inline float dis(const ZRGB<T>& x, const ZRGB<T>& y) { return sqrt(sqrdis(x, y)); }
template <class T> inline float dot(const ZRGB<T>& x, const ZRGB<T>& y) { return x.r*y.r+x.g*y.g+x.b*y.b; }

template <class T> inline ZRGB<T> operator + (const ZRGB<T>& x) { return x; }
template <class T> inline ZRGB<T> operator - (const ZRGB<T>& x) { return ZRGB<T>(-x.r, -x.g, -x.b); }
template <class T> inline ZRGB<T> operator + (const ZRGB<T>& x, float y) { return ZRGB<T> (x.r+y, x.g+y, x.b+y); }
template <class T> inline ZRGB<T> operator - (const ZRGB<T>& x, float y) { return ZRGB<T> (x.r-y, x.g-y, x.b-y); }
template <class T> inline ZRGB<T> operator * (const ZRGB<T>& x, float y) { return ZRGB<T> (x.r*y, x.g*y, x.b*y); }
template <class T> inline ZRGB<T> operator / (const ZRGB<T>& x, float y) { return ZRGB<T> (x.r/y, x.g/y, x.b/y); }
template <class T1, class T2> inline ZRGB<T1> operator + (const ZRGB<T1>& x, const ZRGB<T2>& y) { ZRGB<T1> v(x);  return v+=y; }
template <class T1, class T2> inline ZRGB<T1> operator - (const ZRGB<T1>& x, const ZRGB<T2>& y) { ZRGB<T1> v(x);  return v-=y; }
template <class T1, class T2> inline ZRGB<T1> operator * (const ZRGB<T1>& x, const ZRGB<T2>& y) { ZRGB<T1> v(x);  return v*=y; }
template <class T1, class T2> inline ZRGB<T1> operator / (const ZRGB<T1>& x, const ZRGB<T2>& y) { ZRGB<T1> v(x);  return v/=y; }
template <class T> inline bool operator == (const ZRGB<T>& x, const ZRGB<T>& y) { return (x.r==y.r) && (x.g==y.g) && (x.b==y.b); }
template <class T> inline bool operator != (const ZRGB<T>& x, const ZRGB<T>& y) {return !(x==y); }

template <class T> inline ZRGB<double> acos(const ZRGB<T>& v) { return ZRGB<double>(acos(v.r), acos(v.g), acos(v.b)); }
template <class T> inline ZRGB<double> asin(const ZRGB<T>& v) { return ZRGB<double>(asin(v.r), asin(v.g), asin(v.b)); }
template <class T> inline ZRGB<double> atan(const ZRGB<T>& v) { return ZRGB<double>(atan(v.r), atan(v.g), atan(v.b)); }
template <class T> inline ZRGB<double> cos(const ZRGB<T>& v) { return ZRGB<double>(cos(v.r), cos(v.g), cos(v.b)); }
template <class T> inline ZRGB<double> sin(const ZRGB<T>& v) { return ZRGB<double>(sin(v.r), sin(v.g), sin(v.b)); }
template <class T> inline ZRGB<double> cosh(const ZRGB<T>& v) { return ZRGB<double>(cosh(v.r), cosh(v.g), cosh(v.b)); }
template <class T> inline ZRGB<double> sinh(const ZRGB<T>& v) { return ZRGB<double>(sinh(v.r), sinh(v.g), sinh(v.b)); }
template <class T> inline ZRGB<double> tan(const ZRGB<T>& v) { return ZRGB<double>(tan(v.r), tan(v.g), tan(v.b)); }
template <class T> inline ZRGB<double> tanh(const ZRGB<T>& v) { return ZRGB<double>(tanh(v.r), tanh(v.g), tanh(v.b)); }
template <class T> inline ZRGB<double> exp(const ZRGB<T>& v) { return ZRGB<double>(exp(v.r), exp(v.g), exp(v.b)); }
template <class T> inline ZRGB<double> log(const ZRGB<T>& v) { return ZRGB<double>(log(v.r), log(v.g), log(v.b)); }
template <class T> inline ZRGB<double> log10(const ZRGB<T>& v) { return ZRGB<double>(log10(v.r), log10(v.g), log10(v.b)); }
template <class T> inline ZRGB<double> sqrt(const ZRGB<T>& v) { return ZRGB<double>(sqrt(v.r), sqrt(v.g), sqrt(v.b)); }
template <class T> inline ZRGB<double> pow(const ZRGB<T>& x, const ZRGB<T>& y) 
{ return ZRGB<double>(pow(x.r, y.r), pow(x.g, y.g), pow(x.b, y.b)); }
template <class T> inline ZRGB<double> pow(const ZRGB<T>& v, const T& f) 
{ return ZRGB<double>(pow(v.r, f), pow(v.g, f), pow(v.b, f)); }
template <class T> inline ZRGB<double> pow(const T& f, const ZRGB<T>& v) 
{ return ZRGB<double>(pow(f, v.r), pow(f, v.g), pow(f, v.b)); }

template <class T>
inline std::ostream& operator << (std::ostream& s, const ZRGB<T>& n)
{
	s << n.r << " " << n.g << " " << n.b;
	return s ;
}

template <class T>
inline std::istream& operator >> (std::istream& s, ZRGB<T>& num)
{
	int a, b, c;  s >> a >> b >> c;
	num.r = T(a); num.g = T(b); num.b = T(c);
	return s ;
}

/*
template <class T1, class T2>
inline bool operator < (const ZRGB<T1>& x, const ZRGB<T2>& y) 
{ return x.r < y.r && x.g < y.g && x.b < y.b ; }

template <class T1, class T2>
inline bool operator > (const ZRGB<T1>& x, const ZRGB<T2>& y) 
{ return x.r > y.r && x.g > y.g && x.b > y.b ; }

template <class T1, class T2>
inline bool operator <= (const ZRGB<T1>& x, const ZRGB<T2>& y) 
{ return x.r <= y.r && x.g <= y.g && x.b <= y.b ; }

template <class T1, class T2>
inline bool operator >= (const ZRGB<T1>& x, const ZRGB<T2>& y) 
{ return x.r >= y.r && x.g >= y.g && x.b >= y.b ; }

*/
template <class T> struct RGBType	{ typedef ZRGB<T> type; };

#ifndef _MSC_VER

template <class T> struct FloatType< ZRGB<T> >	{ typedef ZRGB<float> type; };
template <class T> struct ShortType< ZRGB<T> >	{ typedef ZRGB<short> type; };
template <class T> struct BYTEType< ZRGB<T> >	{ typedef ZRGB<BYTE> type; };
template <class T> struct GrayType< ZRGB<T> >	{ typedef T type; };
template <class T> struct RGBType< ZRGB<T> >	{ typedef ZRGB<T> type; };

#else

template <> struct FloatType< ZRGB<BYTE> >	{ typedef ZRGB<float>	type; };
template <> struct FloatType< ZRGB<short> > { typedef ZRGB<float>	type; };
template <> struct FloatType< ZRGB<WORD> >	{ typedef ZRGB<float>	type; };
template <> struct FloatType< ZRGB<int> >	{ typedef ZRGB<float>	type; };
template <> struct FloatType< ZRGB<UINT> >	{ typedef ZRGB<float>	type; };
template <> struct FloatType< ZRGB<float> > { typedef ZRGB<float>	type; };
template <> struct FloatType< ZRGB<double> >{ typedef ZRGB<float>	type; };

template <> struct BYTEType< ZRGB<BYTE> >	{ typedef ZRGB<BYTE>	type; };
template <> struct BYTEType< ZRGB<short> >	{ typedef ZRGB<BYTE>	type; };
template <> struct BYTEType< ZRGB<WORD> >	{ typedef ZRGB<BYTE>	type; };
template <> struct BYTEType< ZRGB<int> >	{ typedef ZRGB<BYTE>	type; };
template <> struct BYTEType< ZRGB<UINT> >	{ typedef ZRGB<BYTE>	type; };
template <> struct BYTEType< ZRGB<float> >	{ typedef ZRGB<BYTE>	type; };
template <> struct BYTEType< ZRGB<double> > { typedef ZRGB<BYTE>	type; };

template <> struct ShortType< ZRGB<BYTE> >	{ typedef ZRGB<short>	type; };
template <> struct ShortType< ZRGB<short> >	{ typedef ZRGB<short>	type; };
template <> struct ShortType< ZRGB<WORD> >	{ typedef ZRGB<short>	type; };
template <> struct ShortType< ZRGB<int> >	{ typedef ZRGB<short>	type; };
template <> struct ShortType< ZRGB<UINT> >	{ typedef ZRGB<short>	type; };
template <> struct ShortType< ZRGB<float> >	{ typedef ZRGB<short>	type; };
template <> struct ShortType< ZRGB<double> > { typedef ZRGB<short>	type; };

template <> struct GrayType< ZRGB<BYTE> >	{ typedef BYTE 	type; };
template <> struct GrayType< ZRGB<short> >	{ typedef short	type; };
template <> struct GrayType< ZRGB<WORD> >	{ typedef WORD 	type; };
template <> struct GrayType< ZRGB<int> >	{ typedef int	type; };
template <> struct GrayType< ZRGB<UINT> >	{ typedef UINT 	type; };
template <> struct GrayType< ZRGB<float> >	{ typedef float	type; };
template <> struct GrayType< ZRGB<double> > { typedef double type; };

template <> struct RGBType< ZRGB<BYTE> >	{ typedef ZRGB<BYTE>	type; };
template <> struct RGBType< ZRGB<short> >	{ typedef ZRGB<short>	type; };
template <> struct RGBType< ZRGB<WORD> >	{ typedef ZRGB<WORD>	type; };
template <> struct RGBType< ZRGB<int> >		{ typedef ZRGB<int>		type; };
template <> struct RGBType< ZRGB<UINT> >	{ typedef ZRGB<UINT>	type; };
template <> struct RGBType< ZRGB<float> >	{ typedef ZRGB<float>	type; };
template <> struct RGBType< ZRGB<double> >	{ typedef ZRGB<double>	type; };

#endif


#ifndef _MSC_VER

template <class T> class Opt< ZRGB<T> >
{
public:
    static double min, max;

	static ZRGB<T> Cast(const ZRGB<float>& v)
	{
		T r, g, b;

		if(v.r > max) r = T(max);
		else if(v.r < min) r = T(min);
		else r = T(v.r);

		if(v.g > max) g = T(max);
		else if(v.g < min) g = T(min);
		else g = T(v.g);

		if(v.b > max) b = T(max);
		else if(v.b < min) b = T(min);
		else b = T(v.b);
		
		return ZRGB<T>(r, g, b);
	}

    static ZRGB<T> Add(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r + y.r;
        float g = x.g + y.g;
        float b = x.b + y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Sub(const ZRGB<float>& x, const ZRGB<float>& y) { return Add(x, -y); }

    static ZRGB<T> Dif(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r - y.r;
        float g = x.g - y.g;
        float b = x.b - y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(Abs(r), Abs(g), Abs(b)));
    }

    static ZRGB<T> Mul(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r * y.r;
        float g = x.g * y.g;
        float b = x.b * y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Div(const ZRGB<float>& x, const ZRGB<float>& y) 
	{ return Mul(x, ZRGB<float>(1.0f/y.r, 1.0f/y.g, 1.0f/y.b)); }
};

template <class T> double Opt< ZRGB<T> >::min = MinValue<T>();
template <class T> double Opt< ZRGB<T> >::max = MaxValue<T>();

#else

template <> class Opt< ZRGB<BYTE> >
{
	typedef BYTE T;
public:
	static double min, max;

	static ZRGB<T> Cast(const ZRGB<float>& v)
	{
		T r, g, b;

		if(v.r > max) r = T(max);
		else if(v.r < min) r = T(min);
		else r = T(v.r);

		if(v.g > max) g = T(max);
		else if(v.g < min) g = T(min);
		else g = T(v.g);

		if(v.b > max) b = T(max);
		else if(v.b < min) b = T(min);
		else b = T(v.b);
		
		return ZRGB<T>(r, g, b);
	}

    static ZRGB<T> Add(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r + y.r;
        float g = x.g + y.g;
        float b = x.b + y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Sub(const ZRGB<float>& x, const ZRGB<float>& y) { return Add(x, -y); }

    static ZRGB<T> Dif(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r - y.r;
        float g = x.g - y.g;
        float b = x.b - y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(Abs(r), Abs(g), Abs(b)));
    }

    static ZRGB<T> Mul(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r * y.r;
        float g = x.g * y.g;
        float b = x.b * y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Div(const ZRGB<float>& x, const ZRGB<float>& y) 
	{ return Mul(x, ZRGB<float>(1.0f/y.r, 1.0f/y.g, 1.0f/y.b)); }
};

template <> class Opt< ZRGB<short> >
{
	typedef short T;
public:
	static double min, max;

	static ZRGB<T> Cast(const ZRGB<float>& v)
	{
		T r, g, b;

		if(v.r > max) r = T(max);
		else if(v.r < min) r = T(min);
		else r = T(v.r);

		if(v.g > max) g = T(max);
		else if(v.g < min) g = T(min);
		else g = T(v.g);

		if(v.b > max) b = T(max);
		else if(v.b < min) b = T(min);
		else b = T(v.b);
		
		return ZRGB<T>(r, g, b);
	}

    static ZRGB<T> Add(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r + y.r;
        float g = x.g + y.g;
        float b = x.b + y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Sub(const ZRGB<float>& x, const ZRGB<float>& y) { return Add(x, -y); }

    static ZRGB<T> Dif(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r - y.r;
        float g = x.g - y.g;
        float b = x.b - y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(Abs(r), Abs(g), Abs(b)));
    }

    static ZRGB<T> Mul(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r * y.r;
        float g = x.g * y.g;
        float b = x.b * y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Div(const ZRGB<float>& x, const ZRGB<float>& y) 
	{ return Mul(x, ZRGB<float>(1.0f/y.r, 1.0f/y.g, 1.0f/y.b)); }
};

template <> class Opt< ZRGB<float> >
{
	typedef float T;
public:
	static double min, max;

	static ZRGB<T> Cast(const ZRGB<float>& v)
	{
		T r, g, b;

		if(v.r > max) r = T(max);
		else if(v.r < min) r = T(min);
		else r = T(v.r);

		if(v.g > max) g = T(max);
		else if(v.g < min) g = T(min);
		else g = T(v.g);

		if(v.b > max) b = T(max);
		else if(v.b < min) b = T(min);
		else b = T(v.b);
		
		return ZRGB<T>(r, g, b);
	}

    static ZRGB<T> Add(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r + y.r;
        float g = x.g + y.g;
        float b = x.b + y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Sub(const ZRGB<float>& x, const ZRGB<float>& y) { return Add(x, -y); }

    static ZRGB<T> Dif(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r - y.r;
        float g = x.g - y.g;
        float b = x.b - y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(Abs(r), Abs(g), Abs(b)));
    }

    static ZRGB<T> Mul(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r * y.r;
        float g = x.g * y.g;
        float b = x.b * y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Div(const ZRGB<float>& x, const ZRGB<float>& y) 
	{ return Mul(x, ZRGB<float>(1.0f/y.r, 1.0f/y.g, 1.0f/y.b)); }
};

template <> class Opt< ZRGB<int> >
{
	typedef int T;
public:
	static double min, max;

	static ZRGB<T> Cast(const ZRGB<float>& v)
	{
		T r, g, b;

		if(v.r > max) r = T(max);
		else if(v.r < min) r = T(min);
		else r = T(v.r);

		if(v.g > max) g = T(max);
		else if(v.g < min) g = T(min);
		else g = T(v.g);

		if(v.b > max) b = T(max);
		else if(v.b < min) b = T(min);
		else b = T(v.b);
		
		return ZRGB<T>(r, g, b);
	}

    static ZRGB<T> Add(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r + y.r;
        float g = x.g + y.g;
        float b = x.b + y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Sub(const ZRGB<float>& x, const ZRGB<float>& y) { return Add(x, -y); }

    static ZRGB<T> Dif(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r - y.r;
        float g = x.g - y.g;
        float b = x.b - y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(Abs(r), Abs(g), Abs(b)));
    }

    static ZRGB<T> Mul(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r * y.r;
        float g = x.g * y.g;
        float b = x.b * y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Div(const ZRGB<float>& x, const ZRGB<float>& y) 
	{ return Mul(x, ZRGB<float>(1.0f/y.r, 1.0f/y.g, 1.0f/y.b)); }
};


template <> class Opt< ZRGB<WORD> >
{
	typedef WORD T;
public:
	static double min, max;

	static ZRGB<T> Cast(const ZRGB<float>& v)
	{
		T r, g, b;

		if(v.r > max) r = T(max);
		else if(v.r < min) r = T(min);
		else r = T(v.r);

		if(v.g > max) g = T(max);
		else if(v.g < min) g = T(min);
		else g = T(v.g);

		if(v.b > max) b = T(max);
		else if(v.b < min) b = T(min);
		else b = T(v.b);
		
		return ZRGB<T>(r, g, b);
	}

    static ZRGB<T> Add(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r + y.r;
        float g = x.g + y.g;
        float b = x.b + y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Sub(const ZRGB<float>& x, const ZRGB<float>& y) { return Add(x, -y); }

    static ZRGB<T> Dif(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r - y.r;
        float g = x.g - y.g;
        float b = x.b - y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(Abs(r), Abs(g), Abs(b)));
    }

    static ZRGB<T> Mul(const ZRGB<float>& x, const ZRGB<float>& y)
    {
        float r = x.r * y.r;
        float g = x.g * y.g;
        float b = x.b * y.b;

		return Opt< ZRGB<T> >::Cast(ZRGB<float>(r, g, b));
    }

    static ZRGB<T> Div(const ZRGB<float>& x, const ZRGB<float>& y) 
	{ return Mul(x, ZRGB<float>(1.0f/y.r, 1.0f/y.g, 1.0f/y.b)); }
};

#endif //_MSC_VER

#endif // ___RGBPIXELS_H__
