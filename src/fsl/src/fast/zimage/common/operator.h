/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares template class for typed arithmatic

#ifndef __OPERATOR_H__
#define __OPERATOR_H__

#include "mydefine.h"

template <class T> struct FloatType	{ typedef float type; };
template <class T> struct ShortType	{ typedef short type; };
template <class T> struct BYTEType	{ typedef unsigned char	type; };
template <class T> struct GrayType	{ typedef T type; };

template <class T> class Opt
{
public:
	static double min, max;

	static T Cast(const float& v)
	{
		T ret;
		if(v > max) ret = T(max);
		else if(v < min) ret = T(min);
		else ret = T(v);
		return ret;
	}

	static T Add(const float& x, const float& y)
	{
		float add = x + y;
		T ret;
		if(add > max) ret = T(max);
		else if(add < min) ret = T(min);
		else ret = T(add);
		return ret;
	}

	static T Sub(const float& x, const float& y) { return Add(x, -y); }

	static T Dif(const float& x, const float& y)
	{
		float add = x - y;
		add = Abs(add);
		T ret;
		if(add > max) ret = T(max);
		else ret = T(add);
		return ret;
	}

	static T Mul(const float& x, const float& y)
	{
		float mul = x * y;
		T ret;
		if(mul > max) ret = T(max);
		else if(mul < min) ret = T(min);
		else ret = T(mul);
		return ret;
	}

	static T Div(const float& x, const float& y) { return Mul(x, 1.0f/y); }
};

#ifndef _MSC_VER

template <class T> double Opt<T>::min = MinValue<T>();
template <class T> double Opt<T>::max = MaxValue<T>();

#endif

template <class T1, class T2>
struct Add
{
	T1 operator () (const T1& x, const T2& y) const { return Opt<T1>::Add(x, y); }
};

template <class T1, class T2>
struct Sub
{
	T1 operator () (const T1& x, const T2& y) const { return Opt<T1>::Sub(x, y); }
};

template <class T>
struct Dif
{
	T operator () (const T& x, const T& y) const { return Opt<T>::Dif(x, y); }
};

template <class T1, class T2>
struct Mul
{
	T1 operator () (const T1& x, const T2& y) const { return Opt<T1>::Mul(x, y); }
};

template <class T1, class T2>
struct Div
{
	T1 operator () (const T1& x, const T2& y) const { return Opt<T1>::Div(x, y); }
};

template <class T>
struct Avg
{
	T operator () (const T& x, const T& y) const { return T((typename FloatType<T>::type(x)+y) / float(2)); }
};

template <class T>
struct AvgNoBackground
{
	T operator () (const T& x, const T& y) const { if(x==T(0)) return y; if(y==T(0)) return x; return T((typename FloatType<T>::type(x)+y) / float(2)); }
};

#endif
