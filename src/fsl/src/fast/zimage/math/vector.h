/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares the general std::vector-based ZVector class

#ifndef __VECTOR_H_
#define __VECTOR_H_

#include <vector>
#include <iostream> // C++ I/O
#include <iomanip>   // setprecision(),...
#include <cmath>

#include "mydefine.h"
#include "vectorD.h"

#undef min
#undef max

// template class ZVector
#undef _VALOP
#define _VALOP(TYPE, LENGTH, RHS) \
	ZVector<TYPE> s(LENGTH); \
	for (UINT i = 0; i < LENGTH; ++i) \
		s[i] = TYPE(RHS); \
	return (s)

#undef _VALGOP1
#define _VALGOP1(RHS) \
	for (UINT i = 0; i < this->size(); ++i) \
		this->begin()[i] = T (this->begin()[i] RHS); \
	return (*this)

template <class T>
class ZVector : public std::vector<T>
{
public:
	typedef T*			iterator;
	typedef const T*	const_iterator;

	// constructors
 	ZVector (UINT size=0, T v = T()) : std::vector<T>(size, v) {}

	template <class T2>
	ZVector (const ZVector<T2>& v) : std::vector<T>(v.size()) 
	{	
		for(UINT i=0; i<this->size(); i++) this->begin()[i] = T(v[i]);
	}
	
	iterator pbegin(){ 
	  return &(*(this->begin()));
	}

	iterator pend(){
	  return &(*(this->end()));
	}

	const_iterator pbegin() const{ 
	  return &(*(this->begin()));
	}

	const_iterator pend() const{
	  return &(*(this->end()));
	}

	template <class T2>
	ZVector& operator = (const ZVector<T2>& v) 
	{	
		resize(v.size());
		for(UINT i=0; i<this->size(); i++) this->begin()[i] = T(v[i]);
		return *this;
	}

	ZVector (const T * b, UINT n) : std::vector<T>(n) 
	{ copy(&b[0], &b[n], this->begin()); }

	ZVector (const ZVector2D& v) : std::vector<T>(2) { this->begin()[0]=v.x, this->begin()[1]=v.y; }
	ZVector (const ZVector3D& v) : std::vector<T>(3) { this->begin()[0]=v.x, this->begin()[1]=v.y, this->begin()[2]=v.z; }
	ZVector (const ZVector4D& v) : std::vector<T>(4) { this->begin()[0]=v.x, this->begin()[1]=v.y, this->begin()[2]=v.z, this->begin()[3]=v.w; }

	ZVector& replace(const ZVector& m, UINT start=0)
	{
		UINT end = Min(start + m.size(), this->size());
		for(UINT i = start; i < end; i++) this->begin()[i] = m[i-start];
		return *this;
	}

	ZVector& operator = (const ZVector2D& v) { *this = ZVector(v); }
	ZVector& operator = (const ZVector3D& v) { *this = ZVector(v); }
	ZVector& operator = (const ZVector4D& v) { *this = ZVector(v); }

	// member functions 
	ZVector& normalize() { *this /= len(*this); return *this; }

	ZVector& operator = (T num) { return fill(num); }
	ZVector& fill(T num){ for(iterator p=pbegin(); p!=pend(); p++) *p = num; return *this; }

	ZVector shift(int n) const
	{
		_VALOP(T, this->size(), 
			(0 < n && int(i) >= n) || (n < 0 && i-n<this->size()) ? this->begin()[i - n] : T(0)); 
	}

	ZVector cshift(int n) const
	{
		if (this->size() != 0)
			if (n < 0)
			{ if ((n += this->size()) < 0) n = this->size() - (-n) % this->size(); }
			else if (int(this->size()) <= n) n %= this->size();

		_VALOP(T, this->size(),
			int(i) >= n ? this->begin()[i - n] : this->begin()[i + this->size() - n]); 
	}

	ZVector apply(T func(T)) const	{_VALOP(T, this->size(), func(this->begin()[i])); }
	ZVector apply(T func(const T&)) const {_VALOP(T, this->size(), func(this->begin()[i])); }

	template <class T2>
	ZVector& operator += (const ZVector<T2>& v)	{ _VALGOP1(+ v[i]); }
	template <class T2>
	ZVector& operator -= (const ZVector<T2>& v)	{ _VALGOP1(- v[i]); }
	template <class T2>
	ZVector& operator *= (const ZVector<T2>& v)	{ _VALGOP1(* v[i]); }
	template <class T2>
	ZVector& operator /= (const ZVector<T2>& v)	{ _VALGOP1(/ v[i]); }

	ZVector& operator += (T v)	{_VALGOP1(+ v); }
	ZVector& operator -= (T v)	{_VALGOP1(- v); }
	ZVector& operator *= (float v)	{_VALGOP1(* v); }
	ZVector& operator /= (float v)	{_VALGOP1(/ v); }
}; // class ZVector

template<class T> inline float sqrlen(const ZVector<T>& v) { return dot(v, v); }
template<class T> inline double len(const ZVector<T>& v) { return sqrt(sqrlen(v)); }
template<class T> inline ZVector<T> norm(const ZVector<T>& v) { return v / len(v); }
template<class T> inline ZVector<T> sub(const ZVector<T>& v, UINT start, UINT len)
{
	if(len > v.size() - start) 
	{	len = v.size() - start; ZWarning("SubVector", "len %d is too large!", len); }
	return(ZVector<T>(v.begin()+start, len));
}
template<class T> inline T dot(const ZVector<T>& l, const ZVector<T>& r)
{ 
	T v = 0; 
	for(typename ZVector<T>::const_iterator p1=l.pbegin(), p2=r.pbegin(); p1!=l.pend(); p1++, p2++) v += *p1 * *p2;
	return v;
}
template<class T> inline T sum(const ZVector<T>& v)
{
	T s = 0;
	for(typename ZVector<T>::const_iterator p=v.pbegin(); p!=v.pend(); p++) s += *p;
	return s; 
}
template<class T> inline T Min(const ZVector<T>& v) { return *(std::min_element(v.begin(), v.end())); }
template<class T> inline T Max(const ZVector<T>& v) { return *(std::max_element(v.begin(), v.end())); }
template<class T> inline float mean(const ZVector<T>& v) { return float(sum(v))/v.size(); }
template<class T> inline ZVector<T> Abs(const ZVector<T>& v) 
{ 
	ZVector<T> r(v.size()); 
	typename ZVector<T>::const_iterator p_v = v.pbegin(); 
	typename ZVector<T>::iterator p_r = r.pbegin(); 
	for(; p_v!=v.pend(); p_v++, p_r++) *p_r = Abs(*p_v); 
	return r;
}
template<class T> inline float sqrdis(const ZVector<T>& l, const ZVector<T>& r) 
{ 
	T d=0; typename ZVector<T>::const_iterator p_l = l.pbegin(), p_r = r.pbegin();
	for(; p_l!=l.pend(); p_l++, p_r++) d += (*p_l-*p_r) * (*p_l-*p_r);
	return d;
}
template<class T> inline float dis(const ZVector<T>& l, const ZVector<T>& r) { return sqrt(sqrdis(l, r)); }

template<class T> inline
	ZVector<T> operator+(const ZVector<T>& v) { return v; }
template<class T> inline
	ZVector<T> operator-(const ZVector<T>& v) {_VALOP(T, v.size(), -v[i]); }
template <class T1, class T2>
inline ZVector<T1> operator+(const ZVector<T1>& l, const ZVector<T2>& r) { ZVector<T1> v(l); return v+=r; }
template <class T1, class T2>
inline ZVector<T1> operator-(const ZVector<T1>& l, const ZVector<T2>& r) { ZVector<T1> v(l); return v-=r; }
template <class T1, class T2>
inline ZVector<T1> operator*(const ZVector<T1>& l, const ZVector<T2>& r) { ZVector<T1> v(l); return v*=r; }
template <class T1, class T2>
inline ZVector<T1> operator/(const ZVector<T1>& l, const ZVector<T2>& r) { ZVector<T1> v(l); return v/=r; }
template<class T> inline
	ZVector<T> operator+(const ZVector<T>& l, T r) { ZVector<T> v(l); return v+=r; }
template<class T> inline
	ZVector<T> operator-(const ZVector<T>& l, T r) { ZVector<T> v(l); return v-=r; }
template<class T> inline
	ZVector<T> operator*(const ZVector<T>& l, float r) { ZVector<T> v(l); return v*=r; }
template<class T> inline
	ZVector<T> operator/(const ZVector<T>& l, float r) { ZVector<T> v(l); return v/=r; }
template<class T> inline
	ZVector<double> acos(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::acos(l[i])); }
template<class T> inline
	ZVector<double> asin(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::asin(l[i])); }
template<class T> inline
	ZVector<double> atan(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::atan(l[i])); }
template<class T> inline
	ZVector<double> atan2(const ZVector<T>& l,
		const ZVector<T>& r)
	{_VALOP(double, l.size(), ::atan2(l[i], r[i])); }
template<class T> inline
	ZVector<double> atan2(const ZVector<T>& l, const T& r)
	{_VALOP(double, l.size(), ::atan2(l[i], r)); }
template<class T> inline
	ZVector<double> atan2(const T& l, const ZVector<T>& r)
	{_VALOP(double, r.size(), ::atan2(l, r[i])); }
template<class T> inline
	ZVector<double> cos(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::cos(l[i])); }
template<class T> inline
	ZVector<double> cosh(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::cosh(l[i])); }
template<class T> inline
	ZVector<double> exp(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::exp(l[i])); }
template<class T> inline
	ZVector<double> log(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::log(l[i])); }
template<class T> inline
	ZVector<double> log10(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::log10(l[i])); }
template<class T> inline
	ZVector<double> pow(const ZVector<T>& l, const ZVector<T>& r)
	{_VALOP(double, l.size(), ::pow(l[i], r[i])); }
template<class T> inline
	ZVector<double> pow(const ZVector<T>& l, T r)
	{_VALOP(double, l.size(), ::pow(l[i], r)); }
template<class T> inline
	ZVector<double> pow(T l, const ZVector<T>& r)
	{_VALOP(double, r.size(), ::pow(l, r[i])); }
template<class T> inline
	ZVector<double> sin(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::sin(l[i])); }
template<class T> inline
	ZVector<double> sinh(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::sinh(l[i])); }
template<class T> inline
	ZVector<double> sqrt(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::sqrt(l[i])); }
template<class T> inline
	ZVector<double> tan(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::tan(l[i])); }
template<class T> inline
	ZVector<double> tanh(const ZVector<T>& l)
	{_VALOP(double, l.size(), ::tanh(l[i])); }


// overloaded insertion operator <<
template <class T>
inline std::ostream& operator << (std::ostream& s, const ZVector<T>& v)
{
	s << v[0]; 
	for(UINT i=1; i<v.size(); i++) s << " " << v[i];

	return s;
}

#endif // __VECTOR_H_
