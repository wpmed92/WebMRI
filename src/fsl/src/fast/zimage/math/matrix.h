/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares the general std::vector-based ZMatrix class

#ifndef __MATRIX_H_  // prevent multiple includes
#define __MATRIX_H_
#include <cstdio>
#include <iomanip>

#include "vector.h"
#include "matrixiter.h"
#include "common.h"
#include "matrix2D.h"
#include "matrix3D.h"
#include "matrix4D.h"

#undef _VALOP
#define _VALOP(TYPE, NROW, NCOL, RHS) \
	ZMatrix<TYPE> mat(NROW, NCOL); \
	typename ZMatrix<TYPE>::row_iterator P_mat = mat.pbegin(); \
	for (; P_mat != mat.pend(); P_mat++) *P_mat = RHS; \
	return (mat)

#undef _VALGOP1
#define _VALGOP1(RHS) \
	for (row_iterator ptr=pbegin(); ptr!=pend(); ptr++) *ptr RHS; \
	return (*this)

enum V2M {row=1, column=2, dia=3};

template <class T>
class ZMatrix	:	public std::vector< ZVector<T> >
{
public:
	typedef ZVector<T>*			row_iterator;
	typedef const ZVector<T>*	const_row_iterator;

	typedef ZMatrixIterator<ZMatrix<T>&, T, row_iterator,typename ZVector<T>::iterator, T&> uiterator;
	typedef ZMatrixIterator<const ZMatrix<T>&, T, const_row_iterator, typename ZVector<T>::const_iterator, const T&> const_uiterator;

	// constructor
	ZMatrix (UINT r=0, UINT c=0, T val=T()) : std::vector< ZVector<T> >(r)
	{ for(row_iterator ptr=pbegin(); ptr!=pend(); ptr++) ptr->resize(c, val); }

	template <class T1>
	ZMatrix (const ZMatrix<T1>& m) : std::vector< ZVector<T> >(m.size())
	{ 
		row_iterator		p_s = pbegin();
		const ZVector<T1>*	p_m = m.pbegin();
		for(; p_s!=pend(); p_s++, p_m++) *p_s = *p_m; 
	}

	row_iterator pbegin(){
	  return &(*(this->begin()));
	}

	row_iterator pend(){
	  return &(*(this->end()));
	}

	const_row_iterator pbegin() const{
	  return &(*(this->begin()));
	}

	const_row_iterator pend() const{
	  return &(*(this->end()));
	}

	template <class T1>
	ZMatrix& operator = (const ZMatrix<T1>& m)
	{ 
		resize(m.NRows(), m.NCols());
		row_iterator		p_s = pbegin();
		const ZVector<T1>*	p_m = m.pbegin();
		for(; p_s!=pend(); p_s++, p_m++) *p_s = *p_m; 
		return *this;
	}

	template <class T1>
	ZMatrix(T1* data, int width, int height) : std::vector< ZVector<T> >(height)
	{
		for(row_iterator ptr = pbegin(); ptr!=pend(); ptr++)
		{
			ptr->resize(width);
			for(int i=0; i<width; i++, data++) (*ptr)[i] = T(*data);
		}
	}

	ZMatrix (const ZVector<T>& v, V2M dir = row);
	ZMatrix (const ZMatrix2D& m) : std::vector< ZVector<T> >() { Import(m.begin(), 2, 2); }
	ZMatrix (const ZMatrix3D& m) : std::vector< ZVector<T> >() { Import(m.begin(), 3, 3); }
	ZMatrix (const ZMatrix4D& m) : std::vector< ZVector<T> >() { Import(m.begin(), 4, 4); }

	template <class T1>
	void Import(T1* data, int width, int height)
	{
		resize(height, width);
		for(row_iterator ptr = pbegin(); ptr!=pend(); ptr++)
		for(int i=0; i<width; i++, data++) (*ptr)[i] = T(*data);
	}

	template <class T1>
	void Export(T1* data, int width, int height)
	{
		row_iterator ptr = pbegin();
		for(; ptr!=pend(); ptr++)
			for(int i=0; i<width; i++, data++) *data = Opt<T1>::Cast((*ptr)[i]);
	}

	UINT NRows () const { return this->size(); }
	UINT NCols () const { return this->begin()[0].size(); }
	UINT NElements () const { return NRows() * NCols(); }

	void RowStatistics (ZVector<double>& mean, ZMatrix<double>& covariance, bool bNoZero=false) const;
	void RowStatistics (ZVector<float>& mean, ZMatrix<float>& covariance, bool bNoZero=false) const
	{
		ZVector<double> m;
		ZMatrix<double> c;
		RowStatistics(m, c, bNoZero);
		mean = m; covariance = c;
	}

	T& operator() (UINT i, UINT j) { return this->begin()[i][j]; }
	const T& operator() (UINT i, UINT j) const { return this->begin()[i][j]; }

    void resize(UINT nrows, UINT ncols, T val=T())
	{
		if (this->size() < nrows) insert(this->end(), nrows - this->size(), ZVector<T>());
		else if (nrows < this->size()) erase(this->begin() + nrows, this->end()); 

		for(row_iterator ptr=pbegin(); ptr!=pend(); ptr++) ptr->resize(ncols, val);
	}

	void GetColumn(ZVector<T>& vec, UINT col) const;
	void SetColumn(const ZVector<T>& vec, UINT col);

	ZMatrix<T>&	normalize(V2M v=column);
	ZMatrix<T>&	Replace(const ZMatrix<T>& m, UINT top, UINT left);

	ZMatrix<float> ColumnScaler() const;

	ZMatrix<T>& operator = (T n) { fill(n); return *this; }
	void	fill(T n){ for(row_iterator ptr=pbegin(); ptr!=pend(); ptr++) ptr->fill(n); }

	ZMatrix& MakeDiag(T v) { int d=Min(NCols(), NRows()); for(int i=0; i<d; i++) this->begin()[i][i] = v; return *this; }
	template <class T2>
	ZMatrix& operator += (const ZMatrix<T2> & r)	{ const_row_iterator p_r=r.pbegin(); _VALGOP1(+= *p_r++); }
	template <class T2>
	ZMatrix& operator -= (const ZMatrix<T2> & r)	{ const_row_iterator p_r=r.pbegin(); _VALGOP1(-= *p_r++); }
	template <class T2>
	ZMatrix& operator *= (const ZMatrix<T2> & r)	{ const_row_iterator p_r=r.pbegin(); _VALGOP1(*= *p_r++); }
	template <class T2>
	ZMatrix& operator /= (const ZMatrix<T2> & r)	{ const_row_iterator p_r=r.pbegin(); _VALGOP1(/= *p_r++); }

	ZMatrix& operator += (T r)	{_VALGOP1(+= r); }
	ZMatrix& operator -= (T r)	{_VALGOP1(-= r); }
	ZMatrix& operator *= (float r)	{_VALGOP1(*= r); }
	ZMatrix& operator /= (float r)	{_VALGOP1(/= r); }

	template <class T2>
	ZMatrix operator * (const ZMatrix<T2>& m) const
	{
		if(NCols() != m.NRows()) { ZError("Matrix *", "Size unmatched!"); return ZMatrix(); }

		ZMatrix r (NRows(), m.NCols()) ;

		row_iterator p_r = r.pbegin(); const_row_iterator ptr = pbegin();
		for (; ptr!=pend(); ptr++, p_r++)
		for (UINT j=0; j<m.NCols(); j++)
		{
			T sum = 0; 	const ZVector<T2>* p_m = m.pbegin();
			for (UINT k=0; p_m!=m.pend(); k++, p_m++) sum += (*ptr)[k] * (*p_m)[j];
			(*p_r)[j] = sum;
		}
		
		return r;
	}

	template <class T2>
	ZVector<T> operator * (const ZVector<T2>& V) const
	{
		if(NCols() != V.size()) { ZError("Matrix *", "Size unmatched!"); return ZVector<T>(); }

		ZVector<T> r (NRows()); 
		typename ZVector<T>::iterator p_r = r.pbegin();
		const_row_iterator ptr = pbegin();
		for (; ptr!=pend(); ptr++, p_r++)
		{
			const T* p_v = V.pbegin();
			for(UINT i=0; p_v!=V.pend(); p_v++, i++) *p_r += (*ptr)[i] * *p_v;
		}

		return r;
	}
}; // class ZMatrix

template<class T> inline float det(const ZMatrix<T>& v)
{
	if(v.NRows() != v.NCols()) 
	{
		ZError("Matrix det", "The matrix is not square!");
		return 0;
	}
	ZMatrix<float> A = v;  bool sign = LUdcmp(A);
	float result=1; for(UINT i=0; i<v.NRows(); i++) result *= A[i][i];
	return sign ? result : -result;
}
template<class T> inline ZMatrix<T> sub(const ZMatrix<T>& v, UINT x, UINT y, UINT w, UINT h)
{
	if(w > v.NCols() - x) 
	{	w = v.NCols() - x; ZWarning("SubMatrix", "w %d is too large!", w); }
	if(h > v.NRows() - y) 
	{	h = v.NRows() - y; ZWarning("SubMatrix", "h %d is too large!", h); }
	ZMatrix<T> s(h, w); 
	typename ZMatrix<T>::row_iterator p_s=s.pbegin();
	typename ZMatrix<T>::const_row_iterator p_v=v.pbegin();
	for(; p_s!=s.pend(); p_s++, p_v++) copy(p_v->begin()+x, p_v->begin()+x+w, p_s->begin());
	return s;
}
template<class T> inline T sum(const ZMatrix<T>& v)
{
	T s = 0;
	for(typename ZMatrix<T>::const_row_iterator p=v.pbegin(); 
	    p!=v.pend(); p++) s += sum(*p);
	return s; 
}
template<class T> inline T Min(const ZMatrix<T>& v) 
{	
	T m = v(0,0);
	typename ZMatrix<T>::const_uiterator p_v(v); 
	for(; p_v.more(); p_v++) if(m > *p_v) m=*p_v;
	return m;
}
template<class T> inline T Max(const ZMatrix<T>& v) 
{	
	T m = v(0,0);
	typename ZMatrix<T>::const_uiterator p_v(v); 
	for(; p_v.more(); p_v++) if(m < *p_v) m=*p_v;
	return m;
}
template<class T> inline float mean(const ZMatrix<T>& v) { return float(sum(v))/v.NElements(); }
template<class T> inline ZVector<T> mean_row(const ZMatrix<T>& v) 
{ 
	ZVector<float> m(v.NCols());
	for(typename ZMatrix<T>::const_row_iterator ptr=v.pbegin(); 
	    ptr!=v.pend(); ptr++) m += *ptr;
	return m / v.NRows();
}
template<class T> inline ZVector<T> mean_col(const ZMatrix<T>& v) 
{ 
	ZVector<T> m(v.NRows());
	typename ZVector<T>::iterator p_m = m.pbegin(); 
	typename ZMatrix<T>::const_row_iterator p_v = v.pbegin();
	for(; p_m!=m.pend(); p_m++, p_v++) *p_m = mean(*p_v);
	return m;
}
template<class T> inline ZMatrix<T> Abs(const ZMatrix<T>& v) 
{ 
	ZMatrix<T> r(v.NRows(), v.NCols()); 
	typename ZMatrix<T>::const_uiterator p_v = v.pbegin(); 
	typename ZMatrix<T>::uiterator p_r = r.pbegin(); 
	for(; p_v.more(); p_v++, p_r++) *p_r = Abs(*p_v); 
	return r;
}
template<class T> inline ZMatrix<T> trans(const ZMatrix<T>& v)
{
	ZMatrix<T> r(v.NCols(), v.NRows());
	for(UINT i=0; i<v.NCols(); i++) for(UINT j=0; j<v.NRows(); j++) r(i, j) = v(j, i);
	return r;
}
template<class T> inline T tracs(const ZMatrix<T>& v)
{
	 int d=Min(v.NCols(), v.NRows()); T s=0;
	 for(int i=0; i<d; i++) s+=v(i, i); 
	 return s; 
}
template<class T> inline ZMatrix<float> cov(const ZMatrix<T>& v)
{
	ZMatrix<float> mm = v; 	ZVector<float> m = mean_row(v);
	for(typename ZMatrix<float>::row_iterator p_mm = mm.pbegin(); p_mm!=mm.pend(); p_mm++) *p_mm -= m;
	ZMatrix<float> cov = trans(mm) * mm;
	return (cov /= v.NRows());
}

template<class T> inline ZMatrix<T> operator+(const ZMatrix<T>& l) { return l; }
template<class T> inline ZMatrix<T> operator-(const ZMatrix<T>& l)
{ typename ZMatrix<T>::const_row_iterator p_l=l.pbegin(); _VALOP(T, l.NRows(), l.NCols(), -*p_l++); }
template <class T1, class T2>
inline ZMatrix<T1> operator + (const ZMatrix<T1>& l, const ZMatrix<T2>& r) { ZMatrix<T1> v(l); return v+=r; }
template <class T1, class T2>
inline ZMatrix<T1> operator - (const ZMatrix<T1>& l, const ZMatrix<T2>& r) { ZMatrix<T1> v(l); return v-=r; }
template <class T1, class T2>
inline ZMatrix<T1> operator / (const ZMatrix<T1>& l, const ZMatrix<T2>& r) { ZMatrix<T1> v(l); return v/=r; }

template<class T> inline ZMatrix<T> operator + (const ZMatrix<T>& l, T r)	{ ZMatrix<T> v(l); return v+=r; }
template<class T> inline ZMatrix<T> operator - (const ZMatrix<T>& l, T r)	{ ZMatrix<T> v(l); return v-=r; }
template<class T> inline ZMatrix<T> operator * (const ZMatrix<T>& l, float r)	{ ZMatrix<T> v(l); return v*=r; }
template<class T> inline ZMatrix<T> operator / (const ZMatrix<T>& l, float r)	{ ZMatrix<T> v(l); return v/=r; }

template<class T> inline
	ZMatrix<double> acos(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::acos(*P_X++)); }
template<class T> inline
	ZMatrix<double> asin(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::asin(*P_X++)); }
template<class T> inline
	ZMatrix<double> atan(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::atan(*P_X++)); }
template<class T> inline
	ZMatrix<double> atan2(const ZMatrix<T>& x, const ZMatrix<T>& y)
	{ 
		typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(), P_Y=y.pbegin(); 
		_VALOP(double, x.NRows(), x.NCols(), ::atan2(*P_X++, *P_Y++)); 
	}
template<class T> inline
	ZMatrix<double> atan2(const ZMatrix<T>& x, const T& y)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::atan2(*P_X++, y)); }
template<class T> inline
	ZMatrix<double> atan2(const T& x, const ZMatrix<T>& y)
	{ typename ZMatrix<T>::const_row_iterator P_Y=y.pbegin(); _VALOP(double, y.NRows(), y.NCols(), ::atan2(x, *P_Y++)); }
template<class T> inline
	ZMatrix<double> cos(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::cos(*P_X++)); }
template<class T> inline
	ZMatrix<double> cosh(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::cosh(*P_X++)); }
template<class T> inline
	ZMatrix<double> exp(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::exp(*P_X++)); }
template<class T> inline
	ZMatrix<double> log(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::log(*P_X++)); }
template<class T> inline
	ZMatrix<double> log10(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::log10(*P_X++)); }
template<class T> inline
	ZMatrix<double> pow(const ZMatrix<T>& x, const ZMatrix<T>& y)
	{ 
		typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(), P_Y=y.pbegin(); 
		_VALOP(double, x.NRows(), x.NCols(), ::pow(*P_X++, *P_Y++)); 
	}
template<class T> inline
	ZMatrix<double> pow(const ZMatrix<T>& x, const T& y)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::pow(*P_X++, y)); }
template<class T> inline
	ZMatrix<double> pow(const T& x, const ZMatrix<T>& y)
	{ typename ZMatrix<T>::const_row_iterator P_Y=y.begin(); _VALOP(double, y.NRows(), y.NCols(), ::pow(x, *P_Y++)); }
template<class T> inline
	ZMatrix<double> sin(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.begin(); _VALOP(double, x.NRows(), x.NCols(), ::sin(*P_X++)); }
template<class T> inline
	ZMatrix<double> sinh(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.begin(); _VALOP(double, x.NRows(), x.NCols(), ::sinh(*P_X++)); }
template<class T> inline
	ZMatrix<double> sqrt(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.begin(); _VALOP(double, x.NRows(), x.NCols(), ::sqrt(*P_X++)); }
template<class T> inline
	ZMatrix<double> tan(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::tan(*P_X++)); }
template<class T> inline
	ZMatrix<double> tanh(const ZMatrix<T>& x)
	{ typename ZMatrix<T>::const_row_iterator P_X=x.pbegin(); _VALOP(double, x.NRows(), x.NCols(), ::tanh(*P_X++)); }
template<class T1, class T2> inline
	ZVector<T1> operator*(const ZVector<T1>& x, const ZMatrix<T2>& y)
{
	ZVector<T1> mat(y.NCols());
	typename ZVector<T1>::iterator P_mat = mat.pbegin();
	for (UINT i=0; P_mat != mat.pend(); P_mat++, i++)
	{
		typename ZVector<T1>::const_iterator P_X = x.pbegin();
		const ZVector<T2>* P_Y = y.pbegin();
		for (; P_X != x.pend(); P_X++, P_Y++)
			*P_mat += *P_X * (*P_Y)[i];
	}
	return (mat);
}

#endif // _MATRIX_H
