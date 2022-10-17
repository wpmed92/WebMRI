/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <iostream>
#include "matrix.h"
#include "common.h"
#include "mathcommon.h"

using namespace std;

template <class T>
ZMatrix<T>::ZMatrix (const ZVector<T>& v, V2M dir)
{
	UINT nrow=1, ncol=1;

	if(dir & row) nrow = v.size();
	if(dir & column) ncol = v.size();

	resize(nrow, ncol);

	if(dir == row)	*(pbegin()) = v;
	else if(dir == column) 
	{
		typename ZVector<T>::const_iterator p_v = v.pbegin();
		row_iterator ptr = pbegin();
		for (; ptr!=pend(); ptr++, p_v++) (*ptr)[0] = *p_v;
	}
	else
	{
		typename ZVector<T>::const_iterator p_v = v.pbegin();
		row_iterator ptr = pbegin();
		for (UINT i=0; ptr!=pend(); ptr++, p_v++, i++) (*ptr)[i] = *p_v;
	}
}

template <class T>
void ZMatrix<T>::GetColumn(ZVector<T>& vec, UINT col) const
{
	if(col >= NCols()) { ZError("Matrix GetColumn", "Out of range!"); return; }
	if(vec.size()!=NRows()) vec.resize(NRows());
	typename ZVector<T>::iterator p_vec = vec.pbegin();
	const_row_iterator ptr=pbegin();
	for(; p_vec!=vec.pend(); p_vec++, ptr++) *p_vec = (*ptr)[col];
}

template <class T>
void ZMatrix<T>::SetColumn(const ZVector<T>& vec, UINT col)
{
	if(col >= NCols()) { ZError("Matrix GetColumn", "Out of range!"); return; }

	if(NRows() != vec.size()) { ZError("Matrix GetColumn", "Size unmatched!"); return; }
	typename ZVector<T>::const_iterator p_vec = vec.pbegin();
	row_iterator ptr=pbegin();
	for(; p_vec!=vec.pend(); p_vec++, ptr++) (*ptr)[col] = *p_vec;
}

template <class T> 
ZMatrix<float> ZMatrix<T>::ColumnScaler() const
{
	ZMatrix<float> r(NCols(), NCols());

	typename ZMatrix<float>::row_iterator p_r = r.pbegin();

	ZVector<T> vec(NRows());
	for(UINT i=0; p_r!=r.pend(); p_r++, i++)
	{
		GetColumn(vec, i); 	float norm = len(vec);
		(*p_r)[i] = (norm > 0) ? (*p_r)[i] = float(1.0f / norm) : 1;
  	}

	return r;
}

template <class T>
ZMatrix<T>& ZMatrix<T>::normalize(V2M dir)
{
	if(dir == row) for(row_iterator ptr=pbegin(); ptr!=pend(); ptr++) ptr->normalize();
	else if(dir == column)
	{
		float norm;  ZVector<T> vec(NRows());
		for(UINT i=0; i<NCols(); i++)
		{
			GetColumn(vec, i);
			if((norm = len(vec)) == 0) continue;
			for(UINT j=0; j<NRows(); j++) pbegin()[j][i] = T(float(pbegin()[j][i]) / norm);
  		}
	}
	return *this;
}

template <class T>
ZMatrix<T>& ZMatrix<T>::Replace(const ZMatrix<T>& m, UINT top, UINT left)
{
	UINT bottom = Min(top + m.NRows(), NRows());
	UINT right = Min(left + m.NCols(), NCols());
	for(UINT i = top; i < bottom; i++)
	for(UINT j = left; j < right; j++)
		this->begin()[i][j] = m[i-top][j-left];

	return *this;
}

template <class T>
void ZMatrix<T>::RowStatistics(ZVector<double>& mean, ZMatrix<double>& covariance, bool bNoZero) const
{
	int spec = NCols();
	mean.resize(spec); covariance.resize(spec, spec);
	mean = 0; covariance = 0; int nor=0;

	const_row_iterator p_r;
	for(p_r = pbegin(); p_r!=pend(); p_r++)
	if(!bNoZero || (*p_r)[0]) 
	{
		mean += *p_r;
		nor++;
	}

	if(nor == 0) return;

	mean /= nor;

	ZVector<double> dif(spec);
	for(p_r = pbegin(); p_r!=pend(); p_r++)
	if(!bNoZero || (*p_r)[0]) 
	{
		dif = *p_r; dif -= mean;

		for(int i=0; i<spec; i++)
		for(int j=0; j<spec; j++)
			covariance(i,j) += dif[i] * dif[j];
	}

	covariance /= nor;
}

template class ZMatrix<double>;
template class ZMatrix<float>;
template class ZMatrix<int>;
