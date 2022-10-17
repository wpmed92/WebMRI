/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// The realization of the ZEigenSystem class.

#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

#include "eigensystem.h"
#include "jacobi.h"
#include "matrix.h"

bool Eigensystem(const ZMatrix<float>& vectors, ZVector<float>& mean, ZMatrix<float>& eigenvectors, ZVector<float>& eigenvalues)
{
	if(vectors.NRows() == 0 || vectors.NCols() == 0) return false;

	int nrows = vectors.NRows();
	int ncols = vectors.NCols();
	ZMatrix<float> mm(nrows, ncols);
	
	mean.clear(); eigenvectors.clear(); eigenvalues.clear();
	mean.resize(ncols);

	int i;
	ZMatrix<float>::const_row_iterator p_vec = vectors.pbegin();
	ZMatrix<float>::row_iterator p_mm = mm.pbegin();
	for(i=0; i<nrows; i++, p_vec++, p_mm++) 
	for(int j=0; j<ncols; j++) mean[j] += (*p_mm)[j] = (*p_vec)[j];
	mean /= nrows;

	p_mm = mm.pbegin();
	for(i=0; i<nrows; i++, p_mm++)
	for(int j=0; j<ncols; j++) (*p_mm)[j] -= mean[j];
	ZMatrix<float> covariance = trans(mm) * mm / (nrows - 1);
	
	return Jacobi(covariance, eigenvalues, eigenvectors);
}

void ZEigenSystem::Create (UINT nbase, UINT nsize)
{
	m_nBase = nbase;
	m_nSizeofBase = nsize;

	m_MeanVector.resize(m_nSizeofBase);
	m_EigenValue.resize(m_nBase);
	m_EigenVectorTrans.resize(m_nBase, m_nSizeofBase);
	m_EigenVector.resize(m_nSizeofBase, m_nBase);
}

bool ZEigenSystem::Create(ZVector<float>* vectors, UINT number)
{
	vector< ZVector<float> > vec(vectors, vectors+number-1);
	return Create(vec);
}

bool ZEigenSystem::Create(vector< ZVector<float> >& vectors)
{
	UINT i, size = vectors[0].size();
	UINT number = vectors.size();
	UINT num = size > number ? number : size;
	for(i=0; i<number; i++) 
		if(size != vectors[i].size()) return false;

	ZMatrix<float>					Corvariance;

	cout << "Computing the mean and difference vectors..." << flush;
	m_nSizeofBase = size;
	m_nBase = number;
	m_MeanVector.resize(size);
	for(i=0; i<number; i++) m_MeanVector += vectors[i];
	m_MeanVector /= number;

	ZMatrix<float> vec_grp_trans(number, size), vec_grp;
	for(i=0; i<number; i++)
	{
		vec_grp_trans[i] = vectors[i] - m_MeanVector;
	}

	vec_grp_trans /= sqrt(float(number));

	vec_grp = trans(vec_grp_trans);
	cout << "done!" << endl;

	cout << "Computing the corvariance matrix...." << flush;
	if(size > number) Corvariance = vec_grp_trans * vec_grp;
	else Corvariance = vec_grp * vec_grp_trans;
	cout << "done!" << endl;

	cout << "Computing eigenvalue and eigenvectors using Jacobi method...." << flush;
	if(Jacobi(Corvariance, m_EigenValue, m_EigenVector) == false)
	{
		cout << "fail to compute eigen system!" << endl;
		return false;
	}
	cout << "done!" << endl;

	cout << "Sorting...." << flush;
	m_EigenVector = trans(m_EigenVector);

	ZVector<float> vec(num);
	for(i=0; i<num; i++) 
	{
		UINT k;
		float p = m_EigenValue[k=i];
		for(UINT j=i+1; j<number; j++) 
			if(m_EigenValue[j] > p) p = m_EigenValue[k=j];
		
		if(k != i)
		{
			m_EigenValue[k] = m_EigenValue[i];
			m_EigenValue[i] = p;

			vec = m_EigenVector[k];
			m_EigenVector[k] = m_EigenVector[i];
			m_EigenVector[i] = vec;
		}
	}
	cout << "done!" << endl;

	cout << "Creating the final eigen system...." << flush;
	m_EigenVector = trans(m_EigenVector);
	if(size > number)
	{
		ZMatrix<float> tmp = vec_grp * m_EigenVector;
		m_EigenVector = tmp;
	}

	m_EigenVector.normalize();
	m_EigenVectorTrans = trans(m_EigenVector);
	cout << "done successfully!" << endl;

	return true;
}

bool ZEigenSystem::Projection(ZVector<float>& vec, ZVector<float>& coefficients) const
{
	if(vec.size() != m_nSizeofBase) return false;

	coefficients.resize(m_nBase);
	coefficients = m_EigenVectorTrans * (vec - m_MeanVector);
	
	return true;
}

void ZEigenSystem::Reconstruction(ZVector<float>& coefficients, ZVector<float>& recons) const
{
	UINT size = (coefficients.size() > m_nBase) ? m_nBase : coefficients.size();
	recons.resize(m_nSizeofBase);

	for(UINT i=0; i<m_nSizeofBase; i++) 
	for(UINT j=0; j<size; j++)
		recons[i] += m_EigenVector[i][j] * coefficients[j];

	recons += m_MeanVector;
}

bool ZEigenSystem::Save(const char* fname) const
{
	fstream	file(fname, ios::out | ios::binary);
	if(!file)
	{
		cout << "Can not open file " << fname << " for writing!" << endl;
		return false;
	}

	file.write((char*)&m_nBase, sizeof(int));
	file.write((char*)&m_nSizeofBase, sizeof(int));

	file.write((char*)(m_MeanVector.pbegin()), m_nSizeofBase*sizeof(float));
	file.write((char*)(m_EigenValue.pbegin()), m_nBase*sizeof(float));
	for(UINT i=0; i<m_nBase; i++)
	{
		file.write((char*)(m_EigenVectorTrans[i].pbegin()), m_nSizeofBase*sizeof(float));
	}
	
	file.close();
	return true;
}

bool ZEigenSystem::Load(const char* fname)
{
	fstream	file(fname, ios::in | ios::binary);
	if(!file)
	{
		cout << "Can not open file " << fname << " for reading!" << endl;
		return false;
	}

	file.read((char*)&m_nBase, sizeof(int));
	file.read((char*)&m_nSizeofBase, sizeof(int));

	if(file.fail())
	{
		cout << "read file error!" << endl;
		return false;
	}

	m_MeanVector.resize(m_nSizeofBase);
	m_EigenValue.resize(m_nBase);
	m_EigenVectorTrans.resize(m_nBase, m_nSizeofBase);

	file.read((char*)(m_MeanVector.pbegin()), m_nSizeofBase*sizeof(float));
	file.read((char*)(m_EigenValue.pbegin()), m_nBase*sizeof(float));
	for(UINT i=0; i<m_nBase; i++)
	{
		file.read((char*)(m_EigenVectorTrans[i].pbegin()), m_nSizeofBase*sizeof(float));
	}
	m_EigenVector = trans(m_EigenVectorTrans);
	
	file.close();
	return true;
}
