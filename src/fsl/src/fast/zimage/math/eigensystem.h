/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// The file declares the ZEigenSystem class.

#ifndef ___EIGENSYSTEM_H__
#define ___EIGENSYSTEM_H__

#include <iostream>  // C++ I/O
#include <iomanip>   // setprecision(),...
#include <cassert>    // assert()
#include <vector>

#include "mydefine.h"
#include "vector.h"
#include "matrix.h"

// template class ZEigenSystem
class ZEigenSystem
{
protected:
	UINT							m_nBase;
	UINT							m_nSizeofBase;

public:
	// data members
	ZVector<float>					m_MeanVector;
	ZVector<float>					m_EigenValue;
	ZMatrix<float>					m_EigenVector;	
	ZMatrix<float>					m_EigenVectorTrans;		// transpose of the real eigenvectors.

	// constructors
	ZEigenSystem (UINT nbase=0, UINT nsize=0) : m_nBase(nbase), m_nSizeofBase(nsize) {}

	// member functions
	bool Create (ZVector<float>* vectors, UINT number);
	bool Create (std::vector< ZVector<float> >& vectors);
	void Create (UINT nbase, UINT nsize);
	bool Projection(ZVector<float>& vec, ZVector<float>& coefficients) const;
	void Reconstruction(ZVector<float>& coefficients, ZVector<float>& recons) const;

	void Clear(void) {m_MeanVector.clear(); m_EigenVector.clear(); m_EigenVectorTrans.clear();} 

	bool Load(const char* fname);
	bool Save(const char* fname) const;

	UINT GetSize (void) const;
};

inline UINT ZEigenSystem::GetSize (void) const
{
	return m_nBase;
}

bool Eigensystem(const ZMatrix<float>& vectors, ZVector<float>& mean, ZMatrix<float>& eigenvectors, ZVector<float>& eigenvalues);

#endif // ___EIGENSYSTEM_H__
