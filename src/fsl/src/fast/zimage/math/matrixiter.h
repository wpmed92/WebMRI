/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares/defines the matrix iterator class.

#ifndef __ZMATRIXITERATOR_H__
#define __ZMATRIXITERATOR_H__

#include "vector.h"

template<class matrix_reference, class T, class row_iterator, class pointer=T*, class reference=T&> 
class ZMatrixIterator
{
	matrix_reference m_Matrix;
	int		m_nWidth;
	row_iterator m_pRow, m_pEnd;

	pointer	m_pCur;
	pointer	m_pRowEnd;

public:
	typedef ZMatrixIterator<matrix_reference, T, row_iterator, pointer, reference> matrix_iterator;

	ZMatrixIterator(matrix_reference matrix)
		: m_Matrix(matrix), m_nWidth(matrix.NCols()), m_pRow(matrix.begin()), m_pEnd(matrix.end()),
		  m_pCur(matrix.begin()->begin()), m_pRowEnd(m_pCur + m_nWidth)
	{
	}

	void reset(void)
	{
		m_pRow = m_Matrix.begin();
		m_pCur = m_pRow->begin();
		m_pRowEnd = m_pRow->end();
		m_pEnd = m_Matrix.end();
	}

	bool more(void) const
	{
		return (m_pRow < m_pEnd);
	}

	void operator ++ ()
	{
		if ( ++m_pCur < m_pRowEnd) return;
		else
		{
			m_pRow ++;
			m_pCur = m_pRow->begin();
			m_pRowEnd = m_pRow->end();
		}
	}

	void operator ++ (int) { operator ++(); }

/*  Both m_pRowLim and dif are undeclared and I have no idea where they come from.
	void operator += (int off) 
	{
		pointer n = m_pCur + off;
		if(n < m_pRowLim) { m_pCur = n; return; }
		
		m_pRow += dif % m_nWidth;
		int ncol = dif / m_nWidth;

		if(m_pRow < m_pEnd)
		{
			m_pCur = m_pRow->begin()+ncol;
			m_pRowEnd = m_pRow->end();
		}
	}
*/
	reference operator *() { return *m_pCur; }

	pointer operator ->() { return m_pCur; }

	operator pointer () { return m_pCur; }

	friend bool operator == (const matrix_iterator& x, const matrix_iterator& y)
	{
		return x.m_pCur == y.m_pCur;
	}

	friend bool operator != (const matrix_iterator& x, const matrix_iterator& y)
	{
		return x.m_pCur != y.m_pCur;
	}
};

#endif // __ZMATRIXITERATOR_H__
