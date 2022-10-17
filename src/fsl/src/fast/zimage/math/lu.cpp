/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "zmath.h"
#include "common.h"
using namespace std;

///////////////////////////////////////////////////////////////////////////
// LUdcmp function: Given a matrix A (N x N), this routine
//                  replaces it by the LU decomposition of a
//                  rowwise permutation of itself. a and n are
//                  input. a is output, arranged with L and U in
//                  the same matrix; indx[1..n] is an output vector
//                  that records the row permutation affected by
//                  the partial pivoting; d is output as +/-1
//                  depending on whether the number of row
//                  interchanges was even or odd, respectively.
//                  This routine is used in combination with LUbksb()
//                  to solve linear equations or invert a matrix.
///////////////////////////////////////////////////////////////////////////

bool LUdcmp(ZMatrix<double>& A, ZVector<int>& indx)
{
	UINT N = A.NRows();
	if(A.NCols() != N) { ZError("Matrix LUdcmp", "Size unmatched!"); return false; }

	UINT	i,imax(0),j,k;
	double	big,dum,sum;
	bool	sign = true;

	indx.resize(N);
	ZVector<double> vv(N);		//vv stores the implicit scaling of each row.

	for (i=0; i<N; i++) 
	{						//Loop over rows to get the implicit scaling information
		big=0.0;
		for (j=0; j<N; j++) 
			if (Abs(A[i][j]) > big) big = Abs(A[i][j]);

		if (big == 0) return false;	//No nonzero largest element.
		vv[i] = 1.0f / big;				//Save the scaling.
	}

	for (j=0; j<N; j++) 
	{				//This is the loop over nCols of Crout's method.
		for (i=0; i<j; i++) 
		{			//This is equation (2.3.12) except for i = j.
			sum = A[i][j];
			for (k=0; k<i; k++) sum -= A[i][k] * A[k][j];
			A[i][j]=sum;
		}
		big=0; //Initialize for the search for largest pivot element.
		for (i=j; i<N; i++) 
		{			//This is i = j of equation (2.3.12) and i = j +1: ::N of equation (2.3.13). 
			sum = A[i][j];
			for (k=0; k<j; k++) sum -= A[i][k] * A[k][j];
			A[i][j]=sum;
			if ( (dum = vv[i] * Abs(sum)) >= big) 
			{
				big = dum;
				imax = i;
			}
		}
		if (j != imax) 
		{							//Do we need to interchange rows?
			for (k=0; k<N; k++) 
			{						//Yes, do so...
				dum = A[imax][k];
				A[imax][k] = A[j][k];
				A[j][k] = dum;
			}
			sign = !sign;			//...and change the parity of d.
			vv[imax] = vv[j];	//Also interchange the scale factor.
		}
		indx[j]=imax;
		if (Abs(A[j][j]) < TINY) A[j][j] = 0;

// If the pivot element is zero the matrix is singular (at least to the precision of the
// algorithm). For some applications on singular matrices, it is desirable to substitute
// TINY for zero.
		if (j != N) 
		{			//Now, Manally, divide by the pivot element.
			dum = 1.0f / (A[j][j]);
			for (i=j+1; i<N; i++) A[i][j] *= dum;
		}
	} //Go back for the next column in the reduction.

	return sign;
}

///////////////////////////////////////////////////////////////////////////
//  LUbksb function: Solves the set of n linear equations A.X = B.
//                   Here A (N x N) is input, not as the matrix A
//                   but rather as its LU decomposition, determined
//                   by the routine LUdcmp(). indx[1..N] is input as
//                   the permutation vector returned by LUdcmp().
//                   b[1..n] is input as the right hand side vector B,
//                   and returns with the solution vector X. a, n, and
//                   indx are not modified by this routine and can be
//                   left in place for successive calls with different
//                   right-hand sides b. This routine takes into account
//                   the possibility that b will begin with many zero
//                   elements, so it it efficient for use in matrix
//                   inversion.
///////////////////////////////////////////////////////////////////////////

void LUbksb (const ZMatrix<double>& A, const ZVector<int>& indx, ZVector<double>& B )
{
	int N = int(A.NRows());
	if(A.NCols() != UINT(N)) { ZError("Matrix LUbksb", "Size unmatched!"); return; }

    int i, ii=-1;  double sum;

    for (i=0; i<N; i++)			// When ii is set to a positive value, it will 
	{
        int ip = indx[i];		// become the index of the first nonvanishing 
        double sum = B[ip];     // element of b. We now do the forward substitution. 
        B[ip] = B[i];			// The only new wrinkle is to unscramble the 
        if (ii>-1)					// permutation as we go. 
            for (int j=ii; j<=i-1; j++) sum -= A[i][j]*B[j];
        else if (sum)			// A nonzero element was encountered, so from now on 
            ii = i;				// we will have to do the sums in the loop above 
        
		B[i]=sum;
    }
    
	for (i=N-1; i>=0; i--)		// Now we do the backsubstitution. 
	{
        sum = B[i];
        for (int j=i+1; j<N; j++) sum -= A[i][j] * B[j];
        B[i] = sum / A[i][i];  // Store a component of the solution vector X. 
    }                      // All done! 
}

void LUbksb (const ZMatrix<float>& A, const ZVector<int>& indx, ZVector<float>& B )
{
	ZMatrix<double> a = A;
	ZVector<double> b = B;
	LUbksb(a, indx, b);
	B = b;
}

bool LUdcmp(ZMatrix<float>& A, ZVector<int>& indx)
{
	ZMatrix<double> a = A;
	bool ret = LUdcmp(a, indx);
	A = a;
	return ret;
}

bool LUdcmp(ZMatrix<float>& A)
{
	ZVector<int> indx;
	return LUdcmp(A, indx);
}

bool LUdcmp(ZMatrix<double>& A)
{
	ZVector<int> indx;
	return LUdcmp(A, indx);
}
