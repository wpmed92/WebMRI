/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "matrix.h"
#include "common.h"
using namespace std;

static double PYTHAG(double a, double b)
{
	double absa = a > 0 ? a : -a;
	double absb = b > 0 ? b : -b;
	if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	else return (absb == 0 ? 0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

//LU decomposition of matrix U. Return ture or false depending on whether the number of row interhanges was 
// even or odd;
bool svdcmp(ZMatrix<double>& U, ZVector<double>& W, ZMatrix<double>& V)
{
#define SIGN(a,b) ((b) >= 0 ? (a) : -(a))
	
	int flag,i,its,j,jj,k,l=0,nm=0;
	double c,f,h,s,x,y,z;
	double anorm=0, g=0, scale=0;

	int m = U.NRows();
	int n = U.NCols();

	ZVector<double> rv1(n);

	//Householder reduction to bidiagonal form.
 
	// 1. Householder reduction to bidiagonal form
	for(i=0; i<n; i++) 
	{
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0;
		if (i < m) 
		{
			for (k=i;k<m;k++) scale += Abs(U[k][i]);
			if (scale != 0.0) 
			{
				for (k=i; k<m; k++) 
				{
					U[k][i] /= scale;
					s += U[k][i]*U[k][i];
				}
				f = U[i][i];
				g = -SIGN(sqrt(s),f);
				h = f * g - s;
	
				U[i][i] = f-g;
				if (i != n) 
				{
					for (j=l; j<n; j++) 
					{
						for (s=0,k=i;k<m;k++) s += U[k][i] * U[k][j];
						f = s/h;
						for (k=i;k<m;k++) U[k][j] += f * U[k][i];
					}
				}
				for (k=i;k<m;k++) U[k][i] *= scale;
			}
		}

		W[i] = scale * g;
		g = s = scale = 0;

		if (i < m && i != n) 
		{
			for (k=l; k<n; k++) scale += Abs(U[i][k]);
			
			if (scale!=0.0) 
			{
				for (k=l; k<n; k++) 
				{
					U[i][k] /= scale;
					s += U[i][k] * U[i][k];
				}
				
				f = U[i][l];
				g = -SIGN(sqrt(s),f);
				h = f * g - s;
				U[i][l] = f-g;
	
				for (k=l; k<n; k++) rv1[k] = U[i][k] / h;
				
				if (i != m) 
				{
					for (j=l; j<m; j++) 
					{
						for (s=0,k=l;k<n;k++) s += U[j][k]*U[i][k];
						for (k=l; k<n; k++) U[j][k] += s * rv1[k];
					}
				}
				for (k=l; k<n; k++) U[i][k] *= scale;
			}
		}
		double tmp = Abs(W[i]) + Abs(rv1[i]);
		anorm = Max(anorm, tmp);
	}

	// 2. Accumulation of right-hand transform V.
	for (i=n-1; i>=0; i--) 
	{
		if (i < n) 
		{
			if (g!=0.0) 
			{
				for(j=l; j<n; j++)
					V[j][i] = (U[i][j] / U[i][l]) / g;		// double div to avoid underflow
				for(j=l; j<n; j++) 
				{
					for(s=0,k=l; k<n; k++) s += U[i][k] * V[k][j];
					for(k=l;k<n;k++) V[k][j] += s * V[k][i];
				}
			}
			for (j=l; j<n; j++) V[i][j] = V[j][i] = 0;
		}
		V[i][i] = 1;
		g = rv1[i];
		l = i;
	}

	// 3. Accumulation of left-hand transform U.
	for(i=n-1; i>=0; i--) 
	{
		l = i+1;
		g=W[i];
		if (i < n-1) 
			for(j=l; j<n; j++) U[i][j] = 0;
		
		if (g!=0.0) 
		{
			g = 1/g;
			if (i != n) 
			{
				for (j=l; j<n; j++) 
				{
					for (s=0,k=l; k<m; k++) s += U[k][i] * U[k][j];
					f = (s/U[i][i]) * g;
					for (k=i; k<m; k++) U[k][j] += f * U[k][i];
				}
			}
			for (j=i; j<m; j++) U[j][i] *= g;
		} 
		else 
			for (j=i;j<m;j++) U[j][i] = 0;
		++U[i][i];
	}

	// 4. Diagonalization of bidiagonal form
	for (k=n-1; k>=0; k--) 
	{ // loop over singular values
	
		for (its=1; its<=30; its++) 
		{	// loop over allowed iterations
			flag=1;
			for (l=k; l>=0; l--) 
			{			// test for splitting
				nm = l-1;
				if((Abs(rv1[l])+anorm) == anorm) 
				{
					flag=0;
					break;
				}
				
				if ((Abs(W[nm])+anorm) == anorm) break;
			}
			if (flag) 
			{
				c=0;
				s=1;
				for (i=l;i<=k;i++) 
				{
					f = s * rv1[i];
					rv1[i] = c*rv1[i];
					if ((Abs(f)+anorm) == anorm) break;
					g = W[i];
					h = PYTHAG(f,g);
					W[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for(j=0; j<m; j++) 
					{
						y = U[j][nm];
						z = U[j][i];
						U[j][nm] = y * c + z * s;
						U[j][i] = z * c - y * s;
					}
				}
			}
			
			z = W[k];
			if (l == k) 
			{				// convergence
				if (z < 0) 
				{				// singular val is made non negative
					W[k] = -z;
					for (j=0; j<n; j++) V[j][k] = -V[j][k];
				}
				break;
			}
			
			if (its == 30) 
			{
				cout << "No convergence in 30 SVDCMP iterations" << endl;
				return false;
			}
			
			x = W[l];					// shift from bottom 2x2 minor
			nm = k-1;
			y = W[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y);
			g = PYTHAG(f,1);
			f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c = s = 1;					// next QR transformation
			for (j=l; j<=nm; j++) 
			{
				i = j + 1;
				g = rv1[i];
				y = W[i];
				h = s*g;
				g = c*g;
				z = PYTHAG(f,h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				
				for(jj=0; jj<n; jj++) 
				{
					x = V[jj][j];
					z = V[jj][i];
					V[jj][j] = x*c+z*s;
					V[jj][i] = z*c-x*s;
				}
				
				z = PYTHAG(f,h);
				W[j] = z;
				if (z!=0.0) 
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = (c*y) - (s*g);
				for (jj=0; jj<m; jj++) 
				{
					y = U[jj][j];
					z = U[jj][i];
					U[jj][j] = y*c+z*s;
					U[jj][i] = z*c-y*s;
				}
			}
			rv1[l] = 0;
			rv1[k] = f;
			W[k] = x;
		}
	}
	return true;
	
#undef SIGN
}

bool _SVD(ZMatrix<double>& U, ZVector<double>& W, ZMatrix<double>& V)
{
    UINT n = U.NCols();						// dimension of vectors in V.
	V.resize(n, n);
	W.resize(n);

	if(svdcmp(U, W, V))
	{
		double max_w = 0.0;
		
		for (UINT k=0; k < n; k++) 
		{	// find max absolute weight
			double weight = Abs(W[k]);		// in diagonal
			if (max_w < weight) max_w = weight;
		}
		
		double lbound_w = max_w * VERYTINY;

		for (UINT l=0; l < n; l++) 	// find near-singular weights
			if (Abs(W[l]) < lbound_w) W[l] = 0;

		return true;
	}
	return false;
}

bool _inv(ZMatrix<double>& U)
{
	ZMatrix<double> V;
	ZVector<double> W;

	if(_SVD(U, W, V) == false)
	{
		ZError("Matrix Inverse", "Not invertable!");
		return false;
	}

	ZMatrix<double> WI(U.NCols(), U.NCols());
	
	ZMatrix<double>::row_iterator p_wi = WI.pbegin();
	ZVector<double>::const_iterator p_w = W.pbegin();
	for(UINT i=0; p_w!=W.pend(); p_w++, p_wi++, i++)
	{
		if(*p_w > VERYTINY) (*p_wi)[i] = 1.0 / *p_w;
	}

	U = V * WI * trans(U);
	return true;
}
