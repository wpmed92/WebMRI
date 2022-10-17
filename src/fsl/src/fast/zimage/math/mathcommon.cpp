/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <iostream> // C++ I/O
#include <cstdlib>
#include "zmath.h"
using namespace std;

ZVector<float> RandVector(UINT r)
{
	ZVector<float> m(r);

	for(UINT i=0; i<r; i++)
	{
		m[i] = float(rand()) / (RAND_MAX + 1);
	}

	return m;
}

ZMatrix<float> RandMatrix(UINT r, UINT c)
{
	ZMatrix<float> m(r, c);

	for(UINT i=0; i<r; i++)
	for(UINT j=0; j<c; j++)
	{
		m(i, j) = float(rand()) / (RAND_MAX + 1);
	}

	return m;
}

ostream& operator << (ostream& s, const ZMatrix<int>& m)
{
	s << endl; 
	for (UINT i=0; i<m.NRows(); i++)
	{
		for (UINT j=0; j<m.NCols(); j++)
			s << setw(10) << m[i][j] << ' ';
		s << endl;
	}
	return s ;
}

ostream& operator << (ostream& s, const ZMatrix<float>& m)
{
	s << endl; 
	for (UINT i=0; i<m.NRows(); i++)
	{
		for (UINT j=0; j<m.NCols(); j++)
		{
			int pre = int(Abs(log10(m[i][j])));
			if(pre > 5) pre = 0;
			else pre = 5 - pre;
			s << setw(10)  << setiosflags(ios::fixed) << setprecision(pre) << m[i][j] << ' ';
		}
		s << endl;
	}
	return s ;
}

ostream& operator << (ostream& s, const ZMatrix<double>& m)
{
	s << endl; 
	for (UINT i=0; i<m.NRows(); i++)
	{
		for (UINT j=0; j<m.NCols(); j++)
		{
			int pre = int(Abs(log10(m[i][j])));
			if(pre > 4) pre = 0;
			else pre = 4 - pre;
			s << setw(10)  << setiosflags(ios::fixed) << setprecision(pre) << m[i][j] ;
		}
		s << endl;
	}
	return s ;
}

ostream& operator << (std::ostream& s, const ZMatrix2D& m)
{
	s << endl; 
	for (int i=0; i<2; i++)
	{
		for (int j=0; j<2; j++)
			s << setw(10) << m(i, j) << ' ';
		s << endl;
	}
	return s ;
}

ostream& operator << (std::ostream& s, const ZMatrix3D& m)
{
	s << endl; 
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
			s << setw(10) << m(i, j) << ' ';
		s << endl;
	}
	return s ;
}

ostream& operator << (std::ostream& s, const ZMatrix4D& m)
{
	s << endl; 
	for (int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)
			s << setw(10) << m(i, j) << ' ';
		s << endl;
	}
	return s ;
}

void laplacian(double* image,double* lap, int nx, int ny)
{
	int index, i, j;

	for(i=1;i<ny-1;i++)
	for(j=1;j<nx-1;j++)
	{
		index = (i*nx) +j;
		lap[index] = (image[(i-1)*nx + (j-1)] + 
						image[(i-1)*nx + (j+1)] + 
						image[(i+1)*nx + (j+1)] + 
						image[(i+1)*nx + (j-1)] +
						(4.0 * (image[index+1] + image[index+nx] + 
						image[index-1] + image[index-nx])) +
						(-20.0 * image[index])) / 6.0;
	}
	for(j=0;j<nx;j++)
	{
		lap[(ny-1)*nx + j] = 0.0;
		lap[j] = 0.0;
	}
	for(i=0;i<ny;i++)
	{
		lap[(i*nx) + nx-1] = 0.0;
		lap[i*nx] = 0.0;
	}
}

void horn_derives(double* image1, double* image2, double* dx, double* dy, double* dt, int nx,int ny)
{
	int index, i, j;

	for(i=0;i<ny-1;i++)
	for(j=0;j<nx-1;j++)
	{
		index = (i*nx) +j;
		dy[index] = ((image1[index+nx] + image2[index+nx] + image1[index+nx+1] + image2[index+nx+1]) 
					-(image1[index] + image2[index] + image1[index+1] + image2[index+1])) / 4.0;
		dx[index] = ((image1[index+1] + image2[index+1] + image1[index+nx+1] + image2[index+nx+1]) 
					-(image1[index] + image2[index] + image1[index+nx] + image2[index+nx])) / 4.0;
		dt[index] = ((image2[index] + image2[index+1] + image2[index+nx] + image2[index+nx+1])  
					-(image1[index] + image1[index+1] + image1[index+nx] + image1[index+nx+1])) / 4.0;
	}

	for(j=0;j<nx;j++)
	{
		index = (ny-1)*nx + j;
		dx[index] = 0.0;
		dy[index] = 0.0;
		dt[index] = 0.0;
	}
	for(i=0;i<ny;i++)
	{
		index = (i*nx) + (nx-1);
		dx[index] = 0.0;
		dy[index] = 0.0;
		dt[index] = 0.0;
	}
}
