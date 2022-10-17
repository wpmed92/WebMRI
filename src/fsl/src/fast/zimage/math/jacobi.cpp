/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "zmath.h"
using namespace std;

bool _Jacobi(ZMatrix<float>& a, ZVector<float>& d, ZMatrix<float>& v)
{
#define ROTATE(a, i, j, k, l, n) g=a[i][j]; h=a[k][l]; a[i][j] = g-s*(h+g*tau); \
								 a[k][l] = h+s*(g-h*tau);

	UINT j,iq,ip,i, n=a.NRows();
	float tresh,theta,tau,t,sm,s,h,g,c;

	ZVector<float> b(n), z(n);
	d.resize(n);
	v.resize(n, n);

	for(ip=0; ip<n; ++ip)  
	{
		// Initialize the identity matrix
		v[ip][ip] = 1;

		// Initailize b and d to be diagonal of a
		b[ip] = d[ip] = a[ip][ip];
	}
  
	for(i=0;i<50;++i)
	{
		// Sum off-diagonal elements
		sm=0;

		for(ip=0; ip<n-1; ip++)
		for(iq=ip+1; iq<n; iq++)
			sm += Abs(a[ip][iq]);
  
		//  convergence, normal return.
		if(sm == 0) return true;
  
		// tresh is different on first three iterations...
		tresh=(i<3) ? float(0.2*sm/(n*n)) : 0 ;
  
		for(ip=0; ip<n-1; ip++)
		{
			for(iq=ip+1; iq<n; iq++)
			{
				g = float(100.0 * Abs(a[ip][iq]));

				// After four sweeps, skip the rotation if the off-diagonal element is small
				// This test is taken directly from the text and looks a little suspect to me...

				if(i > 3 && g < TINY) a[ip][iq] = 0.0 ;
				else if(Abs(a[ip][iq]) > tresh) 
				{
					h = d[iq]-d[ip];
					if(g < TINY)	t = (Abs(a[ip][iq]) > TINY) ? (a[ip][iq])/h : 0; 
					else
					{ 
						theta = (Abs(h) < TINY) ? 0 : 0.5f * h / (a[ip][iq]);
						t = float(1.0 / (Abs(theta) + sqrt(1.0 + theta * theta)));
						if(theta < 0.0) t = -t ; 
					} 
					c = float(1.0 / sqrt(1.0+t*t));
					s = t * c;
					tau = float(s / (1.0+c));
			  
					h = t * a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0;
			  
					for(j=0;j<ip;j++)
					{
						ROTATE(a,j,ip,j,iq,n);
					}
					for(j=ip+1; j<iq; j++)
					{
						ROTATE(a,ip,j,j,iq,n);
					}
					for(j=iq+1;j<n;j++) 
					{
						ROTATE(a,ip,j,iq,j,n);
					}
					for(j=0;j<n;j++)
					{
						ROTATE(v,j,ip,j,iq,n);
					}
				}
			}
		}
		for(ip=0; ip<n; ip++)
		{
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}

	// Failed to converge!!
	return false;
#undef ROTATE
}
