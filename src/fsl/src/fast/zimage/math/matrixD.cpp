/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "matrix2D.h"
#include "matrix3D.h"
#include "matrix4D.h"
using namespace std;

ZMatrix2D&	ZMatrix2D::normalize(int dir)
{
	if(dir==1)	// row
	{
		float n1 = float(sqrt(d11*d11+d12*d12)), n2 = float(sqrt(d21*d21+d22*d22));
		d11 /= n1; d12 /= n1;  d21 /= n2; d22 /= n2; 
	}
	else if(dir==2) //column
	{
		float n1 = float(sqrt(d11*d11+d21*d21)), n2 = float(sqrt(d12*d12+d22*d22));
		d11 /= n1; d12 /= n2; d21 /= n1; d22 /= n2; 
	}
	return *this;
}

ZMatrix3D operator* (const ZMatrix3D& l, const ZMatrix3D& r)
{
	ZMatrix3D		temp;

	temp.d11 = l.d11 * r.d11 + l.d12 * r.d21 + l.d13 * r.d31;
	temp.d12 = l.d11 * r.d12 + l.d12 * r.d22 + l.d13 * r.d32;
	temp.d13 = l.d11 * r.d13 + l.d12 * r.d23 + l.d13 * r.d33;
	temp.d21 = l.d21 * r.d11 + l.d22 * r.d21 + l.d23 * r.d31;
	temp.d22 = l.d21 * r.d12 + l.d22 * r.d22 + l.d23 * r.d32;
	temp.d23 = l.d21 * r.d13 + l.d22 * r.d23 + l.d23 * r.d33;
	temp.d31 = l.d31 * r.d11 + l.d32 * r.d21 + l.d33 * r.d31;
	temp.d32 = l.d31 * r.d12 + l.d32 * r.d22 + l.d33 * r.d32;
	temp.d33 = l.d31 * r.d13 + l.d32 * r.d23 + l.d33 * r.d33;

	return temp;
}

ZMatrix3D&	ZMatrix3D::normalize(int dir)
{
	if(dir==1)	// row
	{
		float n1 = float(sqrt(d11*d11+d12*d12+d13*d13));
		float n2 = float(sqrt(d21*d21+d22*d22+d23*d23));
		float n3 = float(sqrt(d31*d31+d32*d32+d33*d33));
		d11 /= n1; d12 /= n1; d13 /= n1; 
		d21 /= n2; d22 /= n2; d23 /= n2;
		d31 /= n3; d32 /= n3; d33 /= n3;
	}
	else if(dir==2) //column
	{
		float n1 = float(sqrt(d11*d11+d21*d21+d31*d31));
		float n2 = float(sqrt(d12*d12+d22*d22+d32*d32));
		float n3 = float(sqrt(d13*d13+d23*d23+d33*d33));
		d11 /= n1; d12 /= n2; d13 /= n3; 
		d21 /= n1; d22 /= n2; d23 /= n3;
		d31 /= n1; d32 /= n2; d33 /= n3;
	}
	return *this;
}

ZMatrix3D cov(const ZMatrix3D& v) 
{
	float m1 = (v.d11+v.d21+v.d31)/3, m2 = (v.d12+v.d22+v.d32)/2, m3 = (v.d13+v.d23+v.d33)/3;
	ZMatrix3D mm(v.d11-m1, v.d12-m2, v.d13-m3, v.d21-m2, v.d22-m2, v.d23-m3, v.d31-m2, v.d32-m2, v.d33-m3);
	ZMatrix3D cov = trans(mm) * mm;
	return cov /= 3;
}

ZMatrix3D inv(const ZMatrix3D& m)
{
	float v = det(m); if(v < TINY) v = 0;
	if (!v) { ZWarning("Matrix3D inv", "Not invertable!"); 	return ZMatrix3D(); }

	ZMatrix3D temp;

	temp.d11 = (m.d22*m.d33 - m.d23*m.d32) / v;
	temp.d12 = - (m.d12*m.d33 - m.d13*m.d32) / v;
	temp.d13 = (m.d12*m.d23 - m.d13*m.d22) / v;

	temp.d21 = - (m.d21*m.d33 - m.d23*m.d31) / v;
	temp.d22 = (m.d11*m.d33 - m.d13*m.d31) / v;
	temp.d23 = - (m.d11*m.d23 - m.d13*m.d21) / v;

	temp.d31 = (m.d21*m.d32 - m.d22*m.d31) / v;
	temp.d32 = - (m.d11*m.d32 - m.d12*m.d31) / v;
	temp.d33 = (m.d11*m.d22 - m.d12*m.d21) / v;

	return temp; 
}

#undef _VALOP
#define _VALOP(OP) \
	return ZMatrix4D(OP(v.d11), OP(v.d12), OP(v.d13), OP(v.d14), \
					 OP(v.d21), OP(v.d22), OP(v.d23), OP(v.d24), \
					 OP(v.d31), OP(v.d32), OP(v.d33), OP(v.d34), \
					 OP(v.d41), OP(v.d42), OP(v.d43), OP(v.d44))

#undef _VALOPF
#define _VALOPF(RHS) \
	d11 RHS; d12 RHS; d13 RHS; d14 RHS; \
	d21 RHS; d22 RHS; d23 RHS; d24 RHS; \
	d31 RHS; d32 RHS; d33 RHS; d34 RHS; \
	d31 RHS; d32 RHS; d33 RHS; d44 RHS; \
	return *this

#undef _VALOPM
#define _VALOPM(OP, m) \
	d11 OP m.d11; d12 OP m.d12; d13 OP m.d13; d14 OP m.d14; \
	d21 OP m.d21; d22 OP m.d22; d23 OP m.d23; d24 OP m.d24; \
	d31 OP m.d31; d32 OP m.d32; d33 OP m.d33; d34 OP m.d34; \
	d41 OP m.d31; d42 OP m.d32; d43 OP m.d33; d44 OP m.d34; \
	return *this

ZMatrix4D&	ZMatrix4D::normalize(int dir)
{
	if(dir==1)	// row
	{
		float n1 = float(sqrt(d11*d11+d12*d12+d13*d13+d14*d14));
		float n2 = float(sqrt(d21*d21+d22*d22+d23*d23+d24*d24));
		float n3 = float(sqrt(d31*d31+d32*d32+d33*d33+d34*d34));
		float n4 = float(sqrt(d41*d41+d42*d42+d43*d43+d34*d44));
		d11 /= n1; d12 /= n1; d13 /= n1; d14 /= n1;
		d21 /= n2; d22 /= n2; d23 /= n2; d24 /= n2;
		d31 /= n3; d32 /= n3; d33 /= n3; d34 /= n3;
		d41 /= n4; d42 /= n4; d43 /= n4; d44 /= n4;
	}
	else if(dir==2) //column
	{
		float n1 = float(sqrt(d11*d11+d21*d21+d31*d31+d41*d41));
		float n2 = float(sqrt(d12*d12+d22*d22+d32*d32+d42*d42));
		float n3 = float(sqrt(d13*d13+d23*d23+d33*d33+d43*d43));
		float n4 = float(sqrt(d14*d14+d24*d24+d34*d34+d44*d44));
		d11 /= n1; d12 /= n2; d13 /= n3; d14 /= n4;
		d21 /= n1; d22 /= n2; d23 /= n3; d24 /= n4;
		d31 /= n1; d32 /= n2; d33 /= n3; d34 /= n4;
		d41 /= n1; d42 /= n2; d43 /= n3; d44 /= n4;
	}
	return *this;
}

ZMatrix4D& ZMatrix4D::operator += (const ZMatrix4D& m) { _VALOPM(+=, m); }
ZMatrix4D& ZMatrix4D::operator -= (const ZMatrix4D& m) { _VALOPM(-=, m); }
ZMatrix4D& ZMatrix4D::operator *= (const ZMatrix4D& m) { _VALOPM(*=, m); }
ZMatrix4D& ZMatrix4D::operator /= (const ZMatrix4D& m) { _VALOPM(/=, m); }

ZMatrix4D& ZMatrix4D::operator += (float f) { _VALOPF(+= f); }
ZMatrix4D& ZMatrix4D::operator -= (float f) { _VALOPF(-= f); }
ZMatrix4D& ZMatrix4D::operator *= (float f) { _VALOPF(*= f); }
ZMatrix4D& ZMatrix4D::operator /= (float f) { _VALOPF(/= f); }

float det(const ZMatrix4D& m) 
{ 
	float v;  ZMatrix3D t;

	t.d11 = m.d11; t.d12 = m.d12; t.d13 = m.d13;
	t.d21 = m.d21; t.d22 = m.d22; t.d23 = m.d23;
	t.d31 = m.d31; t.d32 = m.d32; t.d33 = m.d33;
	v = m.d44 * det(t);

	t.d13 = m.d14; t.d23 = m.d24; t.d33 = m.d34;
	v -= m.d43 * det(t);

	t.d12 = m.d13; t.d22 = m.d23; t.d32 = m.d33;
	v += m.d42 * det(t);

	t.d11 = m.d12; t.d21 = m.d22; t.d31 = m.d32;
	v -= m.d41 * det(t);

	return v;
}

float sum(const ZMatrix4D& v) 
{ 
	return v.d11+v.d12+v.d13+v.d14+v.d21+v.d22+v.d23+v.d24
		  +v.d31+v.d32+v.d33+v.d34+v.d41+v.d42+v.d43+v.d44; 
}

float Min(const ZMatrix4D& v) { return *(min_element(v.begin(), v.end())); }
float Max(const ZMatrix4D& v) { return *(max_element(v.begin(), v.end())); }
float mean(const ZMatrix4D& v) { return float(sum(v))/16; }

ZVector4D mean_row(const ZMatrix4D& v) 
{
	return ZVector4D((v.d11+v.d21+v.d31+v.d41)/4, (v.d12+v.d22+v.d32+v.d42)/4, 
					 (v.d13+v.d23+v.d33+v.d34)/4, (v.d43+v.d43+v.d43+v.d44)/4); 
}

ZVector4D mean_col(const ZMatrix4D& v) 
{ 
	return ZVector4D((v.d11+v.d12+v.d13+v.d14)/4, (v.d21+v.d22+v.d23+v.d24)/4, 
					 (v.d31+v.d32+v.d33+v.d34)/4, (v.d41+v.d42+v.d43+v.d44)/4); 
}

ZMatrix4D Abs(const ZMatrix4D& v) { _VALOP(Abs); }
ZMatrix4D trans(const ZMatrix4D& v) 
{ 
	return ZMatrix4D(v.d11, v.d21, v.d31, v.d41, v.d12, v.d22, v.d32, v.d42, 
					 v.d13, v.d23, v.d33, v.d43, v.d14, v.d24, v.d34, v.d44); 
}

ZMatrix4D cov(const ZMatrix4D& v)
{
	float m1 = (v.d11+v.d21+v.d31+v.d41)/4, m2 = (v.d12+v.d22+v.d32+v.d42)/4;
	float m3 = (v.d13+v.d23+v.d33+v.d43)/4, m4 = (v.d14+v.d24+v.d34+v.d44)/4;
	ZMatrix4D mm(v.d11-m1, v.d12-m2, v.d13-m3, v.d14-m4, 
				 v.d21-m2, v.d22-m2, v.d23-m3, v.d24-m4,
				 v.d31-m2, v.d32-m2, v.d33-m3, v.d34-m4,
				 v.d41-m2, v.d42-m2, v.d43-m3, v.d44-m4);
	ZMatrix4D cov = trans(mm) * mm;
	return cov /= 4;
}

ZMatrix4D inv(const ZMatrix4D& MM)
{
	int		i, j, k, row=0, col=0;
	int		ind_x[4], ind_y[4];
	double	temp, max;
	ZMatrix4D m = MM;

	for ( k=0 ; k<4 ; k++ )
	{
		max = 0;
		for ( i=k ; i<4 ; i++ )
		{
			for ( j=k ; j<4 ; j++ )
			{
				if ( Abs(max) < Abs (m[i][j]))
				{
					row = i;
					col = j;
					max = m[i][j];
				}
			}
		}
		
		if ( Abs(max) < TINY ) return ZMatrix4D();

		if ( row != k )
		{
			for ( j=0 ; j<4 ; j++ )
			{
				temp = m[row][j];
				m[row][j] = m[k][j];
				m[k][j] = temp;
			}
		}
		if ( col != k )
		{
			for ( i=0 ; i<4 ; i++ )
			{
				temp = m[i][col];
				m[i][col] = m[i][k];
				m[i][k] = temp;
			}
		}
		ind_x[k] = row;
		ind_y[k] = col;

		m[k][k] = 1 / m[k][k];
		for ( j=0 ; j<4 ; j++ )
		{
			if ( j != k )
				m[k][j] = m[k][j] * m[k][k];
		}
		for ( i=0 ; i<4 ; i++ )
		{
			if ( i != k )
			{
				for ( j=0 ; j<4 ; j++ )
				{
					if ( j != k )
						m[i][j] = m[i][j] - m[i][k] * m[k][j];
				}
			}
		}
		for ( i=0 ; i<4 ; i++ )
		{
			if ( i != k )
				m[i][k] = - m[i][k] * m[k][k];
		}
	}

	for ( k=3 ; k>=0 ; k-- )
	{
		row = ind_x[k];
		col = ind_y[k];
		if ( row != k )
		{
			for ( i=0 ; i<4 ; i++ )
			{
				temp = m[i][row];
				m[i][row] = m[i][k];
				m[i][k] = temp;
			}
		}
		if ( col != k )
		{
			for ( j=0 ; j<4 ; j++ )
			{
				temp = m[col][j];
				m[col][j] = m[k][j];
				m[k][j] = temp;
			}
		}
	}

	return m;
}

ZMatrix4D operator + (const ZMatrix4D& v) { return v; }
ZMatrix4D operator - (const ZMatrix4D& v) { _VALOP(-); }
ZMatrix4D operator + (const ZMatrix4D& l, const ZMatrix4D& r) { ZMatrix4D v(l); return v+=r; }
ZMatrix4D operator - (const ZMatrix4D& l, const ZMatrix4D& r) { ZMatrix4D v(l); return v-=r; }
ZMatrix4D operator / (const ZMatrix4D& l, const ZMatrix4D& r) { ZMatrix4D v(l); return v/=r; }
	
ZMatrix4D operator + (const ZMatrix4D& l, float r) { ZMatrix4D v(l); return v+=r; }
ZMatrix4D operator - (const ZMatrix4D& l, float r) { ZMatrix4D v(l); return v-=r; }
ZMatrix4D operator * (const ZMatrix4D& l, float r) { ZMatrix4D v(l); return v*=r; }
ZMatrix4D operator / (const ZMatrix4D& l, float r) { ZMatrix4D v(l); return v/=r; }

ZMatrix4D cos(const ZMatrix4D& v) { _VALOP(cos); }
ZMatrix4D sin(const ZMatrix4D& v) { _VALOP(sin); }
ZMatrix4D tan(const ZMatrix4D& v) { _VALOP(tan); }
ZMatrix4D cosh(const ZMatrix4D& v) { _VALOP(cosh); }
ZMatrix4D sinh(const ZMatrix4D& v) { _VALOP(sinh); }
ZMatrix4D tanh(const ZMatrix4D& v) { _VALOP(tanh); }
ZMatrix4D acos(const ZMatrix4D& v) { _VALOP(acos); }
ZMatrix4D asin(const ZMatrix4D& v) { _VALOP(asin); }
ZMatrix4D atan(const ZMatrix4D& v) { _VALOP(atan); }
ZMatrix4D sqrt(const ZMatrix4D& v) { _VALOP(sqrt); }
ZMatrix4D exp(const ZMatrix4D& v) { _VALOP(exp); }
ZMatrix4D log(const ZMatrix4D& v) { _VALOP(log); }
ZMatrix4D log10(const ZMatrix4D& v) { _VALOP(log10); }

ZMatrix4D operator* (const ZMatrix4D& l, const ZMatrix4D& r)
{
	ZMatrix4D		temp;

	temp.d11 = l.d11 * r.d11 + l.d12 * r.d21 + l.d13 * r.d31 + l.d14 * r.d41;
	temp.d12 = l.d11 * r.d12 + l.d12 * r.d22 + l.d13 * r.d32 + l.d14 * r.d42;
	temp.d13 = l.d11 * r.d13 + l.d12 * r.d23 + l.d13 * r.d33 + l.d14 * r.d43;
	temp.d14 = l.d11 * r.d14 + l.d12 * r.d24 + l.d13 * r.d34 + l.d14 * r.d44;
	temp.d21 = l.d21 * r.d11 + l.d22 * r.d21 + l.d23 * r.d31 + l.d24 * r.d41;
	temp.d22 = l.d21 * r.d12 + l.d22 * r.d22 + l.d23 * r.d32 + l.d24 * r.d42;
	temp.d23 = l.d21 * r.d13 + l.d22 * r.d23 + l.d23 * r.d33 + l.d24 * r.d43;
	temp.d24 = l.d21 * r.d14 + l.d22 * r.d24 + l.d23 * r.d34 + l.d24 * r.d44;
	temp.d31 = l.d31 * r.d11 + l.d32 * r.d21 + l.d33 * r.d31 + l.d34 * r.d41;
	temp.d32 = l.d31 * r.d12 + l.d32 * r.d22 + l.d33 * r.d32 + l.d34 * r.d42;
	temp.d33 = l.d31 * r.d13 + l.d32 * r.d23 + l.d33 * r.d33 + l.d34 * r.d43;
	temp.d34 = l.d31 * r.d14 + l.d32 * r.d24 + l.d33 * r.d34 + l.d34 * r.d44;
	temp.d41 = l.d41 * r.d11 + l.d42 * r.d21 + l.d43 * r.d31 + l.d44 * r.d41;
	temp.d42 = l.d41 * r.d12 + l.d42 * r.d22 + l.d43 * r.d32 + l.d44 * r.d42;
	temp.d43 = l.d41 * r.d13 + l.d42 * r.d23 + l.d43 * r.d33 + l.d44 * r.d43;
	temp.d44 = l.d41 * r.d14 + l.d42 * r.d24 + l.d43 * r.d34 + l.d44 * r.d44;

	return temp;
}

ZVector4D operator* (const ZMatrix4D& l, const ZVector4D& r)
{
	return ZVector4D(l.d11*r.x+l.d12*r.y+l.d13*r.z+l.d14*r.w, 
					 l.d21*r.x+l.d22*r.y+l.d23*r.z+l.d24*r.w,
					 l.d31*r.x+l.d32*r.y+l.d33*r.z+l.d34*r.w,
					 l.d41*r.x+l.d42*r.y+l.d43*r.z+l.d44*r.w);
}

ZVector4D operator* (const ZVector4D& l, const ZMatrix4D& r)
{
	return ZVector4D(l.x*r.d11+l.y*r.d21+l.z*r.d31+l.w*r.d41, 
					 l.x*r.d12+l.y*r.d22+l.z*r.d32+l.w*r.d42,
					 l.x*r.d13+l.y*r.d23+l.z*r.d33+l.w*r.d43,
					 l.x*r.d14+l.y*r.d24+l.z*r.d34+l.w*r.d44);
}
