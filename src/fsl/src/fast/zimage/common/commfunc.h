/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares common functions used by most applications.

#ifndef __COMMFUNC_H__
#define __COMMFUNC_H__

#include <functional>
#include <vector>
#include <algorithm>

#include "mydefine.h"

void Reverse2(void *ptr, int num=1);
void Reverse4(void *ptr, int num=1);
void Reverse8(void *ptr, int num=1);

template <class T> inline T Abs(const T& a) { return a > 0 ? a : -a; }

template <class T>
inline T ipow(const T& x, int y)
{
	if(y==0) return 1;
	int yy=Abs(y); T r=x;
	for(int i=1; i<yy; i++) r *= x;
	return y>0 ? r : T(1.0/r);
}

template <class T>
inline const T& Max(const T& a, const T& b)
{
	return (a<b) ? b : a;
}

template <class T>
inline const T& Min(const T& a, const T& b)
{
	return (a<b) ? a : b;
}

template <class T>
inline T Mix(const T& x, const T& y, float s)
{
	T v = y; v -= x; v *= s; v += x;
	return v;
}

template <class T>
inline T Clip(const T& x, const T& min, const T& max)
{
	if (x < min) return min;
	else if (x > max) return max;
	return x;
}

template <class T>
inline void Reverse(T& val)
{
	int size=sizeof(T);
	if(size == 2)
		std::swap(*(PBYTE(&val)), *(PBYTE(&val)+1));
	else if(size == 4)
	{
		std::swap(*(PBYTE(&val)), *(PBYTE(&val)+3));
		std::swap(*(PBYTE(&val)+1), *(PBYTE(&val)+2));
	}
	else if(size == 8)
	{
		std::swap(*(PBYTE(&val)), *(PBYTE(&val)+7));
		std::swap(*(PBYTE(&val)+1), *(PBYTE(&val)+6));
		std::swap(*(PBYTE(&val)+2), *(PBYTE(&val)+5));
		std::swap(*(PBYTE(&val)+3), *(PBYTE(&val)+4));
	}
}

template<class For> 
inline float Mean(For first, For last)
{
	float sum  = *first;
	int n;
	for (n=1; ++first < last; n++) sum += *first;
	return sum / n;
}

template<class T, class For> 
inline T median(For first, For last)
{
	std::vector<T> buf;
	for(For i=first; ++i<last;) buf.push_back(*i);
	sort(buf.begin(), buf.end());

	return buf[buf.size()/2];
}

template<class For> 
inline void statistics(For first, For last, float& mean, float& variance, bool nozero=false)
{
	mean = variance = 0; int n;
	For ptr = first;
	for (n=0; ptr < last; ptr++) if(!nozero || *ptr) mean += *ptr, n++;
	if(n==0) return;

	mean /= n;
	if(n==1) return;

	ptr = first;
	for (; ptr < last; ptr++) if(!nozero || *ptr) variance += (*ptr - mean) * (*ptr - mean);
	variance /= n-1;
}

template<class For> 
inline float correlation(For x1, For x2, For y1, float mx, float vx, float my, float vy)
{
	if(vx <= 0 || vy <=0 ) return 0;
	float cor=0; int n=0;
	for (For p_x=x1, p_y=y1; p_x!=x2; p_x++, p_y++, n++) cor += (*p_x - mx) * (*p_y - my);
	if(n<2) return 0;
	cor /= n-1;
	return float(cor/(sqrt(vx * vy)));
}

template<class For> 
inline float scorrelation(For x1, For x2, For y1, bool nozero=false)
{
	float mx=0, my=0, vx=0, vy=0, cor=0;
	For p_x=x1, p_y=y1;  int n=0;
	for (; p_x!=x2; p_x++, p_y++) if(!nozero || (*p_x && *p_y)) mx+=*p_x, my+=*p_y, n++;
	if(n==0) return 0;
	mx /= n; my /= n;
	if(--n==0) return 0;
	for (p_x = x1, p_y = y1; p_x!=x2; p_x++, p_y++) 
	if(!nozero || (*p_x && *p_y)) 
	{ float x=(*p_x-mx), y=(*p_y-my); vx+=x*x, vy+=y*y; cor+=x*y; }
	return float(cor/(sqrt(vx * vy)));
}

#endif
