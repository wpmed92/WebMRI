/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file defines those general/basic constants used by the whole library
// and should be used by all the applications.

#ifndef __MYDEFINE_H__
#define __MYDEFINE_H__

#ifdef _MSC_VER
#pragma warning( disable : 4710 )  
#endif

#define TINY			1e-10
#define VERYTINY		1e-20
#define HUGHTINY		1e-30
#define LARGE			1e10
#define VERYHUGE		1e100

#ifdef _WINDOWS
	#include <afx.h>
#else
	#ifndef _WINDOWS_
		typedef unsigned char BYTE;
		typedef unsigned char* PBYTE;
		typedef unsigned short WORD;
		typedef unsigned int UINT;
		typedef float* PFLOAT;
	#endif 
#endif 

#undef	NULL
#define NULL	0

#define LOWORD(l)           ((WORD)(l))
#define HIWORD(l)           ((WORD)(((DWORD)(l) >> 16) & 0xFFFF))
#define LOBYTE(w)           ((BYTE)(w))
#define HIBYTE(w)           ((BYTE)(((WORD)(w) >> 8) & 0xFF))

#define CEIL(num,den)		(((num)+(den)-1)/(den) )
#define ODD(x)				((x) & 0x01)
#define EVEN(x)				(!ODD(x))

template <class T> inline double MaxValue() { return 255; }
template <> inline double MaxValue<unsigned char>() { return 255; }
template <> inline double MaxValue<char>() { return 127; }
template <> inline double MaxValue<short>() { return 32767.0; }
template <> inline double MaxValue<unsigned short>() { return 65535.0; }
template <> inline double MaxValue<int>() { return 2147483647.0; }
template <> inline double MaxValue<unsigned int>() { return 4294967295.0; }
template <> inline double MaxValue<float>() { return 2147483647.0; }
template <> inline double MaxValue<double>() { return 2147483647.0; }

template <class T> inline double MinValue() { return 0; }
template <> inline double MinValue<BYTE>() { return 0; }
template <> inline double MinValue<char>() { return -128.0; }
template <> inline double MinValue<short>() { return -32768.0; }
template <> inline double MinValue<WORD>() { return 0; }
template <> inline double MinValue<int>() { return -2147483648.0; }
template <> inline double MinValue<UINT>() { return 0; }
template <> inline double MinValue<float>() { return -2147483648.0; }
template <> inline double MinValue<double>() { return -2147483648.0; }

#endif //__MYDEFINE_H__
