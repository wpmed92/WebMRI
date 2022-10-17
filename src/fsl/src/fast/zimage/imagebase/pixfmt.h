/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file defines the constant values used to describe standard
// pixel types used in images.  It also provides a templated-function
// to get a constant describing a given pixel type.

#ifndef __PIXELFMT_H__
#define __PIXELFMT_H__

#include "rgbpixels.h"

// This macro is used to identify unreferenced parameters to functions (so
// that we don't get level 4 warnings).
#ifndef Unreferenced
#define Unreferenced(x)		(x)
#endif // !Unreferenced


////////////////////////////////////////////////////////////////////////////
//  
// @enum EPixFmt |
// Standard pixel types
//
////////////////////////////////////////////////////////////////////////////
enum EPixFmt
{
	// @emem
	// Pixel type is not known.
	epixfmtUnknown			= 0x0000,

	//------------------------------------------------------------------
	// @group Numeric Types

	// @emem
	// Values for this type are <t char> or <t unsigned char> values.
	epixfmtChar				= 0x0001,

	// @emem
	// Values for this type are <t short> or <t unsigned short> values.
	epixfmtShort				= 0x0002,

	// @emem
	// Values for this type are <t int> or <t unsigned int> values.
	epixfmtInt				= 0x0003,

	// @emem
	// Values for this type are <t long> or <t unsigned long> values.
	epixfmtLong				= 0x0004,

	// @emem
	// Values for this type are <t float> values.
	epixfmtFloat			= 0x0005,

	// @emem
	// Values for this type are <t double> values.
	epixfmtDouble			= 0x0008,

	// @emem
	// Mask for the numeric type information.
	epixfmtNumericTypeMask	= 0x003f,
	//------------------------------------------------------------------
	// @group Flag for Signed or Unsigned Integer Types

	// @emem
	// Indicates that this type is <t signed>.  Not used with
	// <e EPixFmt::epixfmtFloat> or  <e EPixFmt::epixfmtDouble>.
	epixfmtSigned			= 0x0000,

	// @emem
	// Indicates that this type is <t unsigned>.  Not used with
	// <e EPixFmt::epixfmtFloat> or  <e EPixFmt::epixfmtDouble>.
	epixfmtUnsigned			= 0x0040,

	// @emem
	// Mask for the <t signed> of <t unsigned> integer type information.
	epixfmtUnsignedMask		= 0x0040,
	//------------------------------------------------------------------
	// @group Unsigned Integer Types

	// @emem
	// Values for this type are <t unsigned char> values.
	epixfmtUChar			= epixfmtChar | epixfmtUnsigned,
	epixfmtByte				= epixfmtChar | epixfmtUnsigned,

	// @emem
	// Values for this type are <t unsigned short> values.
	epixfmtUShort			= epixfmtShort | epixfmtUnsigned,

	// @emem
	// Values for this type are <t unsigned int> values.
	epixfmtUInt				= epixfmtInt | epixfmtUnsigned,

	// @emem
	// Values for this type are <t unsigned long> values.
	epixfmtULong				= epixfmtLong | epixfmtUnsigned,


	//------------------------------------------------------------------
	// @group Structure of Information in Pixel

	// @emem
	// The pixel consists of a single value of this type.
	epixfmtGray				= 0x0000,


	// @emem
	// The pixel consists of a structure containing red, green, and
	// blue values of this type.
	epixfmtRGB				= 0x0100,


	// @emem
	// Mask for the pixel structure information.
	epixfmtStructureMask	= 0x0100,
	//------------------------------------------------------------------
	// @group Common Grayscale Pixel Types

	// @emem
	// Grayscale <c char> pixels
	epixfmtGrayChar = epixfmtGray | epixfmtSigned | epixfmtChar,

	// @emem
	// Grayscale <c short> pixels
	epixfmtGrayShort = epixfmtGray | epixfmtSigned | epixfmtShort,

	// @emem
	// Grayscale <c int> pixels
	epixfmtGrayInt = epixfmtGray | epixfmtSigned | epixfmtInt,

	// @emem
	// Grayscale <c long> pixels
	epixfmtGrayLong = epixfmtGray | epixfmtSigned | epixfmtLong,

	// @emem
	// Grayscale <c BYTE> pixels
	epixfmtGrayByte = epixfmtGray | epixfmtUnsigned | epixfmtChar,

	// @emem
	// Grayscale <c unsigned short> pixels
	epixfmtGrayUShort = epixfmtGray | epixfmtUnsigned | epixfmtShort,

	// @emem
	// Grayscale <c unsigned int> pixels
	epixfmtGrayUInt = epixfmtGray | epixfmtUnsigned | epixfmtInt,

	// @emem
	// Grayscale <c unsigned long> pixels
	epixfmtGrayULong = epixfmtGray | epixfmtUnsigned | epixfmtLong,

	// @emem
	// Grayscale <c float> pixels
	epixfmtGrayFloat = epixfmtGray | epixfmtFloat,

	// @emem
	// Grayscale <c double> pixels
	epixfmtGrayDouble = epixfmtGray | epixfmtDouble,

	//------------------------------------------------------------------
	// @group Common Color (RGB) Pixel Types

	// @emem
	// Color (RGB) <c char> pixels
	epixfmtRGBChar = epixfmtRGB | epixfmtSigned | epixfmtChar,

	// @emem
	// Color (RGB) <c short> pixels
	epixfmtRGBShort = epixfmtRGB | epixfmtSigned | epixfmtShort,

	// @emem
	// Color (RGB) <c int> pixels
	epixfmtRGBInt = epixfmtRGB | epixfmtSigned | epixfmtInt,

	// @emem
	// Color (RGB) <c long> pixels
	epixfmtRGBLong = epixfmtRGB | epixfmtSigned | epixfmtLong,

	// @emem
	// Color (RGB) <c BYTE> pixels
	epixfmtRGBByte = epixfmtRGB | epixfmtUnsigned | epixfmtChar,

	// @emem
	// Color (RGB) <c unsigned short> pixels
	epixfmtRGBUShort = epixfmtRGB | epixfmtUnsigned | epixfmtShort,

	// @emem
	// Color (RGB) <c unsigned int> pixels
	epixfmtRGBUInt = epixfmtRGB | epixfmtUnsigned | epixfmtInt,

	// @emem
	// Color (RGB) <c unsigned long> pixels
	epixfmtRGBULong = epixfmtRGB | epixfmtUnsigned | epixfmtLong,

	// @emem
	// Color (RGB) <c float> pixels
	epixfmtRGBFloat = epixfmtRGB | epixfmtFloat,

	// @emem
	// Color (RGB) <c double> pixels
	epixfmtRGBDouble = epixfmtRGB | epixfmtDouble

};
 
inline UINT PixelSize(EPixFmt pixel)
{
	int num = (pixel & epixfmtStructureMask) == epixfmtRGB ? 3 : 1;
	int size = pixel & epixfmtNumericTypeMask;
	if(size < 6 && size > 2) size = 4;

	return num * size;
}

// Default:  Unknown pixel type.
template<class TPixel>
inline EPixFmt PixFmtGetTPixel(const TPixel)
{
	return epixfmtUnknown;
}

// Grayscale values
inline EPixFmt PixFmtGetTPixel(const signed char)
{
	return epixfmtGrayChar;
}

inline EPixFmt PixFmtGetTPixel(const signed short)
{
	return epixfmtGrayShort;
}

inline EPixFmt PixFmtGetTPixel(const signed int)
{
	return epixfmtGrayInt;
}

inline EPixFmt PixFmtGetTPixel(const signed long)
{
	return epixfmtGrayLong;
}

inline EPixFmt PixFmtGetTPixel(const unsigned char)
{
	return epixfmtGrayByte;
}

inline EPixFmt PixFmtGetTPixel(const unsigned short)
{
	return epixfmtGrayUShort;
}

inline EPixFmt PixFmtGetTPixel(const unsigned int)
{
	return epixfmtGrayUInt;
}

inline EPixFmt PixFmtGetTPixel(const unsigned long)
{
	return epixfmtGrayULong;
}

inline EPixFmt PixFmtGetTPixel(const float)
{
	return epixfmtGrayFloat;
}

inline EPixFmt PixFmtGetTPixel(const double)
{
	return epixfmtGrayDouble;
}


// RGB color values
inline EPixFmt PixFmtGetTPixel(const ZRGB<signed char>)
{
	return epixfmtRGBChar;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<signed short>)
{
	return epixfmtRGBShort;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<signed int>)
{
	return epixfmtRGBInt;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<signed long>)
{
	return epixfmtRGBLong;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<unsigned char>)
{
	return epixfmtRGBByte;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<unsigned short>)
{
	return epixfmtRGBUShort;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<unsigned int>)
{
	return epixfmtRGBUInt;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<unsigned long>)
{
	return epixfmtRGBULong;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<float>)
{
	return epixfmtRGBFloat;
}

inline EPixFmt PixFmtGetTPixel(const ZRGB<double>)
{
	return epixfmtRGBDouble;
}

template <class TPixel>
struct Intensity
{
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (TPixel& pixel, UINT channel=0){ return pixel; }
	const pixeltype& operator() (const TPixel& pixel, UINT channel=0) const { return pixel; }
};

#ifndef _MSC_VER

template <class TPixel>
struct Intensity< ZRGB<TPixel> >
{
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (ZRGB<TPixel>& pixel, UINT channel=0) { return pixel(channel); }
	pixeltype operator() (const ZRGB<TPixel>& pixel, UINT channel=0) const { return pixel(channel); }
};

#else

template <>
struct Intensity< ZRGB<BYTE> >
{
	typedef BYTE	TPixel;
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (ZRGB<TPixel>& pixel, UINT channel=0) { return pixel(channel); }
	pixeltype operator() (const ZRGB<TPixel>& pixel, UINT channel=0) const { return pixel(channel); }
};

template <>
struct Intensity< ZRGB<short> >
{
	typedef short	TPixel;
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (ZRGB<TPixel>& pixel, UINT channel=0) { return pixel(channel); }
	pixeltype operator() (const ZRGB<TPixel>& pixel, UINT channel=0) const { return pixel(channel); }
};

template <>
struct Intensity< ZRGB<WORD> >
{
	typedef WORD	TPixel;
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (ZRGB<TPixel>& pixel, UINT channel=0) { return pixel(channel); }
	pixeltype operator() (const ZRGB<TPixel>& pixel, UINT channel=0) const { return pixel(channel); }
};

template <>
struct Intensity< ZRGB<int> >
{
	typedef int	TPixel;
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (ZRGB<TPixel>& pixel, UINT channel=0) { return pixel(channel); }
	pixeltype operator() (const ZRGB<TPixel>& pixel, UINT channel=0) const { return pixel(channel); }
};

template <>
struct Intensity< ZRGB<UINT> >
{
	typedef UINT	TPixel;
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (ZRGB<TPixel>& pixel, UINT channel=0) { return pixel(channel); }
	pixeltype operator() (const ZRGB<TPixel>& pixel, UINT channel=0) const { return pixel(channel); }
};

template <>
struct Intensity< ZRGB<float> >
{
	typedef float	TPixel;
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (ZRGB<TPixel>& pixel, UINT channel=0) { return pixel(channel); }
	pixeltype operator() (const ZRGB<TPixel>& pixel, UINT channel=0) const { return pixel(channel); }
};

template <>
struct Intensity< ZRGB<double> >
{
	typedef double	TPixel;
	typedef TPixel	pixeltype;
	
	pixeltype& operator() (ZRGB<TPixel>& pixel, UINT channel=0) { return pixel(channel); }
	pixeltype operator() (const ZRGB<TPixel>& pixel, UINT channel=0) const { return pixel(channel); }
};

#endif

#endif // __PIXELFMT_H__
