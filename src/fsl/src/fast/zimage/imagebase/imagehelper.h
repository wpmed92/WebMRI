/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file provides helper functions for the image class.

#ifndef __IMAGEHELPER_H__
#define __IMAGEHELPER_H__

#include <fstream>
#include <vector>
#include <algorithm>
#include "image.h"

const char* PixelTypeDescrip(EPixFmt pixfmt);
inline const char* PixelType(const ZImageBase& image) { return PixelTypeDescrip(image.PixFmt()); }

inline bool is8bit(const ZImageBase& image) { return (image.PixFmt() & epixfmtNumericTypeMask) == epixfmtChar; }
inline bool is16bit(const ZImageBase& image) { return (image.PixFmt() & epixfmtNumericTypeMask) == epixfmtShort; }
inline bool is32bit(const ZImageBase& image) { return (image.PixFmt() & epixfmtNumericTypeMask) == epixfmtInt; }
inline bool isFloat(const ZImageBase& image) { return (image.PixFmt() & epixfmtNumericTypeMask) >= epixfmtFloat; }
inline bool isDouble(const ZImageBase& image) { return (image.PixFmt() & epixfmtNumericTypeMask) == epixfmtDouble; }
inline bool isUnsigned(const ZImageBase& image) { return (image.PixFmt() & epixfmtUnsignedMask)  == epixfmtUnsigned; }

ZImageBase*	TypeCopyFrom(const ZImageBase* image);
ZImageBase*	FloatTypeCopyFrom(const ZImageBase* image);
ZImageBase*	ShortTypeCopyFrom(const ZImageBase* image);
ZImageBase*	ByteTypeCopyFrom(const ZImageBase* image);
ZImageBase* GrayTypeCopyFrom(const ZImageBase* image);
ZImageBase* RGBTypeCopyFrom(const ZImageBase* image);
ZImageBase* RGBToGray(const ZImageBase* image);
ZImageBase* GrayToRGB(const ZImageBase* image, std::string lut="");
bool LUTexsited(std::string lut);

ZImageBase* ToByteImage(const ZImageBase* image);
ZImageBase* ToShortImage(const ZImageBase* image);
ZImageBase* ToFloatImage(const ZImageBase* image);

void SplitRGBChannel(const ZImageBase& image, ZImageBase& r, ZImageBase& g, ZImageBase& b);
void CombineRGBChannel(ZImageBase& image, const ZImageBase& r, const ZImageBase& g, const ZImageBase& b);

void exp(ZGrayFloatImage& x);
void log(ZGrayFloatImage& x);
void log10(ZGrayFloatImage& x);
void sqrt(ZGrayFloatImage& x);
void exp(ZRGBFloatImage& x);
void log(ZRGBFloatImage& x);
void log10(ZRGBFloatImage& x);
void sqrt(ZRGBFloatImage& x);

ZImageBase* exp(ZImageBase& image); 
ZImageBase* log(ZImageBase& image); 
ZImageBase* log10(ZImageBase& image);
ZImageBase* sqrt(ZImageBase& image);

#endif
