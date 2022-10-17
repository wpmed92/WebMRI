/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares the typed image classes.  These classes
// support two-dimensional, optionally multi-slice,
// arrays of pixels of any C++ or user-defined type.
//
// The size of an image is set when an image is constructed.  Subimages
// which share the same memory can be constructed from images.

#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <algorithm>

#include "imagebase.h"
#include "imageiter.h"

template <class T> class ZImage;

typedef ZImage<signed char>		ZGrayCharImage;
typedef ZImage<short>			ZGrayShortImage;
typedef ZImage<int>				ZGrayIntImage;
typedef ZImage<long>			ZGrayLongImage;
typedef ZImage<unsigned char>	ZGrayByteImage;
typedef ZImage<unsigned char>	ZGrayUCharImage;
typedef ZImage<unsigned short>	ZGrayUShortImage;
typedef ZImage<unsigned int>	ZGrayUIntImage;
typedef ZImage<unsigned long>	ZGrayULongImage;
typedef ZImage<float>			ZGrayFloatImage;
typedef ZImage<double>			ZGrayDoubleImage;


typedef ZImage< ZRGB<char> >			ZRGBCharImage;
typedef ZImage< ZRGB<short> >			ZRGBShortImage;
typedef ZImage< ZRGB<int> >				ZRGBIntImage;
typedef ZImage< ZRGB<long> >			ZRGBLongImage;
typedef ZImage< ZRGB<unsigned char> >	ZRGBByteImage;
typedef ZImage< ZRGB<unsigned char> >	ZRGBUCharImage;
typedef ZImage< ZRGB<unsigned short> >	ZRGBUShortImage;
typedef ZImage< ZRGB<unsigned int> >	ZRGBUIntImage;
typedef ZImage< ZRGB<unsigned long> >	ZRGBULongImage;
typedef ZImage< ZRGB<float> >			ZRGBFloatImage;
typedef ZImage< ZRGB<double> >			ZRGBDoubleImage;

////////////////////////////////////////////////////////////////////////////
//  
// @class
//  
// This templated class defines an image whose pixels values are of
// the type specified.  Users of the library can use this class to
// define images of any desired pixel type.
//
// Standard image types are defined in the <l ImageTypes\.h> file.
//  
// @tcarg class | TPixel | Pixel type.
//  
// @base public | ZImageBase
//
// @xref <c ZImageBase>
//
////////////////////////////////////////////////////////////////////////////

template <class TPixel>
class ZImage : public ZImageBase
{
public:
	typedef ZIterator<ZImage<TPixel>&, TPixel, TPixel*, TPixel&> iterator;
	typedef ZIterator<const ZImage<TPixel>&, TPixel, const TPixel*, const TPixel&>	const_iterator;

	iterator begin() { return iterator(*this); }
	const_iterator begin() const { return const_iterator(*this); }
	iterator end() { return iterator(*this, (TPixel*) m_ImageROI.end); }
	const_iterator end() const { return const_iterator(*this, (const TPixel*) m_ImageROI.end); }

    //------------------------------------------------------------------
	// @group Constructors and Assignement Operator

	// Default constructor.  The image will not be useable until another
	// image is assigned to it.
	ZImage(void) : ZImageBase(PixFmtGetTPixel(TPixel())) {}

	// Construct an image of a given size, allocating memory needed
	// for the image or using memory provided.
	ZImage(UINT width, UINT height, UINT depth = 1, PBYTE pbData = 0)
		:	ZImageBase(width, height, depth, pbData, PixFmtGetTPixel(TPixel())) {}
	
	// Construct an image from another image.
	ZImage(const ZImage<TPixel>& imageSrc) : ZImageBase(imageSrc) {}

	bool Create(UINT width, UINT height, UINT depth = 1, PBYTE pbData = 0)
	{ return Allocate(width, height, depth, pbData); }

	ZImage& operator = (const ZImage& refImage)
	{ 
		ZImageBase::operator = ((const ZImageBase&) refImage);
		return *this; 
	}

	//------------------------------------------------------------------
	// @group Pixel Access

	// Find the address of the pixel with the specified coordinates in ROI.
	const TPixel *PPixel(UINT x, UINT y, UINT z=0) const { return (const TPixel*)PbPixel(x, y, z); }
	TPixel *PPixel(UINT x, UINT y, UINT z=0) { return (TPixel*) PbPixel(x, y, z); }

	// Find the address of the pixel with the specified coordinates in the image.
	TPixel *IPPixel(UINT x, UINT y, UINT z=0) { return (TPixel*) IPbPixel(x, y, z); }
	const TPixel *IPPixel(UINT x, UINT y, UINT z=0) const { return (TPixel*) IPbPixel(x, y, z); }

	// Given the coordinates of a point in the image ROI, find the value
	// of the pixel with the specified coordinates.
	const TPixel& Pixel(UINT x, UINT y, UINT z=0) const { return *PPixel(x, y, z); }
	TPixel& Pixel(UINT x, UINT y, UINT z=0) { return *PPixel(x, y, z); }

	// Given the coordinates of a point in the image, find the value
	// of the pixel with the specified coordinates.
	const TPixel& IPixel(UINT x, UINT y, UINT z=0) const { return *IPPixel(x, y, z); }
	TPixel& IPixel(UINT x, UINT y, UINT z=0) { return *IPPixel(x, y, z); }

	TPixel *ImageOrigin() const { return (TPixel*) m_ImageROI.begin; }
	TPixel *ImageEnd() const { return (TPixel*) m_ImageROI.end; }
	TPixel *BufferOrigin() const { return (TPixel*) m_pImageBuf; }

	bool	GetCoordinate(const TPixel* pixel, UINT& x, UINT& y, UINT& z) const 
	{
		if((pixel < (TPixel*) m_pImageBuf) || (pixel > (TPixel*) m_pImageLim)) return false;

		int pos = pixel - (TPixel*) m_pImageBuf;

		x = UINT(pos % m_nWidth);
		y = UINT((pos % PixelsPerSlice()) / m_nWidth);
		z = UINT(pos / PixelsPerSlice());

		return true;
	}
	//------------------------------------------------------------------
	// @group Miscellaneous

//	void	GaussianDownSample(ZImage<TPixel>& sampled, int ratio=2);

	// Set each pixel in the image rectangle to the specified value.
	void	FillPixels(TPixel value)
	{ for (iterator pixel(*this); pixel.more(); pixel++) *pixel = value; }

	
	void	Inflate(int left, int right, int top, int bottom, int front, int back, bool bFill, float val);
	void	Inflate(int left, int right, int top, int bottom, bool bFill, float val);
	void	Inflate2D(int length, bool bFill=true, float val=0);
	void	Inflate(int length, bool bFill=true, float val=0);

	TPixel	Bilinear(float x, float y, float z=0) const;
	float	GrayBilinear(float x, float y, float z=0) const { return float(Bilinear(x, y, z)); }

	void	SendPixelsTo(ZImageBase&) const;
	void	operator >> (ZImageBase&) const;
	void	CopyPixelsTo(ZImageBase&, float dstmin=0, float dstmax=0) const;

#ifdef _WINDOWS
	void	operator >> (ZGDImage&) const;
	void	CopyPixelsTo(ZGDImage&, bool full=true) const;
#endif

	void	GetNeighbor(UINT i, UINT j, UINT k, TPixel* neighbor, int order, TPixel* curptr=0) const;
	void	GetNeighbor(std::vector<TPixel>& neighbor, UINT i, UINT j, UINT k, int order, const TPixel* curptr=0) const;

	void	MinMax(float& min, float& max, bool noback=false, UINT channel = 0) const;
	float	Minimum(bool noback=false, UINT channel = 0) const { float min=0, max=0; MinMax(min, max); return min; }
	float	Maximum(bool noback=false, UINT channel = 0) const { float min=0, max=0; MinMax(min, max); return max; }
	float	Median(bool noback=false, UINT channel = 0) const;
	UINT	Statistics(float& mean, float& stddev, bool noback=false, UINT channel=0) const;
	void	Histogram(UINT* hist, UINT size, bool noback=false, UINT channel=0, float min=0, float max=0) const;
	void	CGV(int& x, int& y, int& z) const;		//Center of Gravity

	TPixel&	operator () (UINT x, UINT y, UINT z=0) { return *PPixel(x,y,z); }
	const TPixel& operator () (UINT x, UINT y, UINT z=0) const { return *PPixel(x,y,z); }
	float	IGetIntensity(UINT x, UINT y, UINT z=0) const;
	ZRGB<float>	IGetRGBIntensity(UINT x, UINT y, UINT z=0) const;
	void	ISetIntensity(const ZRGB<float>& color, UINT x, UINT y, UINT z=0);
	float	GetIntensity(UINT x, UINT y, UINT z=0) const
	{ return IGetIntensity(x+Left(), y+Top(), z+Front()); }
	ZRGB<float>	GetRGBIntensity(UINT x, UINT y, UINT z=0) const
	{ return IGetRGBIntensity(x+Left(), y+Top(), z+Front()); }
	void	SetIntensity(const ZRGB<float>& color, UINT x, UINT y, UINT z=0)
	{ ISetIntensity(color, x+Left(), y+Top(), z+Front()); }

	void Negative();
	ZImageBase*	Reverse() const;

	void operator += (float val)
	{ for(iterator p(*this); p.more(); p++) *p = Opt<TPixel>::Add(*p, val); }

	void operator -= (float val)
	{ for(iterator p(*this); p.more(); p++) *p = Opt<TPixel>::Sub(*p, val); }
	
	void operator *= (float val)
	{ for(iterator p(*this); p.more(); p++) *p = Opt<TPixel>::Mul(*p, val); }

	void operator /= (float val)
	{ for(iterator p(*this); p.more(); p++) *p = Opt<TPixel>::Div(*p, val); }

	template <class TPixel2>
	void operator += (const ZImage<TPixel2>& image)
	{
		if(Size() != image.Size()) { ZError("ZImage::+=", "Size unmatched (%d vx. %d)!", Size(), image.Size()); return; }
		std::transform(begin(), end(), image.begin(), begin(), ::Add<TPixel, TPixel2>());
	}

	template <class TPixel2>
	void operator -= (const ZImage<TPixel2>& image)
	{
		if(Size() != image.Size()) { ZError("ZImage::-=", "Size unmatched (%d vx. %d)!", Size(), image.Size()); return; }
		std::transform(begin(), end(), image.begin(), begin(), ::Sub<TPixel, TPixel2>());
	}

	template <class TPixel2>
	void operator *= (const ZImage<TPixel2>& image)
	{
		if(Size() != image.Size()) { ZError("ZImage::*=", "Size unmatched (%d vx. %d)!", Size(), image.Size()); return; }
		std::transform(begin(), end(), image.begin(), begin(), ::Mul<TPixel, TPixel2>());
	}

	template <class TPixel2>
	void operator /= (const ZImage<TPixel2>& image)
	{
		if(Size() != image.Size()) { ZError("ZImage::/=", "Size unmatched (%d vx. %d)!", Size(), image.Size()); return; }
		std::transform(begin(), end(), image.begin(), begin(), ::Div<TPixel, TPixel2>());
	}

	void operator += (const ZImageBase& image) 
	{ 
		if(image.PixFmt() == PixFmt()) operator += ((ZImage<TPixel>&)image);
		else { ZImage<TPixel> y; image.CopyPixelsTo(y); operator += (y); }
	}
	void operator -= (const ZImageBase& image)
	{ 
		if(image.PixFmt() == PixFmt()) operator -= ((ZImage<TPixel>&)image);
		else { ZImage<TPixel> y; image.CopyPixelsTo(y); operator -= (y); }
	}
	void operator *= (const ZImageBase& image);
	void operator /= (const ZImageBase& image);
	void operator &= (const ZGrayByteImage& mask) 
	{
		iterator p_image(*this); ZImage<BYTE>::const_iterator p_mask(mask);
		for(; p_image.more(); p_image++, p_mask++) if(*p_mask == 0) *p_image = 0;
	}

	void Diff(const ZImageBase& image);
	void Brighter(const ZImageBase& image);
	void Darker(const ZImageBase& image, bool noback=false);
	void Average(const ZImageBase& image, bool noback=false);
	void Mask(const ZImageBase& image);
	void Superimpose(const ZImageBase& image);
	void Mosaic(const ZImageBase& image, UINT size=20);

	void Mirror();
	void Flip();
	void Rotate(float angle, bool keepsize=true);
	void Rescale(UINT width, UINT height, int depth=-1);
	void Rescale(float scale);
	void IsotropicRescale(ZImageBase& iso);

	void GaussianNoise(float mean, float stddev)
	{	
		set_noraml_rand_seed(-10);
		for(iterator pixel(*this); pixel.more(); pixel++) 
		*pixel = Opt<TPixel>::Add(*pixel, normal_random(mean, stddev));
	}

	void RandomNoise(float min, float max) 
	{ 	for(iterator pixel(*this); pixel.more(); pixel++) 
		*pixel = Opt<TPixel>::Add(*pixel, rrand(min, max)); 
	}

	void RayleighNoise(float stddev)
	{ for(iterator pixel(*this); pixel.more(); pixel++) 
		*pixel = Opt<TPixel>::Add(*pixel, rayleigh_random(stddev));
	}

	void MedianFilter(int radius);

	void NewDisplayRange(float contrast, float brightness) const;
	void IdealDisplaySetting(float& contrast, float& brightness) const;
	void AdjustContrast(float contrast, float brightness);

	void Boundary(ZImage<BYTE>& edge) const;
};

#endif // __IMAGE_H__
