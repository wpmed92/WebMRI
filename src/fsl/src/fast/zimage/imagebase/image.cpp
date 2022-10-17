/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#define __SGI_STL_INTERNAL_RELOPS

#include <cmath>
#include <functional>

#include "image.h"
#include "imagehelper.h"
using namespace std;

#define sendpixels(T) \
		ZImage<T>& dst = (ZImage<T>&)imgdst; \
		T dmin = T(imgdst.RangeMin()), dmax = T(imgdst.RangeMax()); \
		const_iterator p(*this); ZImage<T>::iterator pr(dst); \
		for(; p.more(); p++, pr++) *pr = (*p > dmax) ? dmax : ((*p < dmin) ? dmin : T(*p)); \
		break

#define copypixels(T) \
		ZImage<T>::iterator dst((ZImage<T>&)imgdst); \
		for(; dst.more(); src++, dst++) \
		*dst = (*src > srcmax) ? T(dstmax) : ((*src < srcmin) ? T(dstmin) : T((*src - srcmin) * scale + dstmin)); \
		break

#define copyrgbpixels(T) \
		ZImage< ZRGB<T> >::iterator dst((ZImage< ZRGB<T> >&)imgdst); \
		for(; dst.more(); src++, dst++) \
		*dst = (*src > srcmax) ? T(dstmax) : ((*src < srcmin) ? T(dstmin) : T((*src - srcmin) * scale + dstmin)); \
		break

#define colorcopyrgbpixels(T) \
		ZImage< ZRGB<T> >::iterator dst((ZImage< ZRGB<T> >&)imgdst); \
		for(; dst.more(); src++, dst++) \
		{	ZRGB<float> ss =*src; \
			dst->r = (ss.r > srcmax) ? T(dstmax) : ((ss.r < srcmin) ? T(dstmin) : T((ss.r - srcmin) * scale + dstmin)); \
			dst->g = (ss.g > srcmax) ? T(dstmax) : ((ss.g < srcmin) ? T(dstmin) : T((ss.g - srcmin) * scale + dstmin)); \
			dst->b = (ss.b > srcmax) ? T(dstmax) : ((ss.b < srcmin) ? T(dstmin) : T((ss.b - srcmin) * scale + dstmin)); \
		} \
		break


////////////////////////////////////////////////////////////////////////////////////

template <class TPixel>
void ZImage<TPixel>::Inflate(int left, int right, int top, int bottom, int front, int back,
							bool bFill, float val)
{
	if((left || top || right || bottom || front || back) == 0) return;

	TPixel fillvalue = TPixel(val);

	int nwidth = m_nWidth + left + right;
	int nheight = m_nHeight + top + bottom;
	int ndepth = m_nDepth + front + back;
	if (nwidth <= 0 || nheight <= 0 && ndepth <= 0) return;

	UINT	nslicesize = nwidth*nheight;
	UINT	size = nslicesize*ndepth;

	UINT	slicesize = PixelsPerSlice();

	TPixel* ImageBuf = new TPixel[size];

	if(bFill)
	{
		TPixel* ptr = ImageBuf;
		for(UINT i=0; i<size; i++, ptr++) *ptr = fillvalue;
	}

	UINT srcwidth = nwidth, srcheight = nheight, srcdepth = ndepth;
	if(left*right<0) srcwidth -= left>0 ? left : right;
	if(top*bottom<0) srcheight -= top>0 ? top : bottom;
	if(front*back<0) srcdepth -= front>0 ? front : back;

	UINT width = (srcwidth > m_nWidth) ? m_nWidth : srcwidth;
	UINT height = (srcheight > m_nHeight) ? m_nHeight : srcheight;
	UINT depth = (srcdepth > m_nDepth) ? m_nDepth : srcdepth;
	UINT linesize = width * m_cbPixel;

	UINT srcleft, dstleft, srctop, dsttop, srcfront, dstfront;
	srcleft = dstleft = srctop = dsttop = srcfront = dstfront = 0; 

	if(left > 0) dstleft = left;
	else srcleft = UINT(-left);

	if(top > 0) dsttop = top;
	else srctop = UINT(-top);

	if(front > 0) dstfront = front;
	else srcfront = UINT(-front);
	
	TPixel*	dst = ImageBuf + dstfront * nslicesize;
	TPixel*	src = IPPixel(srcleft, srctop, srcfront);

	for(UINT k=0; k<depth; k++, dst+=nslicesize, src+=slicesize)
	{
		TPixel* srcptr = src;
		TPixel* dstptr = dst + dsttop * nwidth + dstleft;

		for(UINT j=0; j<height; j++, dstptr+=nwidth, srcptr+=m_nWidth) 
			memcpy(dstptr, srcptr, linesize);

		if(UINT(nwidth) == width && UINT(nheight) == height) continue;
		
		if(bFill) continue;

		if(dstleft != 0)
		{
			dstptr = dst + dsttop*nwidth;
			for(UINT j=0; j<height; j++, dstptr+=nwidth)
			{
				srcptr = dstptr + 2 * dstleft -1;
				for(UINT i=0; i<dstleft; i++, srcptr--) 
					dstptr[i] = *srcptr;
			}
		}

		UINT dstright = right > 0 ? right : 0;
		if(dstright != 0)
		{
			dstptr = dst + dsttop*nwidth + nwidth - dstright;
			for(UINT j=0; j<height; j++, dstptr+=nwidth)
			{
				srcptr = dstptr - 1;
				for(UINT i=0; i<dstright; i++, srcptr--) 
					dstptr[i] = *srcptr;
			}
		}

		if(dsttop != 0)
		{
			dstptr = dst;
			srcptr = dstptr + (2*dsttop - 1)*nwidth;
			for(UINT j=0; j<dsttop; j++, dstptr+=nwidth, srcptr-=nwidth) 
				memcpy(dstptr, srcptr, nwidth*m_cbPixel);
		}

		UINT	dstbottom = bottom > 0 ? bottom : 0;
		if(dstbottom != 0)
		{
			dstptr = dst + (nheight - dstbottom)*nwidth;
			srcptr = dstptr - nwidth;
			for(UINT j=0; j<dsttop; j++, dstptr+=nwidth, srcptr-=nwidth) 
				memcpy(dstptr, srcptr, nwidth*m_cbPixel);
		}
	}

	if(!bFill)
	{
		dst = ImageBuf;
		src = ImageBuf + (2*front-1) * nslicesize;

		int kk;
		for(kk=0; kk<front; kk++, dst+=nslicesize, src-=nslicesize)
		{
			memcpy(dst, src, nslicesize*m_cbPixel);
		}

		dst = ImageBuf + (ndepth - 1) * nslicesize;
		src = dst - (2*back-1) * nslicesize;

		for(kk=0; kk<back; kk++, dst-=nslicesize, src+=nslicesize)
		{
			memcpy(dst, src, nslicesize*m_cbPixel);
		}
	}

	Allocate(nwidth, nheight, ndepth, (PBYTE)ImageBuf);
}

template <class TPixel>
void ZImage<TPixel>::Inflate(int left, int right, int top, int bottom, bool bFill, float val)
{
	Inflate(left, right, top, bottom, 0, 0, bFill, val);
}

template <class TPixel>
void ZImage<TPixel>::Inflate2D(int length, bool bFill, float val)
{
	Inflate(length, length, length, length, 0, 0, bFill, val);
}

template <class TPixel>
void ZImage<TPixel>::Inflate(int length, bool bFill, float val)
{
	Inflate(length, length, length, length, length, length, bFill, val);
}

template <class TPixel>
void ZImage<TPixel>::operator *= (const ZImageBase& image)
{ 
	if(image.PixFmt() == PixFmt()) operator *= ((const ZImage<TPixel>&)image);
	if(isFloat(image))
	{
		if(image.isColor()) operator *= ((const ZRGBFloatImage&)image);
		else 
		{
			if((image.PixFmt() & epixfmtNumericTypeMask) == epixfmtDouble)
				operator *= ((const ZGrayDoubleImage&)image);
			else operator *= ((const ZGrayFloatImage&)image);
		}
	}
	else { ZImage<TPixel> x; image.CopyPixelsTo(x);  operator *=(x); }
}

template <class TPixel>
void ZImage<TPixel>::operator /= (const ZImageBase& image)
{ 
	if(image.PixFmt() == PixFmt()) operator /= ((const ZImage<TPixel>&)image);
	if(isFloat(image))
	{
		if(image.isColor()) operator /= ((const ZRGBFloatImage&)image);
		else 
		{
			if((image.PixFmt() & epixfmtNumericTypeMask) == epixfmtDouble)
				operator /= ((const ZGrayDoubleImage&)image);
			else operator /= ((const ZGrayFloatImage&)image);
		}
	}
	else { ZImage<TPixel> x; image.CopyPixelsTo(x);  operator /=(x); }
}

template <class TPixel>
void ZImage<TPixel>::Diff(const ZImageBase& image)
{
	if(Size() != image.Size()) { ZError("ZImage::Diff", "Size unmatched!"); return; }

	ZImage<TPixel>* x = (ZImage<TPixel>*)&image;
	if(image.PixFmt() != PixFmt()) { x = new ZImage<TPixel>; image.SendPixelsTo(*x); }

	std::transform(begin(), end(), x->begin(), begin(), Dif<TPixel>());

	if(image.PixFmt() != PixFmt()) delete x;
}

template <class TPixel>
void ZImage<TPixel>::Brighter(const ZImageBase& image)
{
	if(Size() != image.Size()) { ZError("ZImage::Brighter", "Size unmatched!"); return; }

	ZImage<TPixel>* x = (ZImage<TPixel>*)&image;
	if(image.PixFmt() != PixFmt()) { x = new ZImage<TPixel>; image.SendPixelsTo(*x); }

	iterator p1(*this), p2(*x);
	for(; p1.more(); p1++, p2++) *p1 = Max(*p1, *p2);

	if(image.PixFmt() != PixFmt()) delete x;
}

template <class TPixel>
void ZImage<TPixel>::Darker(const ZImageBase& image, bool noback)
{
	if(Size() != image.Size()) { ZError("ZImage::Darker", "Size unmatched!"); return; }

	ZImage<TPixel>* x = (ZImage<TPixel>*)&image;
	if(image.PixFmt() != PixFmt()) { x = new ZImage<TPixel>; image.SendPixelsTo(*x); }

	iterator p1(*this), p2(*x);
	if(noback) { for(; p1.more(); p1++, p2++)	if(*p1!=TPixel(0) && *p2!=TPixel(0)) *p1 = Min(*p1, *p2); }
	else for(; p1.more(); p1++, p2++) *p1 = Min(*p1, *p2);

	if(image.PixFmt() != PixFmt()) delete x;
}

template <class TPixel>
void ZImage<TPixel>::Average(const ZImageBase& image, bool noback)
{
	if(Size() != image.Size()) { ZError("ZImage::Average", "Size unmatched!"); return; }

	ZImage<TPixel>* x = (ZImage<TPixel>*)&image;
	if(image.PixFmt() != PixFmt()) { x = new ZImage<TPixel>; image.SendPixelsTo(*x); }

	if(noback) std::transform(begin(), end(), x->begin(), begin(), AvgNoBackground<TPixel>());
	else std::transform(begin(), end(), x->begin(), begin(), Avg<TPixel>());

	if(image.PixFmt() != PixFmt()) delete x;
}

template <class TPixel>
void ZImage<TPixel>::Mask(const ZImageBase& image)
{
	if(Size() != image.Size()) { ZError("ZImage::Mask", "Size unmatched!"); return; }

	iterator p_o(*this);
	for(UINT k=0; k<Depth(); k++)
	for(UINT j=0; j<Height(); j++)
	for(UINT i=0; i<Width(); i++, p_o++)
	{
		if(!image.GetIntensity(i, j, k)) *p_o = 0;
	}
}

template <class TPixel>
void ZImage<TPixel>::Superimpose(const ZImageBase& image)
{
	if(Size() != image.Size()) { ZError("ZImage::Superimpose", "Size unmatched!"); return; }

	ZImage<TPixel>* x = (ZImage<TPixel>*)&image;
	if(image.PixFmt() != PixFmt()) { x = new ZImage<TPixel>; image.SendPixelsTo(*x); }

	iterator p_o(*this);
	const_iterator p_i(*x);
	for(; p_i.more(); p_i++, p_o++)
	{
		if(*p_i) *p_o = *p_i;
	}

	if(image.PixFmt() != PixFmt()) delete x;
}

template <class TPixel>
void ZImage<TPixel>::Mosaic(const ZImageBase& image, UINT size)
{
	if(Size() != image.Size()) { ZError("ZImage::Mosaic", "Size unmatched!"); return; }

	ZImage<TPixel>* x = (ZImage<TPixel>*)&image;
	if(image.PixFmt() != PixFmt()) { x = new ZImage<TPixel>; image.CopyPixelsTo(*x); }
	
	TPixel* ptrd = ImageOrigin();
	TPixel* ptrs = x->ImageOrigin();
	UINT  dsize = size*2;
	for(UINT k=0; k<Depth(); k++)
	{
		for(UINT j=0; j<Height(); j++)
		for(UINT i=0; i<Width(); i++, ptrs++, ptrd++)
		{
			if(((i%dsize<size)&&(j%dsize<size)) || ((i%dsize>=size)&&(j%dsize>=size)) )
				;
			else
				*ptrd = *ptrs;
		}
	}

	if(image.PixFmt() != PixFmt()) delete x;
}

template <class TPixel>
void ZImage<TPixel>::Negative()
{
	iterator ptr(*this); 
	if(is8bit(*this)) for(; ptr.more(); ptr++) 
		*ptr = TPixel(255)-*ptr;
	else
	{
		float min=0, max=256; MinMax(min, max);
		for(; ptr.more(); ptr++) *ptr = TPixel(max) - *ptr + TPixel(min);
	}
}

template <class TPixel>
ZImageBase*	ZImage<TPixel>::Reverse() const
{
	ZImageBase* pImage = FloatTypeCopyFrom(this);
	SendPixelsTo(*pImage);

	if(isColor())
	{
		ZRGBFloatImage::iterator p_f(*(ZRGBFloatImage*)pImage);
		for(; p_f.more(); p_f++) 
		{
			if(p_f->r) p_f->r = 1.0f / p_f->r;
			if(p_f->g) p_f->g = 1.0f / p_f->g;
			if(p_f->b) p_f->b = 1.0f / p_f->b;
		}
	}
	else
	{
		ZGrayFloatImage::iterator p_f(*(ZGrayFloatImage*)pImage);
		for(; p_f.more(); p_f++) if(*p_f) *p_f = 1.0f / *p_f;
	}

	return pImage;
}

template <class TPixel>
float ZImage<TPixel>::IGetIntensity(UINT x, UINT y, UINT z) const
{
	switch(m_epixfmt)
	{
		case epixfmtGrayByte:	return float(*IPbPixel(x, y, z));
		case epixfmtRGBByte:	return float(*((const ZRGBByte*)IPbPixel(x, y, z)));
		case epixfmtGrayShort:	return float(*((short*)IPbPixel(x, y, z)));
		case epixfmtRGBShort:	return float(*((const ZRGBShort*)IPbPixel(x, y, z)));
		case epixfmtGrayUShort:	return float(*((WORD*)IPbPixel(x, y, z)));
		case epixfmtRGBUShort:	return float(*((const ZRGBUShort*)IPbPixel(x, y, z)));
		case epixfmtGrayInt:	return float(*((int*)IPbPixel(x, y, z)));
		case epixfmtGrayFloat:	return float(*((float*)IPbPixel(x, y, z)));
		case epixfmtRGBFloat:	return float(*((const ZRGBFloat*)IPbPixel(x, y, z)));
		case epixfmtGrayDouble:	return float(*((double*)IPbPixel(x, y, z)));
		default: 
		ZError("ZImage::IGetIntensity", "Unsupported Pixel Format %s!", PixelTypeDescrip(m_epixfmt));
	}
	return 0;
}

template <class TPixel>
ZRGB<float> ZImage<TPixel>::IGetRGBIntensity(UINT x, UINT y, UINT z) const
{
	switch(m_epixfmt)
	{
		case epixfmtGrayByte:	return *IPbPixel(x, y, z);
		case epixfmtRGBByte:	return *((const ZRGBByte*)IPbPixel(x, y, z));
		case epixfmtGrayShort:	return *((short*)IPbPixel(x, y, z));
		case epixfmtRGBShort:	return *((const ZRGBShort*)IPbPixel(x, y, z));
		case epixfmtGrayUShort:	return *((WORD*)IPbPixel(x, y, z));
		case epixfmtRGBUShort:	return *((const ZRGBUShort*)IPbPixel(x, y, z));
		case epixfmtGrayInt:	return *((int*)IPbPixel(x, y, z));
		case epixfmtGrayFloat:	return *((float*)IPbPixel(x, y, z));
		case epixfmtRGBFloat:	return *((const ZRGBFloat*)IPbPixel(x, y, z));
		case epixfmtGrayDouble:	return *((double*)IPbPixel(x, y, z));
		default: 
		ZError("ZImage::IGetRGBIntensity", "Unsupported Pixel Format %s!", PixelTypeDescrip(m_epixfmt));
	}

	return 0;
}

template <class TPixel>
void ZImage<TPixel>::ISetIntensity(const ZRGB<float>& color, UINT x, UINT y, UINT z)
{
	if(x<m_nWidth && y<m_nHeight && z<m_nDepth)
	switch(m_epixfmt)
	{
		case epixfmtGrayByte:	*IPbPixel(x, y, z) = BYTE(color.r); return;
		case epixfmtRGBByte:	*((ZRGB<BYTE>*)IPbPixel(x, y, z)) = color; return;
		case epixfmtGrayShort:	*((short*)IPbPixel(x, y, z)) = short(color.r); return;
		case epixfmtRGBShort:	*((ZRGB<short>*)IPbPixel(x, y, z)) = color; return;
		case epixfmtGrayUShort:	*((WORD*)IPbPixel(x, y, z)) = WORD(color.r); return;
		case epixfmtRGBUShort:	*((ZRGB<WORD>*)IPbPixel(x, y, z)) = color; return;
		case epixfmtGrayInt:	*((int*)IPbPixel(x, y, z)) = int(color.r); return;
		case epixfmtGrayFloat:	*((float*)IPbPixel(x, y, z)) = color.r; return;
		case epixfmtRGBFloat:	*((ZRGB<float>*)IPbPixel(x, y, z)) = color; return;
		case epixfmtGrayDouble:	*((double*)IPbPixel(x, y, z)) = color.r; return;
		default: 
		ZError("ZImage::ISetIntensity", "Unsupported Pixel Format %s!", PixelTypeDescrip(m_epixfmt));
		return;
	}
}

///////////////////////////////////////////////////////////////////////////
// Copy pixels from this image to another image without intensity adjustment.
///////////////////////////////////////////////////////////////////////////
template <class TPixel>
void ZImage<TPixel>::SendPixelsTo(ZImageBase& imgdst) const
{
	if(imgdst.PixFmt() == epixfmtUnknown || imgdst.PixFmt() == PixFmt())
	{
		imgdst = *(const ZImageBase *)this;
		return;
	}

	imgdst.SetPixelDim(m_ImageAtt.m_pixdim[0], m_ImageAtt.m_pixdim[1], m_ImageAtt.m_pixdim[2]);

	imgdst.Resize(Width(), Height(), Depth());

	switch (imgdst.PixFmt())
	{
		default: ZError("ZImage::SendPixelsTo", "Unsupported pixel format %d!", imgdst.PixFmt()); return;
		case epixfmtGrayByte:	{ sendpixels(BYTE); }
		case epixfmtGrayShort:	{ sendpixels(short); } 
		case epixfmtGrayInt:	{ sendpixels(int); } 
		case epixfmtGrayFloat:	{ sendpixels(float); }
		case epixfmtGrayDouble:	{ sendpixels(double); } 
		case epixfmtGrayUShort:	{ sendpixels(WORD); } 

		case epixfmtRGBByte:	{ sendpixels(ZRGBByte); }
		case epixfmtRGBShort:	{ sendpixels(ZRGBShort); }
		case epixfmtRGBFloat:	{ sendpixels(ZRGBFloat); }
		case epixfmtRGBUShort:	{ sendpixels(ZRGBUShort); }
	}
}

template <class TPixel>
void ZImage<TPixel>::operator >> (ZImageBase& imgdst) const
{
	ZImage<TPixel>::CopyPixelsTo(imgdst);
}

///////////////////////////////////////////////////////////////////////////
// Copy pixels from this image to another image with intensity adjustment
// if they have different pixel type or it is required to do so.
///////////////////////////////////////////////////////////////////////////
template <class TPixel>
void ZImage<TPixel>::CopyPixelsTo(ZImageBase& imgdst, float dstmin, float dstmax) const
{
	if((imgdst.PixFmt() == epixfmtUnknown || imgdst.PixFmt() == PixFmt()) && dstmin == 0 && dstmax == 0)
	{
		imgdst = *(const ZImageBase *)this;
		return;
	}

	imgdst.SetPixelDim(m_ImageAtt.m_pixdim[0], m_ImageAtt.m_pixdim[1], m_ImageAtt.m_pixdim[2]);

	imgdst.Resize(Width(), Height(), Depth());

	float minIntSrc=0, maxIntSrc=256; MinMax(minIntSrc, maxIntSrc);

	float minRangeDst = isFloat(imgdst) ? (isFloat(*this) ? minIntSrc : RangeMin()) : imgdst.RangeMin();
	float maxRangeDst = isFloat(imgdst) ? (isFloat(*this) ? maxIntSrc : RangeMax()) : imgdst.RangeMax();

	float srcmin = m_ImageAtt.m_DispMin, srcmax = m_ImageAtt.m_DispMax;
	if((srcmin == srcmax) || (isFloat(imgdst)))
	{
		if(isColor())
		{
			float min1=0, min2=0, min3=0, max1=256, max2=256, max3=256;
			MinMax(min1, max1, false, 1);
			MinMax(min1, max1, false, 2);
			MinMax(min1, max1, false, 3);
			srcmin = minIntSrc = Min(ZRGB<float>(min1, min2, min3));
			srcmax = maxIntSrc = Max(ZRGB<float>(max1, max2, max3));
		}
		else { srcmin = minIntSrc;  srcmax = maxIntSrc; }
	}

	if(dstmin == 0 && dstmax == 0) 
	{
		dstmin = minRangeDst; dstmax = maxRangeDst;
		if(!isFloat(*this))
		{
			if(maxRangeDst >= RangeMax()) dstmax = srcmax, dstmin = srcmin;
		}
	}

	if(dstmin < minRangeDst) dstmin = minRangeDst;
	if(dstmax > maxRangeDst) dstmax = maxRangeDst;

	float scale = (srcmax==srcmin) ? 1 : (dstmax - dstmin) / (srcmax - srcmin);
	if(scale < 0) scale = -scale;

	const_iterator src(*this);
	switch (imgdst.PixFmt())
	{
		default: ZError("ZImage::CopyPixelsTo", "Unsupported pixel format %d!", imgdst.PixFmt()); return;
		case epixfmtGrayByte:	{ copypixels(BYTE); } 
		case epixfmtGrayShort:	{ copypixels(short); } 
		case epixfmtGrayInt:	{ copypixels(int); } 
		case epixfmtGrayFloat:	{ copypixels(float); }
		case epixfmtGrayDouble:	{ copypixels(double); } 
		case epixfmtGrayUShort:	{ copypixels(WORD); } 

		case epixfmtRGBByte:	{ if(isColor()) { colorcopyrgbpixels(BYTE); } else { copyrgbpixels(BYTE); }	}
		case epixfmtRGBShort:	{ if(isColor()) { colorcopyrgbpixels(short); } else { copyrgbpixels(short);} }
		case epixfmtRGBFloat:	{ if(isColor()) { colorcopyrgbpixels(float); } else { copyrgbpixels(float);} }
		case epixfmtRGBUShort:	{ if(isColor()) { colorcopyrgbpixels(WORD); } else { copyrgbpixels(WORD);	} }
	}
}

#ifdef _WINDOWS

template <class TPixel>
void ZImage<TPixel>::CopyPixelsTo(ZGDImage& gdimage, bool full) const
{
	ImageRect backROI = m_ImageROI;

	float mind  = m_ImageAtt.m_DispMin, maxd = m_ImageAtt.m_DispMax;
	if(m_ImageAtt.m_DispMin == -1 && m_ImageAtt.m_DispMax == -1)
	{
		if(is8bit(*this)) { m_ImageAtt.m_DispMin = 0; m_ImageAtt.m_DispMax = 255; }
		else
		{
			FullROI();
			MinMax(m_ImageAtt.m_DispMin, m_ImageAtt.m_DispMax);
		}
	}

	m_ImageROI = backROI;

	GotoSlice(Front());
	if(full) FullSlice();
	gdimage.DeleteContent();
	
	UINT awidth,width;

	int PixelSize = 3;
	if(isColor()) awidth = (Width() * 24 + 31) /32 * 4;
	else awidth = (Width() + 3) /4 *4, PixelSize = 1;

	typedef typename BYTEType<TPixel>::type byte_type;
	ZImage<byte_type> gimage;

	CopyPixelsTo(gimage, 0, 255);

	PBYTE buf = new BYTE[awidth * Height()];
	PBYTE	src = gimage.GetBuffer(), dst = buf;
	for(UINT j=0; j<Height(); j++, src += Width()*PixelSize, dst+=awidth)
	{
		memcpy(dst, src, Width()*PixelSize);
	}

	FillBitmapInfo(gdimage.GetBmpInfoAddress());

	gdimage.SetColor(isColor());
	if(isColor()) width = Width() * 3;
	else width = Width();

	if(isColor())
	{
		PBYTE	ptrDst = buf;
		for(UINT j=0; j<Height(); j++, ptrDst+=awidth)
		for(UINT i=0; i<width; i+=3)
		{
			BYTE pixel = ptrDst[i];
			ptrDst[i] = ptrDst[i+2];
			ptrDst[i+2] = pixel;
		}
	}
	gdimage.SetBuffer(buf);
	m_ImageAtt.m_DispMin = mind; m_ImageAtt.m_DispMax = maxd;
	m_ImageROI = backROI;
}

template <class TPixel>
void ZImage<TPixel>::operator >> (ZGDImage& gdimage) const
{
	CopyPixelsTo(gdimage);
}

#endif

template <class TPixel>
void ZImage<TPixel>::MinMax(float& min, float& max, bool noback, UINT channel) const
{
	const_iterator	pixel(*this);
	Intensity<TPixel> intensity;

	min = 1e10, max = -1e10;
	for (; pixel.more(); pixel++)
	{
		float val = intensity(*pixel, channel);
		if(noback && *pixel == TPixel(m_ImageAtt.m_background)) continue;
		if(min > val) min = val;
		if(max < val) max = val;
	}
	if(min==1e10) min = max = 0;
}

template <class TPixel>
float ZImage<TPixel>::Median(bool noback, UINT channel) const
{
	const_iterator	pixel(*this);
	vector<float> buf;
	Intensity<TPixel> intensity;

	for (; pixel.more(); pixel++)
	{
		float val = intensity(*pixel, channel);
		if(noback && val == m_ImageAtt.m_background) continue;
		
		buf.push_back(val);
	}

	if(buf.size() == 0) return 0;

	sort(buf.begin(), buf.end());
	return buf[buf.size()/2];
}


template <class TPixel>
UINT ZImage<TPixel>::Statistics(float& mean, float& stddev, bool noback, UINT channel) const
{
	UINT nbpixels = 0;
	
	const_iterator	pixel(*this);
	Intensity<TPixel> intensity;

	for (mean=stddev=0; pixel.more(); pixel++)
	{
		float val = intensity(*pixel, channel);
		if(noback && val == m_ImageAtt.m_background) continue;

		mean += val;
		nbpixels ++;
	}

	if(nbpixels == 0) return 0;
	
	mean /= nbpixels;

	for (pixel.reset(); pixel.more(); pixel++)
	{
		float val = Intensity<TPixel>()(*pixel, channel);
		if(noback && val == m_ImageAtt.m_background) continue;

		stddev += (val - mean) * (val - mean);
	}
	stddev = float(sqrt(stddev / nbpixels));

	return nbpixels;
}

template <class TPixel>
void ZImage<TPixel>::Histogram(UINT* hist, UINT size, bool noback, UINT channel, float min, float max) const
{
	bool scale = false;
	if(min == 0 && max == 0) 
	{
		if(is8bit(*this)) max = 256;
		else MinMax(min, max, noback, channel);
	}
	else scale = true;

	float nor = (max == min) ? 0 : float(size) / (max - min);

	if(nor < 1 || isFloat(*this)) scale = true;

	memset(hist, 0, size*sizeof(UINT));
	const_iterator	pixel(*this);
	for (; pixel.more(); pixel++)
	{
		float val = Intensity<TPixel>()(*pixel, channel);
		if(noback && val == m_ImageAtt.m_background) continue;
		if(val > max) continue;

		if(!scale) hist[int(val-min)] ++;
		else 
		{
			float tmp = (val - min) * nor;
			int i = (tmp < size) ? int(tmp) : size-1;
			hist[i] ++;
		}
	}
}

template <class TPixel>
void ZImage<TPixel>::GetNeighbor(UINT i, UINT j, UINT k, TPixel* neighbor, int order, TPixel* curptr) const
{
	int width = m_nWidth;		//just for converting m_nWidth to signed int.
	int slice = PixelsPerSlice();
	const TPixel*	ptr = (curptr==0) ? PPixel(i, j, k) : curptr;
	i += Left();
	j += Top();
	k += Front();

	if(m_nDepth>1) neighbor += 13;
	else neighbor += 4;

	*neighbor = *ptr;
	
	neighbor[-1] = (i)								? ptr[-1] : *ptr;
	neighbor[1]  = (i<m_nWidth-1)					? ptr[1] : *ptr;
	neighbor[-3] = (j)								? ptr[-width] : *ptr;
	neighbor[3]  = (j<m_nHeight-1)					? ptr[width] : *ptr;
	neighbor[-4] = (i&&j)							? ptr[-width-1] : neighbor[-3];
	neighbor[-2] = ((i<m_nWidth-1)&&j)				? ptr[-width+1] : neighbor[-3];
	neighbor[2]  = (i&&(j<m_nHeight-1))				? ptr[width-1] : neighbor[3];
	neighbor[4]  = ((i<m_nWidth-1)&&(j<m_nHeight-1))? ptr[width+1] : neighbor[3];

	if(m_nDepth>1)
	{
		neighbor[-9] = (k)							? ptr[-slice] : *ptr;
		neighbor[9]	 = (k<m_nDepth-1)				? ptr[slice] : *ptr;

		if(order > 1)
		{
			neighbor[-13] = (i&&j&&k)				? ptr[-1-width-slice] : *ptr;
			neighbor[-12] = (j&&k)					? ptr[-width-slice] : *ptr;
			neighbor[-11] = ((i<m_nWidth-1)&&j&&k)	? ptr[1-width-slice] : *ptr;
			neighbor[-10] = (i&&k)					? ptr[-1-slice] : *ptr;
			neighbor[-8]  = ((i<m_nWidth-1)&&k)		? ptr[1-slice] : *ptr;
			neighbor[-7]  = (i&&(j<m_nHeight-1)&&k)	? ptr[-1+width-slice] : *ptr;
			neighbor[-6]  = ((j<m_nHeight-1)&&k)	? ptr[width-slice] : *ptr;
			neighbor[-5]  = ((i<m_nWidth-1)&&(j<m_nHeight-1)&&k) ? ptr[1+width-slice] : *ptr;

			neighbor[13] = (i&&j&&(k<m_nDepth-1))				? ptr[-1-width+slice] : *ptr;
			neighbor[12] = (j&&(k<m_nDepth-1))					? ptr[-width+slice] : *ptr;
			neighbor[11] = ((i<m_nWidth-1)&&j&&(k<m_nDepth-1))	? ptr[1-width+slice] : *ptr;
			neighbor[10] = (i&&(k<m_nDepth-1))					? ptr[-1+slice] : *ptr;
			neighbor[8]  = ((i<m_nWidth-1)&&(k<m_nDepth-1))		? ptr[1+slice] : *ptr;
			neighbor[7]  = (i&&(j<m_nHeight-1)&&(k<m_nDepth-1))	? ptr[-1+width+slice] : *ptr;
			neighbor[6]  = ((j<m_nHeight-1)&&(k<m_nDepth-1))	? ptr[width+slice] : *ptr;
			neighbor[5]  = ((i<m_nWidth-1)&&(j<m_nHeight-1)&&(k<m_nDepth-1)) ? ptr[1+width+slice] : *ptr;
		}
	}
}

template <class TPixel>
void ZImage<TPixel>::GetNeighbor(vector<TPixel>& neighbor, UINT i, UINT j, UINT k, int order, const TPixel* curptr) const
{
	int width = m_nWidth;		//just for converting m_nWidth to signed int.
	int slice = PixelsPerSlice();
	const TPixel*	ptr = (curptr) ? curptr : PPixel(i, j, k);
	i += Left();
	j += Top();
	k += Front();

	neighbor.clear();

	neighbor.push_back(*ptr);
	if(i>0) neighbor.push_back(ptr[-1]);
	if(i<m_nWidth-1) neighbor.push_back(ptr[1]);
	if(j>0) neighbor.push_back(ptr[-width]);
	if(j<m_nHeight-1) neighbor.push_back(ptr[width]);

	if(order>1)
	{
		if(i && j) neighbor.push_back(ptr[-width-1]);
		if(i<m_nWidth-1 && j) neighbor.push_back(ptr[-width+1]);
		if(i && j<m_nHeight-1) neighbor.push_back(ptr[width-1]);
		if((i<m_nWidth-1)&&(j<m_nHeight-1)) neighbor.push_back(ptr[width+1]);
	}

	if(m_nDepth>1)
	{
		if(k>0) neighbor.push_back(ptr[-slice]);
		if(k<m_nDepth-1) neighbor.push_back(ptr[slice]);

		if(order > 1)
		{
			if(i&&j&&k) neighbor.push_back(ptr[-1-width-slice]);
			if(j&&k) neighbor.push_back(ptr[-width-slice]);
			if((i<m_nWidth-1)&&j&&k) neighbor.push_back(ptr[1-width-slice]);
			if(i&&k) neighbor.push_back(ptr[-1-slice]);
			if((i<m_nWidth-1)&&k) neighbor.push_back(ptr[1-slice]);
			if(i&&(j<m_nHeight-1)&&k) neighbor.push_back(ptr[-1+width-slice]);
			if((j<m_nHeight-1)&&k) neighbor.push_back(ptr[width-slice]);
			if((i<m_nWidth-1)&&(j<m_nHeight-1)&&k) neighbor.push_back(ptr[1+width-slice]);

			if(i&&j&&(k<m_nDepth-1)) neighbor.push_back(ptr[-1-width+slice]);
			if(j&&(k<m_nDepth-1)) neighbor.push_back(ptr[-width+slice]);
			if((i<m_nWidth-1)&&j&&(k<m_nDepth-1)) neighbor.push_back(ptr[1-width+slice]);
			if(i&&(k<m_nDepth-1)) neighbor.push_back(ptr[-1+slice]);
			if((i<m_nWidth-1)&&(k<m_nDepth-1)) neighbor.push_back(ptr[1+slice]);
			if(i&&(j<m_nHeight-1)&&(k<m_nDepth-1)) neighbor.push_back(ptr[-1+width+slice]);
			if((j<m_nHeight-1)&&(k<m_nDepth-1)) neighbor.push_back(ptr[width+slice]);
			if((i<m_nWidth-1)&&(j<m_nHeight-1)&&(k<m_nDepth-1)) neighbor.push_back(ptr[1+width+slice]);
		}
	}
}

template <class TPixel>
TPixel ZImage<TPixel>::Bilinear(float x, float y, float z) const
{
	if(x<0) x=0;  if(y<0) y=0; if(z<0) z=0;
	if(x>m_nWidth-1) x=m_nWidth-1; 
	if(y>m_nHeight-1) y=m_nHeight-1;
	if(z>m_nDepth-1) z=m_nDepth-1;
	
	UINT ix = int(x);
	UINT iy = int(y);
	UINT iz = int(z);
	float dx = x - ix;
	float dy = y - iy;
	float dz = z - iz;

	const TPixel*	ptr = IPPixel(ix, iy, iz);
	typename FloatType<TPixel>::type v0 = *ptr, v1, v2, v3, v4, v5, v6, v7;
	if(ix != m_nWidth-1 && iy != m_nHeight-1)
	{
		v1 = ptr[1];
		v2 = ptr[1+m_nWidth];
		v3 = ptr[m_nWidth];
	}
	else 
	{
		if(ix == m_nWidth-1 && iy != m_nHeight-1)
		{
			v1 = *ptr;
			v3 = v2 = ptr[m_nWidth];
		}
		else if(ix != m_nWidth-1 && iy == m_nHeight-1)
		{
			v2 = v1 = ptr[1];
			v3 = *ptr;
		}
		else v1 = v2 = v3 = *ptr;
	}

	if(m_nDepth==1)
	{
		return TPixel(v0 + (v1-v0)*dx + (v3-v0)*dy + (v0-v3-v1+v2)*dx*dy);
	}
	
	if(iz != m_nDepth-1)
	{
		ptr += m_nSliceSize;
		v4 = *ptr;
		if(ix != m_nWidth-1 && iy != m_nHeight-1)
		{
			v5 = ptr[1];
			v6 = ptr[1+m_nWidth];
			v7 = ptr[m_nWidth];
		}
		else 
		{
			if(ix == m_nWidth-1 && iy != m_nHeight-1)
			{
				v5 = *ptr;
				v7 = v6 = ptr[m_nWidth];
			}
			else if(ix != m_nWidth-1 && iy == m_nHeight-1)
			{
				v6 = v5 = ptr[1];
				v7 = *ptr;
			}
			else v5 = v6 = v7 = *ptr;
		}
	}
	else
	{
		v4 = v0;
		v5 = v1;
		v6 = v2;
		v7 = v3;
	}

	return TPixel(v0 + (v1-v0)*dx + (v3-v0)*dy + (v4-v0)*dz 
			+ (v0-v3-v1+v2)*dx*dy + (v0-v1+v5-v4)*dx*dz + (v0-v3+v7-v4)*dy*dz
			+ (v1+v3-v0-v2+v4-v7-v5+v6)*dx*dy*dz);
}

template <class TPixel>
void ZImage<TPixel>::CGV(int& x, int& y, int& z) const
{
	int n=0; 
	x = Width()/2; y = Height()/2; z = Depth()/2;

	TPixel*	ptr = ImageOrigin();

	for(UINT k=0; k<Depth(); k++, ptr+=PixelsPerSlice())
	{
		TPixel* ptr1 = ptr;
		for(UINT j=0; j<Height(); j++, ptr1+=ImageWidth())
		{
			TPixel* ptr2 = ptr1;
			for(UINT i=0; i<Width(); i++, ptr2++)
			{
				if(*ptr2 > 0)
				{
					x += i; y += j; z += k; n++;
				}
			}
		}
	}
	if(n>0)	{ x /= n; y /= n; z /= n; }
}

template <class TPixel>  
void ZImage<TPixel>::Mirror()
{
	TPixel *Ori = ImageOrigin();
	
	int	slicesize = PixelsPerSlice();

	for(UINT k=0; k<Depth(); k++, Ori+=slicesize)
	{
		TPixel* ptr = Ori;
		for(UINT j=0; j<Height(); j++, ptr+=ImageWidth())
		for(UINT i=0; i<Width()/2; i++)
		{
			TPixel tmp = ptr[i];
			ptr[i] = ptr[Width()-i-1];
			ptr[Width()-i-1] = tmp;
		}
	}
}

template <class TPixel>
void ZImage<TPixel>::Flip()
{
	TPixel *Ori = ImageOrigin();

	UINT size = BytesPerROIRow();

	vector<BYTE> tmp(size);

	int		slicesize = PixelsPerSlice();
	int		rowsize = BytesPerRow();
	
	for(UINT k=0; k<Depth(); k++, Ori+=slicesize)
	{
		PBYTE	ptr1 = PBYTE(Ori);
		PBYTE	ptr2 = PBYTE(Ori) + (Height() - 1) * rowsize;
		for(UINT j=0; j<Height()/2; j++, ptr1+=rowsize, ptr2-=rowsize)
		{
			copy(ptr1, ptr1+size, tmp.begin());
			copy(ptr2, ptr2+size, ptr1);
			copy(tmp.begin(), tmp.end(), ptr2);
		}
	}
}

template <class TPixel>
void ZImage<TPixel>::Rotate(float angle, bool keepsize)
{
	int width = Width(), height = Height();

	angle = 360 - angle;
	angle *= (3.14159265f / 180.0f);

	float cos_angle = float(cos(angle));
	float sin_angle = float(sin(angle));

	float center_row = float(height)/2;
	float center_col = float(width)/2;

	int new_width = width, new_height = height;
	if(!keepsize)
	{
		int		minx, miny, maxx, maxy;
		maxx = maxy = 0;
		minx = miny = 65536;

		for(int j=0; j<height; j+=height-1)
		for(int i=0; i<width; i+=width-1)
		{
			float xx = cos_angle * ((double)i-center_col) - sin_angle * ((double)j-center_row) + center_col;
			float yy = cos_angle * ((double)j-center_row) + sin_angle * ((double)i-center_col) + center_row;

			int ix, iy;
			if(xx > 0) ix = (int)(xx + 0.5);
			else ix = (int)(xx - 0.5);
			if(yy > 0) iy = (int)(yy + 0.5);
			else iy = (int)(yy - 0.5);

			if(minx>ix) minx = ix;
			if(maxx<ix) maxx = ix;
			if(miny>iy) miny = iy;
			if(maxy<iy) maxy = iy;
		}

		new_width = maxx - minx + 1;
		new_height = maxy - miny + 1;
	}

	float cr = float(new_height)/2;
	float cc = float(new_width)/2;

	int front = Front(), back = Back();
	ZImage<TPixel> image(new_width, new_height, Depth());

	ImageRect roi = m_ImageROI;
	TPixel* ptr = image.ImageOrigin();
	for(int k=front; k<=back; k++)
	{
		GotoSlice(k);
		for(int j=0; j<new_height; j++)
		{
			float relative_row = j - cr;
			for(int i=0; i<new_width; i++, ptr++)
			{
				float relative_col = i - cc;

				float rotate_row = cos_angle * relative_row + sin_angle * relative_col;
				float rotate_col = cos_angle * relative_col - sin_angle * relative_row;
				
				if((rotate_col < -center_col) || (rotate_col > width-center_col-1)) continue;
				if((rotate_row < -center_row) || (rotate_row > height-center_row-1)) continue;

				rotate_row += center_row;
				rotate_col += center_col;

				*ptr = Bilinear(rotate_col, rotate_row, k);
			}
		}
	}

	m_ImageROI = roi;
	if(keepsize) Replace(image);
	else { image.SetAttribute(GetAttribute()); 	*this = image; }
}

template <class TPixel>
void ZImage<TPixel>::Rescale(UINT width, UINT height, int depth)
{
	if(depth == -1) depth = m_nDepth;

	if(width == m_nWidth && height == m_nHeight && depth == int(m_nDepth)) return;

	float stepx = float(m_nWidth) / width;
	float stepy = float(m_nHeight) / height;
	float stepz = float(m_nDepth) / depth;

	float sx,sy,sz;
	GetPixelDim(sx, sy, sz);
	sx *= stepx;
	sy *= stepy;
	sz *= stepz;

	ZImage<TPixel> res(width, height, depth);
	res.SetPixelDim(sx, sy, sz);

	int ix, iy, iz, z;
	UINT x, y;
	float fx, fy, fz;

	for (fz=0.0, z=0; z<depth; z++, fz+=stepz)
	for (fy=0.0, y=0; y<height; y++, fy+=stepy)
	for (fx=0.0, x=0; x<width; x++, fx+=stepx)
	{
		ix = int(fx);  iy = int(fy);  iz = int(fz);
		if (IfInsideImage(ix,iy,iz)) 
			res(x,y,z) = Bilinear(fx,fy,fz);
	}

	*this = res;
}

template <class TPixel>
void ZImage<TPixel>::Rescale(float scale)
{
	UINT width = UINT(m_nWidth * scale);
	UINT height = UINT(m_nHeight * scale);
	UINT depth = UINT(m_nDepth * scale);

	if(depth == 0) depth = 1;

	Rescale(width, height, depth);
}

template <class TPixel>
void ZImage<TPixel>::IsotropicRescale(ZImageBase& output)
{
	int ix, iy, iz, x, y, z;
	float fx, fy, fz, stepx, stepy, stepz;

	float sx,sy,sz;
	GetPixelDim(sx, sy, sz);
	float scale = 1;

	stepx = scale / sx;
	stepy = scale / sy;
	stepz = scale / sz;

	int nz = ImageDepth()==1 ? 1 : int (float(ImageDepth())  / stepz);
	int ny = int (float(ImageHeight()) / stepy);
	int nx = int (float(ImageWidth()) / stepx);

	ZImage<TPixel> iso(nx, ny, nz);

	for (fz=0.0, z=0; z<nz; z++, fz+=stepz)
	for (fy=0.0, y=0; y<ny; y++, fy+=stepy)
	for (fx=0.0, x=0; x<nx; x++, fx+=stepx)
	{
		ix = (int) fx;  iy = (int) fy;  iz = (int) fz;
		if (IfInsideImage(ix,iy,iz)) 
			iso(x,y,z) = Bilinear(fx,fy,fz);
	}

	iso.SetAttribute(GetAttribute());
	iso >> output;
	output.SetAttribute(iso.GetAttribute());
	output.SetPixelDim(scale);
}

//contrast -1 -- 1; brightness : -1 -- 1
template <class TPixel>
void ZImage<TPixel>::NewDisplayRange(float contrast, float brightness) const
{
	if(contrast <= -1) contrast = 0.99f;
	if(contrast >= 1) contrast = 0.99f;
	if(brightness < -1) brightness = -1;
	if(brightness > 1) brightness = 1;

	float orange = 255, imin = 0;
	if(!is8bit(*this)) 
	{
		ImageRect backROI = m_ImageROI;
		FullROI();
		float min=0, max=256; MinMax(min, max);
		orange  =  max - (imin = min);
		m_ImageROI = backROI;
	}

	float nrange = orange * (1 - contrast);
	
	m_ImageAtt.m_DispMin = imin - orange * brightness / 2; 
	m_ImageAtt.m_DispMax = m_ImageAtt.m_DispMin + nrange;
}

//contrast -1 -- 1; brightness : -1 -- 1
template <class TPixel>
void ZImage<TPixel>::IdealDisplaySetting(float& contrast, float& brightness) const
{
	contrast = brightness = 0;
	float min = 1000000, max = -1000000;
	ImageRect roi = FullROI();
	const_iterator p_image(*this);
	if(is8bit(*this)) min = 0, max = 255;
	else
	{
		for(; p_image.more(); p_image++)
		{
			float v = *p_image;
			if(min > v) min = v;
			if(max < v) max = v;
		}

		if(max - min < 0.0001) return;
	}

	float nor = 255.0f / (max - min);

	int hist[256];
	memset(hist, 0, 256 * sizeof(int));
	for(p_image.reset(); p_image.more(); p_image++)
	{
		hist[int(nor * (*p_image - min))]++;
	}

	int i; float sum = 0;
	for(i=0; i<256; i++)
	{
		sum += hist[i];
		if(sum / m_nSize > 0.002) break;
	}

	float imin = float(i) / nor + min;

	for(sum=0, i=255; i>0; i--)
	{
		sum += hist[i];
		if(sum / m_nSize > 0.002) break;
	}

	float imax = float(i) / nor + min;

	float orange = max - min;

	brightness = (min - imin) * 2 / orange;
	contrast = 1 - (imax - imin) / orange;
	
	SetROI(roi);
}

template <class TPixel>
void ZImage<TPixel>::AdjustContrast(float contrast, float brightness)
{
	ZImage<TPixel> image = *this;
	image.m_ImageAtt.m_DispMin = RangeMin();
	image.m_ImageAtt.m_DispMax = RangeMax();
	float orange = image.m_ImageAtt.m_DispMax - image.m_ImageAtt.m_DispMin;

	float nrange = orange * (1 - contrast);
	
	image.m_ImageAtt.m_DispMin -= orange * brightness / 2; 
	image.m_ImageAtt.m_DispMax = image.m_ImageAtt.m_DispMin + nrange;

	image.CopyPixelsTo(*this, RangeMin(), RangeMax());
}

template <class TPixel>
void ZImage<TPixel>::MedianFilter(int radius)
{
	int radius2 = 2*radius;
	int center = (radius2+1) * (radius2+1) / 2;

	Intensity<TPixel> intensity;
	typedef typename Intensity<TPixel>::pixeltype pixeltype;
	vector<pixeltype> array((radius2+1) * (radius2+1));

	int end_chan = isColor() ? 4 : 2;

	ZImage<TPixel>	imgcpy(*this);

	imgcpy.Inflate2D(radius);

	imgcpy.SetROI(0, Width()-1, 0, Height()-1, 0, Depth()-1);

	typename ZImage<TPixel>::iterator in(imgcpy), out(*this);

	for(;in.more(); in++, out++)
	{
		for(int c=1; c<end_chan; c++)
		{
			TPixel* neighbour = in;
			for(int index=0, j=0; j<=radius2; j++,neighbour+=imgcpy.ImageWidth())
			for(int i=0; i<=radius2; i++, index++)
				array[index] = intensity(neighbour[i], c);

			sort(array.begin(), array.end());

			intensity(*out, c) = array[center];
		}
	}
}

template <class TPixel>
void ZImage<TPixel>::Boundary(ZImage<BYTE>& edge) const
{
	edge.Create(Width(), Height(), Depth());

	edge.SetROI(1, Width()-2, 1, Height()-2, 0, Depth()-1);
	ImageRect rect = GetROI();
	SetROI(Left()+1, Right()-1, Top()+1, Bottom()-1, Front(), Back());

	const_iterator p_image(*this);
	ZGrayByteImage::iterator p_edge(edge);
	
	for(; p_edge.more(); p_edge++, p_image++)
	{
		*p_edge = 255;
		if(*p_image == TPixel(0)) *p_edge = 0;
		else if((*(p_image-1)!=0) && (*(p_image+1)!=0) && (*(p_image-m_nWidth)!=0) && (*(p_image+m_nWidth)!=0))
			*p_edge = 0;
	}
	SetROI(rect);
	edge.FullROI();
}

/*
template <class TPixel>
void ZImage<TPixel>::GaussianDownSample(ZImage<TPixel>& sampled, int ratio)
{
	UINT i,j,k;
	float *ptr1,*ptr2;

	float *kernel = new float[(2*radius+1)*(2*radius+1)], sum;
	sum=0.0f;

	for(int i=-radius;i<=radius;i++)
	for(int j=-radius;j<=radius;j++)
	{
		*(kernel+(i+radius)*(2*radius+1)+j+radius) = float(pow(2.71, -(i*i+j*j)/(2*sigma*sigma)));
		sum+=*(kernel+(i+radius)*(2*radius+1)+j+radius);
	}

	int nx=Width()/ratio, ny=Height()/ratio, nz = (Depth()==1) ? 1 : Depth()/ratio;
	sampled.Create(nx, ny, nz);

	ZImage<TPixel> tmp1(nx, Height(), Depth());
	ZImage<TPixel> tmp2(nx, ny, Depth());

	int xy1 = Height() * nx;

	//row filter
	ptr1=ImageOrigin();
	ptr2=tmp1.ImageOrigin();
	
	for(k=0; k<Depth(); k++, ptr1+=PixelsPerSlice(), ptr2+=xy1)
	{
		TPixel *p1 = ptr1, *p2 = ptr2;
		for(j=0; j<Height(); j++, p1+=m_nWidth, p2+=nx)
		{
			for(i=1; i<Width()-1; i+=ratio)
				p2[i/2] = p1[i-1] * 0.25 + p1[i] * 0.5 + p1[i+1] * 0.25;

		if(image.xsize % 2 == 0) 
		  p2[sampledimage.xsize-1] = p1[image.xsize-1] * 0.75 + p1[image.xsize-2] * 0.25;
		  }
	}

	//column filter
	ptr1=tmp1;
	ptr2=tmp2;
	for(k=0;k<image.zsize;k++,ptr1+=xy1,ptr2+=sampledimage.xysize)
	{
	  float *p1 = ptr1, *p2 = ptr2;
	  for(i=0; i<sampledimage.xsize; i++,p1++, p2++)
	{
	  float *pp1 = p1+sampledimage.xsize, *pp2 = p2;
	  for(j=1; j<image.ysize-1; j+=2, pp1+=sampledimage.xsize*2, pp2+=sampledimage.xsize)
		*pp2 = pp1[-sampledimage.xsize] * 0.25 + 
		  pp1[0] * 0.5 + 
		  pp1[sampledimage.xsize] * 0.25;
	  if(image.ysize % 2 == 0) 
		*pp2 = *pp1 * 0.75 + pp1[-sampledimage.xsize] * 0.25;
	}
	}


	//slice filter
	ptr1=tmp2;
	ptr2=sampledimage.data;
	for(j=0; j<sampledimage.ysize; j++)
	for(i=0; i<sampledimage.xsize; i++, ptr1++, ptr2++)
	{
	  float *p1 = ptr1+sampledimage.xysize, *p2 = ptr2;
	  for(k=1;k<image.zsize-1;k+=2,p1+=sampledimage.xysize*2,p2+=sampledimage.xysize)
	{
	  *p2 = p1[-sampledimage.xysize] * 0.25 + p1[0] * 0.5 + p1[sampledimage.xsize] * 0.25;
	}
  
	  if(image.zsize % 2 == 0) 	
	*p2 = *p1 * 0.75 + p1[-sampledimage.xysize] * 0.25;
	}

	delete []tmp1;
	delete []tmp2;
}
*/

#ifdef _MSC_VER

template <> double Opt< ZRGB<char> >::min = -128; //MinValue<char>();
template <> double Opt< ZRGB<char> >::max = 127; //MaxValue<char>();
template <> double Opt< ZRGB<BYTE> >::min = 0; //MinValue<BYTE>();
template <> double Opt< ZRGB<BYTE> >::max = 255; //MaxValue<BYTE>();
template <> double Opt< ZRGB<short> >::min = -32768.0; //MinValue<short>();
template <> double Opt< ZRGB<short> >::max = 32767; //MaxValue<short>();
template <> double Opt< ZRGB<WORD> >::min = 0; //MinValue<WORD>();
template <> double Opt< ZRGB<WORD> >::max = 65536; // MaxValue<WORD>();
template <> double Opt< ZRGB<int> >::min = -2147483648.0; //MinValue<int>();
template <> double Opt< ZRGB<int> >::max = 2147483647.0; //MaxValue<int>();
template <> double Opt< ZRGB<UINT> >::min = 0; //MinValue<UINT>();
template <> double Opt< ZRGB<UINT> >::max = 4294967295.0; //MaxValue<UINT>();
template <> double Opt< ZRGB<float> >::min = -2147483648.0; //MinValue<float>();
template <> double Opt< ZRGB<float> >::max = 2147483647.0; //MaxValue<float>();
template <> double Opt< ZRGB<double> >::min = -2147483648.0; //MinValue<double>();
template <> double Opt< ZRGB<double> >::max = 2147483647.0; //MaxValue<double>();

#endif

/*
template class ZImage<BYTE>;
template class ZImage<short>;
template class ZImage<WORD>;
template class ZImage<float>;
template class ZImage<int>;
template class ZImage<double>;
template class ZImage< ZRGB<BYTE> >;
template class ZImage< ZRGB<short> >;
template class ZImage< ZRGB<WORD> >;
template class ZImage< ZRGB<float> >;
*/
