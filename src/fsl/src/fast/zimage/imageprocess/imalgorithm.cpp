/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */


#include <cmath>
#include <vector>

#include "mydefine.h"
#include "common.h"
#include "imalgorithm.h"

using namespace std;

void AutoRange(const ZImageBase& image, int& left, int& right, int& top, int& bottom, int& front, int& back)
{
	bottom = 0; top = int(image.Height()-1);
	right = 0;	left = int(image.Width()-1);
	back = 0;	front = int(image.Depth()-1);

	ZGrayByteImage bimage;
	image.SendPixelsTo(bimage);
	int i, j, k;
	BYTE* ptr = bimage.GetBuffer();
	for(k=0; k<int(image.Depth()); k++)
	for(j=0; j<int(image.Height()); j++)
	for(i=0; i<int(image.Width()); i++, ptr++)
	if(*ptr)
	{
		if(i < left) left = i;
		if(i > right) right = i;
		if(j < top) top = j;
		if(j > bottom) bottom = j;
		if(k < front) front = k;
		if(k > back) back = k;
	}
	front += image.Front(); back += image.Front();
	top += image.Top(); bottom += image.Top();
	left += image.Left(); right += image.Left();
}

bool AutoCrop(ZImageBase& image)
{
	image.FullROI();

	int left, right, top, bottom, front, back;
	AutoRange(image, left, right, top, bottom, front, back);

	if(left == right || top == bottom) return false;
	
	if(((left | top | front) == 0) && (right == int(image.Width()-1))
								   && (bottom == int(image.Height()-1))
								   && (back == int(image.Depth()-1)))
		return false;

	image.Inflate(-left, right-image.Width()+1, -top, bottom-image.Height()+1, -front, back-image.Depth()+1, true, 0);
	
	return true;
}

bool AutoRange(const vector<ZImageBase*>& images, int& left, int& right, int& top, int& bottom, int& front, int& back)
{
	if(images.size() == 0) 
	{
		ZError("AutoRange", "No images were given!");
		return false;
	}

	UINT i;
	for(i=0; i<images.size(); i++) images[i]->FullROI();

	UINT width=images[0]->Width(), height=images[0]->Height(), depth=images[0]->Depth();

	for(i=1; i<images.size(); i++)
	{
		if(width != images[i]->Width() || height!=images[i]->Height() || depth!=images[i]->Depth())
		{
			ZError("AutoRange", "Image %d has a different size!", i);
			return false;
		}
	}

	left = top = front = 100000;
	right = bottom = back = 0;
	for(i=0; i<images.size(); i++)
	{
		int l, r, t, b, f, bk;
		AutoRange(*images[i], l, r, t, b, f, bk);
		if(l < left) left = l;
		if(t < top) top = t;
		if(f < front) front = f;
		if(r > right) right = r;
		if(b > bottom) bottom = b;
		if(bk > back) back = bk;
	}

	return true;
}

bool AutoCrop(vector<ZImageBase*>& images)
{
	if(images.size() == 1) return AutoCrop(*images[0]);

	for(UINT c=0; c<images.size(); c++) images[c]->FullROI();

	int left, right, top, bottom, front, back;
	if(!AutoRange(images, left, right, top, bottom, front, back)) return false;

	if(left == right || top == bottom) return false;
	
	if(((left | top | front) == 0) && 
		(right == int(images[0]->Width()-1)) && (bottom == int(images[0]->Height()-1)) && (back == int(images[0]->Depth()-1))) return false;

	int width = images[0]->Width(), height = images[0]->Height(), depth = images[0]->Depth();
	for(UINT i=0; i<images.size(); i++)
	{
		images[i]->Inflate(-left, right+1-width, -top, bottom+1-height, -front, back+1-depth, true, 0);
	}

	return true;
}

void RemoveOutlier(ZImageBase& image, float low, float high, bool keep)
{
	switch (image.PixFmt())
	{
		default: 
			ZWarning("RemoveOutlier", "Unsupported pixel format %d!", image.PixFmt()); 
			break;
		case epixfmtGrayByte:	RemoveOutlierT((ZGrayByteImage&)image, low, high, keep); break;
		case epixfmtGrayShort:	RemoveOutlierT((ZGrayShortImage&)image, low, high, keep); break;
		case epixfmtGrayInt:	RemoveOutlierT((ZGrayIntImage&)image, low, high, keep); break;
		case epixfmtGrayFloat:	RemoveOutlierT((ZGrayFloatImage&)image, low, high, keep); break;

		case epixfmtRGBByte:	RemoveOutlierT((ZRGBByteImage&)image, low, high, keep); break;
		case epixfmtRGBShort:	RemoveOutlierT((ZRGBShortImage&)image, low, high, keep); break;
		case epixfmtRGBFloat:	RemoveOutlierT((ZRGBFloatImage&)image, low, high, keep); break;
	}
}

void RemovePercent(ZImageBase& image, float low, float high, bool keep)
{
	float min=0, max=256; image.MinMax(min, max);
	float mean=0, stddev=0;
	UINT nbpixels = image.Statistics(mean, stddev, true);

	UINT hist[256];
	image.Histogram(hist, 256, true);

	int l = int(float(nbpixels)*low), h = int(float(nbpixels)*(1.0f-high));

	float nor = (max - min) / 255.0f;
	int i=0, inxlow=0, inxhigh=0;
	for(i=0; i<128; i++) if((inxlow += hist[i]) > l) break;
	low = nor * i + min;
	for(i=255; i>80; i--) if((inxhigh += hist[i]) > h) break;
	high = nor * i + min;
	RemoveOutlier(image, low, high, keep);
}

void FillPattern(ZGrayByteImage& image, UINT size)
{
	PBYTE begin = image.Origin();
	UINT  dsize = size*2;
	BYTE  fg = 150, bg = 240;

	float x, y, z;
	image.GetPixelDim(x, y, z);
	int zsize = (z == x) ? size : int(float(size) * x / z);

	for(UINT k=0; k<image.Depth(); k++, begin+=image.BytesPerSlice())
	{
		if(k % zsize == 0) swap(fg, bg);
 		PBYTE ptr = begin;
		for(UINT j=0; j<image.Height(); j++)
		for(UINT i=0; i<image.Width(); i++, ptr++)
		{
			if(((i%dsize<size)&&(j%dsize<size)) || ((i%dsize>=size)&&(j%dsize>=size)) )
				*ptr = fg;
			else
				*ptr = bg;
		}
	}
}

