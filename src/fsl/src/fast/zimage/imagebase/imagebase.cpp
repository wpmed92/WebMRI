/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// The realization of the ZImageBase class.

#include <fstream>
 
#include "imagebase.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////
//
// @mfunc:(IMPL) | ZImageBase | ZImageBase |
//
// Constructors.  The default constructor does not allocate memory to
// store image pixel values.
//
// @syntax ZImageBase(EPixFmt epixfmt, UINT cbPixel)
// @syntax ZImageBase(UINT width, UINT height, UINT depth, EPixFmt epixfmt, ColorType ctype, BaseColor color);
// @syntax ZImageBase(const ZImageBase& refimage)
//
////////////////////////////////////////////////////////////////////////////

ZImageBase::ZImageBase(EPixFmt epixfmt)
	:	m_ValidStamp(123456789), m_pImageBuf(0), m_pImageLim(0), 
		m_nWidth(0), m_nHeight(0), m_nDepth(0), m_nSliceSize(0), m_nSize(0),
		m_nBytesPerRow(0), m_nBytesPerSlice(0), m_epixfmt(epixfmt),	
		m_cbPixel(PixelSize(epixfmt))
{
	m_ImageAtt.Reset();
}

ZImageBase::ZImageBase(UINT width, UINT height, UINT depth, PBYTE img, EPixFmt epixfmt)
	:	m_ValidStamp(123456789), m_pImageBuf(0), m_pImageLim(0), 
		m_nWidth(0), m_nHeight(0), m_nDepth(0), m_nSliceSize(0), m_nSize(0),
		m_nBytesPerRow(0), m_nBytesPerSlice(0), m_epixfmt(epixfmt),	
		m_cbPixel(PixelSize(epixfmt))
{
	m_ImageAtt.Reset();
	Allocate(width, height, depth, img);
}

ZImageBase::ZImageBase(const ZImageBase& refimage)
	:	m_ValidStamp(123456789), m_pImageBuf(0), m_pImageLim(0), 
		m_nWidth(0), m_nHeight(0), m_nDepth(0), m_nSliceSize(0), m_nSize(0),
		m_nBytesPerRow(0), m_nBytesPerSlice(0),	m_epixfmt(refimage.m_epixfmt),
		m_cbPixel(refimage.m_cbPixel)
{
	if(refimage.IfValid() && refimage.GetBuffer()!=0)
	{
		Allocate(refimage.Width(), refimage.Height(), refimage.Depth());
		for(UINT k=0; k < Depth(); k++)
		{
			PBYTE	ptrDst = PbSlice(k);
			PBYTE	ptrSrc = refimage.PbSlice(k);
			for(UINT j=0; j<Height(); j++, ptrSrc+=refimage.BytesPerRow(), ptrDst+=BytesPerRow())
				memcpy(ptrDst, ptrSrc, BytesPerRow());
		}
		SetAttribute(refimage.GetAttribute());
	}
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc  | ZImageBase | ~ZImageBase |
//
// Memory allocation is only done in the base class, so the
// destructor does not need to be virtual.
//
////////////////////////////////////////////////////////////////////////////
ZImageBase::~ZImageBase(void)
{
	CleanUp();
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc void | ZImageBase | Allocate |
//
// Allocate or reallocate the memory block for this image.  Most
// users of the library will not need to use this method.
//
// @syntax bool Allocate(UINT width, UINT height, UINT depth, PBYTE img);
//
// @parm | PBYTE img |
// Reference to the image whosh information should be assigned to this
// image.
//
// The image is always to be set to zero. The Clear function can also 
// be called to clear the image to zero.
////////////////////////////////////////////////////////////////////////////
bool ZImageBase::Allocate(UINT width, UINT height, UINT depth, PBYTE img)
{
	m_nSize= width*height*depth;
	if(m_nSize == 0) 
	{
		ZError("ImageBase Allocate", "Zero size given! (%d x %d x %d)", width, height, depth);
		return false;
	}

	UINT size = m_nSize * BytesPerPixel();
	delete []m_pImageBuf;

	if(img != NULL) m_pImageBuf = img;
	else 
	{
		m_pImageBuf = new BYTE[size];
		memset(m_pImageBuf, 0, size);
	}

	if(m_pImageBuf == NULL) return false;

	m_nWidth = width;
	m_nHeight = height;
	m_nDepth = depth;
	m_nSliceSize = m_nWidth * m_nHeight;
	m_nBytesPerRow = m_nWidth * m_cbPixel;
	m_nBytesPerSlice = m_nSliceSize * m_cbPixel;
	m_pImageLim = m_pImageBuf + size;

	SetROI(0, width-1, 0, height-1, 0, depth-1);
	
	return true;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc bool | ZImageBase | SetROI |
//
// Set new coordinates for the ROI.  
//
// @syntax bool SetROI(UINT left, UINT right, UINT top, UINT bottom, UINT front, UINT back)
//
// @parm UINT | left |
// The coordinate which will refer to the left of this image after this
// method is called.
//
// @parm UINT | top |
// The coordinate which will refer to the top of this image after this
// method is called.
//
// @parm UINT | front |
// The coordinate which will refer to the front slice of this image after this
// method is called.
//
// @parm UINT | right |
// The coordinate which will refer to the right of this image after this
// method is called.
//
// @parm UINT | bottom |
// The coordinate which will refer to the bottom of this image after this
// method is called.
//
// @parm UINT | back |
// The coordinate which will refer to the back slice of this image after this
// method is called.
//
////////////////////////////////////////////////////////////////////////////
bool ZImageBase::SetROI(UINT left, UINT right, UINT top, UINT bottom, UINT front, UINT back) const
{
	if(left >= m_nWidth) left = 0;
	if(top >= m_nHeight) top = 0;
	if(front >= m_nDepth) front = 0;
	if(right >= m_nWidth) right = m_nWidth - 1;
	if(bottom >= m_nHeight) bottom = m_nHeight - 1;
	if(back >= m_nDepth) back = m_nDepth - 1;

	if(left > right) swap(left, right);
	if(top > bottom) swap(top, bottom);
	if(front > back) swap(front, back);

	m_ImageROI.left = left;
	m_ImageROI.right = right;
	m_ImageROI.top = top;
	m_ImageROI.bottom = bottom;
	m_ImageROI.front = front;
	m_ImageROI.back = back;
	m_ImageROI.width = right - left + 1;
	m_ImageROI.height = bottom - top + 1;
	m_ImageROI.depth = back - front + 1;
	m_ImageROI.begin = m_pImageBuf + front * BytesPerSlice() + top * BytesPerRow() + left * m_cbPixel;
	m_ImageROI.end = m_pImageBuf + back * BytesPerSlice() + bottom * BytesPerRow() + (right+1) * m_cbPixel;
	m_ImageROI.SliceSize = m_ImageROI.width * m_ImageROI.height;
	m_ImageROI.Size = m_ImageROI.SliceSize * m_ImageROI.depth;
	m_ImageROI.BytesPerROIRow = m_ImageROI.width * m_cbPixel;
	m_ImageROI.BytesPerROISlice = m_ImageROI.SliceSize * m_cbPixel;

	return true;
}

// coordinate of a pixel given its address in the ROI
void ZImageBase::PixelOffset(PBYTE add, int& x, int& y, int& z) const
{
	if(add < m_pImageBuf || add >= m_pImageLim) x = y = z = -1;
	else
	{
		int offset = add - m_pImageBuf;
		z = offset / BytesPerSlice() - Front();
		y = offset % BytesPerSlice() / BytesPerRow() - Top();
		x = offset % BytesPerRow() / BytesPerPixel() - Left();

		if(x < 0) x = 0;
		if(y < 0) y = 0;
		if(z < 0) z = 0;
	}
}

// coordinate of a pixel given its address in the image
void ZImageBase::IPixelOffset(PBYTE add, int& x, int& y, int& z) const
{
	if(add < m_pImageBuf || add >= m_pImageLim) x = y = z = -1;
	else
	{
		int offset = add - m_pImageBuf;
		z = offset / BytesPerSlice();
		y = offset % BytesPerSlice() / BytesPerRow();
		x = offset % BytesPerRow() / BytesPerPixel();
	}
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc void | ZImageBase | Resize |
//
// Resize the Image. The only difference from Allocate is if the size is
// same as before, it will do nothing.
//
// @syntax Resize(UINT width, UINT height, UINT depth);
//
////////////////////////////////////////////////////////////
void ZImageBase::Resize(UINT width, UINT height, UINT depth)
{
	if(m_nWidth == width && m_nHeight == height && m_nDepth == depth) return;

	Allocate(width, height, depth);
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc void | ZImageBase | Replace |
//
// Replace this image with the new image. Basically, it will copy the all
// pixels of the new image to the old one. 
//
// @parm const | ZImageBase& img |
// Reference to the image whosh information is used to replace this image. The 
// pixel type has to be same as this image.
//
// @parm | UINT left, UINT top, UINT front |
// The start point of the replacing (within ROI). Pixels before this point will 
// not be modified.
//
// @syntax Replace(const ZImageBase& img, UINT left, UINT top, UINT front);
//
////////////////////////////////////////////////////////////
void ZImageBase::Replace(const ZImageBase& img, UINT left, UINT top, UINT front)
{
	if(!img.IfValid() || img.GetBuffer() == 0) return;

	if(img.PixFmt() != PixFmt()) 
	{
		ZWarning("ZImageBase::Replace", "Pixel formats are different! Nothing to do.");
		return;
	}

	if(left >= Width() || top >= Height() || front >= Depth()) 
	{
		ZWarning("ZImageBase::Replace", "Out of range!");
		return;
	}

	UINT width = Min(img.Width(), Width() - left);
	UINT height = Min(img.Height(), Height() - top);
	UINT depth = Min(img.Depth(), Depth() - front);

	PBYTE	simg = img.Origin();
	PBYTE	dimg = PbPixel(left, top, front);

	for(UINT k = 0; k < depth; k++, simg+=img.BytesPerSlice(), dimg+=BytesPerSlice())
	{
		PBYTE	sptr = simg;
		PBYTE	dptr = dimg;

		for(UINT j = 0; j < height; j++, sptr+=img.BytesPerRow(), dptr+=BytesPerRow())
			memcpy(dptr, sptr, width*BytesPerPixel());
	}
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc:(IMPL) ZImageBase& | ZImageBase | operator= |
//
// Assignment operator.
//
// @parm const | ZImageBase& refImage |
// Reference to the image whosh information should be assigned to this
// image. The pixel type of the two images has to be the same or 
// the pixel type of this image is epixfmtNil;
//
////////////////////////////////////////////////////////////////////////////
ZImageBase& ZImageBase::operator = (const ZImageBase& refimage)
{
	if(m_epixfmt != refimage.m_epixfmt && m_epixfmt != epixfmtUnknown)
	{
		ZError("Image =", "The pixel type is unmatched!");
		return *this;
	}

	if(refimage.IfValid() && refimage.GetBuffer()!=0)
	{
		Allocate(refimage.Width(), refimage.Height(), refimage.Depth());
		for(UINT k=0; k < Depth(); k++)
		{
			PBYTE	ptrDst = PbSlice(k);
			PBYTE	ptrSrc = refimage.PbSlice(k);
			for(UINT j=0; j<Height(); j++, ptrSrc+=refimage.BytesPerRow(), ptrDst+=BytesPerRow())
				memcpy(ptrDst, ptrSrc, BytesPerRow());
		}
		SetAttribute(refimage.GetAttribute());
	}

	return *this;
}

bool ZImageBase::Clone(const ZImageBase& refimage)
{
	if(m_epixfmt != refimage.m_epixfmt) return false;

	if(refimage.IfValid() && refimage.GetBuffer()!=0)
	{
		Allocate(refimage.ImageWidth(), refimage.ImageHeight(), refimage.ImageDepth());
		memcpy(m_pImageBuf, refimage.GetBuffer(), MemorySize());
		SetROI(refimage.Left(), refimage.Right(), refimage.Top(), 
			   refimage.Bottom(), refimage.Front(), refimage.Back());

		SetAttribute(refimage.GetAttribute());

		return true;
	}
	
	return false;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc void | ZImageBase | GetXYSlice |
//
// Get one slice(a 2D image) in X-Y plane from the 3D images.
//
// @parm | ZImageBase& img |
// Reference to the image which will hold the slice image.
//
// @parm | UINT z |
// the slice number.
//
// @syntax GetXYSlice(ZImageBase& img, UINT z);
//
////////////////////////////////////////////////////////////
void ZImageBase::GetXYSlice(ZImageBase& img, UINT z) const
{
	if(img.PixFmt() != PixFmt()) return;

	img.Resize(Width(), Height());

	PBYTE	simg = PbSlice(z);
	PBYTE	dimg = img.GetBuffer();

	int		bytes = BytesPerROIRow();
	for(UINT j=0; j<Height(); j++, simg+=BytesPerRow(), dimg+=bytes)
		memcpy(dimg, simg, bytes);

	img.SetAttribute(m_ImageAtt);
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc void | ZImageBase | GetYZSlice |
//
// Get one column(a 2D image) in Y-Z plane from the 3D images.
//
// @parm | ZImageBase& img |
// Reference to the image which will hold the column image.
//
// @parm | UINT x |
// the column number.
//
// @syntax GetYZSlice(ZImageBase& img, UINT x);
//
////////////////////////////////////////////////////////////
void ZImageBase::GetYZSlice(ZImageBase& img, UINT x) const
{
	if(img.PixFmt() != PixFmt()) return;

	img.Resize(Height(), Depth());

	PBYTE	simg = PbPixel(x,0)+(Height()-1)*BytesPerRow();
	PBYTE	dimg = img.GetBuffer()+(Depth()-1)*Height()*m_cbPixel;

	int		bytes = img.BytesPerRow();
	for(UINT k=0; k<Depth(); k++, simg+=BytesPerSlice(), dimg-=bytes)
	{
		PBYTE	sptr = simg;
		PBYTE	dptr = dimg;
		for(UINT j=0; j<Height(); j++, sptr-=BytesPerRow(), dptr+=m_cbPixel)
		for(UINT i=0; i<m_cbPixel; i++) 
			dptr[i] = sptr[i];
	}

	img.SetAttribute(m_ImageAtt);
	img.SetPixelDim(m_ImageAtt.m_pixdim[1], m_ImageAtt.m_pixdim[2], m_ImageAtt.m_pixdim[0]);
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc void | ZImageBase | GetXZSlice |
//
// Get one slice(a 2D image) in X-Z plane from the 3D images.
//
// @parm | ZImageBase& img |
// Reference to the image which will hold the row image.
//
// @parm | UINT y |
// the row number.
//
// @syntax GetXZSlice(ZImageBase& img, UINT y);
//
////////////////////////////////////////////////////////////
void ZImageBase::GetXZSlice(ZImageBase& img, UINT y) const
{
	if(img.PixFmt() != PixFmt()) return;

	img.Resize(Width(), Depth());

	PBYTE	simg = PbColZeroOfRow(y);
	PBYTE	dimg = img.GetBuffer()+(Depth()-1)*Width()*m_cbPixel;

	int		bytes = BytesPerROIRow();
	for(UINT k=0; k<Depth(); k++, simg+=BytesPerSlice(), dimg-=bytes)
		memcpy(dimg, simg, bytes);

	img.SetAttribute(m_ImageAtt);
	img.SetPixelDim(m_ImageAtt.m_pixdim[0], m_ImageAtt.m_pixdim[2], m_ImageAtt.m_pixdim[1]);
}

void ZImageBase::GetYZVolume(ZImageBase& img) const
{
	if(img.PixFmt() != PixFmt()) return;

	img.Resize(Height(), Depth(), Width());

	for(UINT x=0; x<Width(); x++)
	{
		PBYTE	simg = PbPixel(x,0)+(Height()-1)*BytesPerRow();
		PBYTE	dimg = img.PbSlice(x)+(Depth()-1)*Height()*m_cbPixel;

		for(UINT k=0; k<Depth(); k++, simg+=BytesPerSlice(), dimg-=img.BytesPerRow())
		{
			PBYTE	sptr = simg;
			PBYTE	dptr = dimg;
			for(UINT j=0; j<Height(); j++, sptr-=BytesPerRow(), dptr+=m_cbPixel)
			for(UINT i=0; i<m_cbPixel; i++) 
				dptr[i] = sptr[i];
		}
	}

	img.SetAttribute(m_ImageAtt);
	img.SetPixelDim(m_ImageAtt.m_pixdim[1], m_ImageAtt.m_pixdim[2], m_ImageAtt.m_pixdim[0]);
}

void ZImageBase::GetXZVolume(ZImageBase& img) const
{
	if(img.PixFmt() != PixFmt()) return;

	img.Resize(Width(), Depth(), Height());

	for(int y=Height()-1; y>=0; y--)
	{
		PBYTE	simg = PbPixel(0,y);
		PBYTE	dimg = img.PbSlice(Height()-1-y)+(Depth()-1)*Width()*m_cbPixel;

		for(UINT k=0; k<Depth(); k++, simg+=BytesPerSlice(), dimg-=BytesPerRow())
			memcpy(dimg, simg, BytesPerRow());
	}

	img.SetAttribute(m_ImageAtt);
	img.SetPixelDim(m_ImageAtt.m_pixdim[0], m_ImageAtt.m_pixdim[2], m_ImageAtt.m_pixdim[2]);
}
////////////////////////////////////////////////////////////////////////////
//
// @mfunc bool | ZImageBase | FillBitmapInfoHeader |
//
// Fill a Windows <BITMAPINFOHEADER> structure with information
// about this image.  
//
// @this const
//
// @parm BITMAPINFO * | bmpinfo |
// Pointer to the Windows <t BITMAPINFOHEADER> structure to be filled
// with information about this image.
//
// @parmopt UINT | cbitPixel | 0 |
// If non-zero, this specifies  the width of displayed pixels, in bits.
// Most callers can just pass in 0 to have the image calculate the
// number of bits to use from the pixel type.
//
////////////////////////////////////////////////////////////////////////////
bool ZImageBase::FillBitmapInfo(BITMAPINFO** pbmpinfo) const
{
	ZImagePalette	pal;
	UINT cbitPixel = isColor() ? 24 : 8; 

	UINT nColors = isColor() ? 0 : 256; 

	BITMAPINFO*	pbmp;
	pbmp = *pbmpinfo = (BITMAPINFO*) new char[sizeof(BITMAPINFOHEADER) + sizeof(PALETTE_ENTRY) * nColors];

	pbmp->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	pbmp->bmiHeader.biWidth = Width();
	pbmp->bmiHeader.biHeight = - (int)Height();
	pbmp->bmiHeader.biPlanes = 1;
	pbmp->bmiHeader.biBitCount = (WORD) cbitPixel;
	pbmp->bmiHeader.biCompression = BI_RGB;
	pbmp->bmiHeader.biSizeImage = 0;
	pbmp->bmiHeader.biXPelsPerMeter = 3150;
	pbmp->bmiHeader.biYPelsPerMeter = 3150;
	pbmp->bmiHeader.biClrUsed = nColors;
	pbmp->bmiHeader.biClrImportant = nColors;

	if(!isColor())
	{
		pal.CreateGrayPalette();
		memcpy(pbmp->bmiColors, pal.GetPalette(), 256*sizeof(PALETTE_ENTRY));
	}

	return true;
}

void ZImageBase::FillBitmapInfoHeader(BITMAPINFOHEADER& bmpinfoheader) const
{
	bmpinfoheader.biSize = sizeof(BITMAPINFOHEADER);
	bmpinfoheader.biWidth = Width();
	bmpinfoheader.biHeight = - (int)Height();
	bmpinfoheader.biPlanes = 1;
	bmpinfoheader.biBitCount = isColor() ? 24 : 8;
	bmpinfoheader.biCompression = BI_RGB;
	bmpinfoheader.biSizeImage = 0;
	bmpinfoheader.biXPelsPerMeter = 3150;
	bmpinfoheader.biYPelsPerMeter = 3150;
	bmpinfoheader.biClrUsed = isColor() ? 0 : 256;
	bmpinfoheader.biClrImportant = isColor() ? 0 : 256;
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc void | ZImageBase | Clear |
//
// Clear the image pixels by filling all bytes in the pixel memory
// with zero values.  This method could cause problems if used with
// some used-defined types.
//
////////////////////////////////////////////////////////////////////////////
void ZImageBase::Clear(void)
{
	if(!IfValid()) return;

	UINT	cbRow = BytesPerROIRow();
	PBYTE	ori = Origin();

	for(UINT k=0; k<Depth(); k++, ori+=BytesPerSlice())
	{
		PBYTE	ptr = ori;
		for(UINT j=0; j<Height(); j++, ptr+=BytesPerRow()) memset(ptr, 0, cbRow);
	}
}

double ZImageBase::RangeMax(void) const
{
	switch (m_epixfmt & (epixfmtNumericTypeMask | epixfmtUnsignedMask))
	{
	case epixfmtChar:	return 127;
	case epixfmtShort:	return 32767;
	case epixfmtInt:	return 2147483647.0;
	case epixfmtLong:	return 2147483647.0;
	case epixfmtFloat:	return 65535.0;
	case epixfmtDouble:	return 65535.0;
	case epixfmtUChar:	return 255;
	case epixfmtUShort: return 65535.0;
	case epixfmtUInt:	return 4294967295.0;
	case epixfmtULong:	return 4294967295.0;
	}

	return 0;
}

double ZImageBase::RangeMin(void) const
{
	int fmt = m_epixfmt & (epixfmtNumericTypeMask | epixfmtUnsignedMask);
	if(fmt == epixfmtShort || fmt == epixfmtFloat || fmt == epixfmtDouble)
		return -32767;
	return 0;
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | DeleteSlice |
//
// Delete one slice; Return the current slice NO.
// 
// @syntax UINT DeleteSlice(UINT nslice)
//
////////////////////////////////////////////////////////////////////////////
UINT ZImageBase::DeleteSlice()
{
	UINT curslice = Front();
	if(m_nDepth==1) return curslice;

	m_nDepth --;
	PBYTE newbuf = new BYTE[m_nSliceSize * m_nDepth * m_cbPixel];

	if(curslice > 0)	memcpy(newbuf, m_pImageBuf, BytesPerSlice()*curslice);

	if(curslice < m_nDepth)
	{
		memcpy(newbuf+BytesPerSlice()*curslice, m_pImageBuf + BytesPerSlice()*(curslice+1), BytesPerSlice()*(m_nDepth-curslice));
	}

	Allocate(m_nWidth, m_nHeight, m_nDepth, newbuf);

	if(curslice == m_nDepth) curslice --;

	SetROI(Left(), Right(), Top(), Bottom(), curslice, curslice);
	return curslice;
}
