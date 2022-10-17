/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares the base image classes.
// Typed image classes, defined in ImageOf.h, are derived
// from these base classes
//
// The base classes support two-dimensional, optionally multi-Slice,
// arrays of pixels of any size (given in bytes).  The typed image
// classes support two-dimensional, optionally multi-Slice,
// arrays of pixels of any C++ or user-defined type.

#ifndef ___IMAGEBASE_H__
#define ___IMAGEBASE_H__

#ifdef _WINDOWS
#include "gdimage.h"
#endif

#ifdef _MSC_VER
#pragma warning( disable : 4710 4100)  
#endif

#include "mydefine.h"
#include "pixfmt.h"
#include "bitmap.h"
#include "palette.h"
#include "common.h"

struct ImageRect
{
	UINT		left, top, front, right, bottom, back;
	UINT		width, height, depth;
	PBYTE		begin,end;
	// slicesize of the Image ROI
	UINT		SliceSize;

	// size of the Image ROI
	UINT		Size;

	// Bytes per ROI row
	UINT		BytesPerROIRow;
	
	// Bytes per ROI row
	UINT		BytesPerROISlice;

	ImageRect() 
	{ 
		left = top = front = right = bottom = back = width = height = depth = 0; 
		begin = end = 0;
		SliceSize = Size = BytesPerROIRow = BytesPerROISlice = 0;
	}
};

template <class T> class ZImage;

class ZImageBase
{
	struct _ImageAttribute
	{
		float	m_pixdim[3];
		float	m_background;
		float	m_DispMin, m_DispMax;

		_ImageAttribute() { Reset(); }
		void Reset()
		{
			m_pixdim[0] = m_pixdim[1] = m_pixdim[2] = 1;
			m_background = 0;
			m_DispMin = m_DispMax = -1;
		}
	};

	unsigned long	m_ValidStamp;		//123456789 for valid object
protected:
	PBYTE		m_pImageBuf;
	PBYTE		m_pImageLim;

	//Region of Interests struct
	mutable ImageRect		m_ImageROI;


	// width of the Image Memory Block
	UINT	m_nWidth;

	// height of the Image Memory Block
	UINT	m_nHeight;

	// depth of the Image Memory Block
	UINT	m_nDepth;

	// slicesize of the Image Memory Block
	UINT	m_nSliceSize;

	// ImageSize of the Image Memory Block
	UINT	m_nSize;

	// Bytes per image row
	UINT	m_nBytesPerRow;

	// Bytes per image slice
	UINT	m_nBytesPerSlice;

	// The type of this pixel (if it is a standard pixel type).
	EPixFmt m_epixfmt;

	// The size of image pixels in bytes.
	UINT	m_cbPixel;

	mutable _ImageAttribute	m_ImageAtt;

	bool	Allocate(UINT width, UINT height, UINT depth, PBYTE img = NULL);
public:

	//------------------------------------------------------------------
	// @group Constructor, Destructor, Assign

	ZImageBase(EPixFmt epixfmt = epixfmtUnknown);
	ZImageBase(UINT width, UINT height, UINT depth, PBYTE img=NULL, EPixFmt epixfmt = epixfmtGrayByte);
	ZImageBase(const ZImageBase& refimage);

	virtual ~ZImageBase(void);

	void	CleanUp(void);
	ZImageBase& operator=(const ZImageBase& refImage);
	
	bool	Clone(const ZImageBase& image);

	//------------------------------------------------------------------
	// @group Image attribute

	EPixFmt		PixFmt(void) const;				// return the pixel format
	void		SetPixelDim(float dim) const;
	void		SetPixelDim(float xdim, float ydim, float zdim) const;
	void		GetPixelDim(float& xdim, float& ydim, float& zdim) const;

	bool		IfValid(void) const;			// if the image has buffer
	bool		isColor() const { return (m_epixfmt & epixfmtStructureMask)  == epixfmtRGB; } 

	float		GetBackground() const { return 	m_ImageAtt.m_background; }
	void		SetBackground(float back) const { m_ImageAtt.m_background = back; }

	_ImageAttribute GetAttribute() const { return m_ImageAtt; } 
	void SetAttribute(const _ImageAttribute& att) const { m_ImageAtt = att; }

	// return the max, min and complement (0) value of this pixel format in double
	double		RangeMax(void) const;
	double		RangeMin(void) const;
	void		GetDisplayRange(float& min, float& max) const { min = m_ImageAtt.m_DispMin; max = m_ImageAtt.m_DispMax; }
	void		SetDisplayRange(float min, float max) const { m_ImageAtt.m_DispMin = min; m_ImageAtt.m_DispMax = max; }

	//------------------------------------------------------------------
	// @group Image Memory modify

	void		Resize(UINT width, UINT height, UINT depth=1);
	void		Replace(const ZImageBase&, UINT left=0, UINT top=0, UINT front=0);
	void		Clear(void);

	//------------------------------------------------------------------
	// @group Image buffer size
	
	PBYTE		GetBuffer(void) const;			// Get the image buffer
	PBYTE		GetLimit(void) const;			// Get the image buffer
	UINT		ImageWidth(void) const;			// Get the Image Width
	UINT		ImageHeight(void) const;		// Get the Image Height
	UINT		ImageDepth(void) const;			// Get the Image Depth
	UINT		PixelsPerSlice(void) const;		// Get the Image Slice Size
	UINT		ImageSize(void) const;			// Get the Image Size
	UINT		MemorySize(void) const;			// Get the Memory Size

	UINT		BytesPerPixel(void) const;
	UINT		BytesPerRow(void) const;
	UINT		BytesPerSlice(void) const;

	//------------------------------------------------------------------
	// @group Image ROI Size

	UINT		Width(void) const;				// Get the ROI Width
	UINT		Height(void) const;				// Get the ROI Height
	UINT		Depth(void) const;				// Get the ROI Depth
	UINT		SliceSize(void) const;			// Get the ROI Slice Size
	UINT		Size(void) const;				// Get the ROI Size

	UINT		BytesPerROIRow(void) const;
	UINT		BytesPerROISlice(void) const;

	
	//------------------------------------------------------------------
	// @group Image ROI coordinates

	UINT		Left(void) const;
	UINT		Top(void) const;
	UINT		Front(void) const;
	UINT		Right(void) const;
	UINT		Bottom(void) const;
	UINT		Back(void) const;

	bool		IfInside(int x, int y, int z=0) const;
	bool		IfInsideImage(int x, int y, int z=0) const;

	//------------------------------------------------------------------
	// @group Address operator for current ROI

	// the address of the first pixel of current ROI
	const PBYTE Origin(void) const;			
	PBYTE		Origin(void);
	const PBYTE End(void) const;			
	PBYTE		End(void);

	// address of a pixel given its coordinate in ROI
	const PBYTE PbPixel(UINT x, UINT y, UINT z=0) const;	
	PBYTE		PbPixel(UINT x, UINT y, UINT z=0);

	// address of a pixel given its coordinate of this image
	const PBYTE	IPbPixel(UINT x, UINT y, UINT z=0) const;
	PBYTE		IPbPixel(UINT x, UINT y, UINT z=0);

	// coordinate of a pixel given its address in the ROI
	void PixelOffset(PBYTE add, int& x, int& y, int& z) const;
	void IPixelOffset(PBYTE add, int& x, int& y, int& z) const;

	// the address of the first pixel in a slice of current ROI
	const PBYTE PbSlice(UINT z) const;
	PBYTE		PbSlice(UINT z);

	// Given a row containing pixels in the image, find the address of
	// the pixel in the first slice of the point in column zero of the
	// row specified.
	const PBYTE PbColZeroOfRow(UINT y) const;
	PBYTE		PbColZeroOfRow(UINT y);

	ImageRect	GetROI() const;
	bool		SetROI(UINT left, UINT right, UINT top, UINT bottom, UINT front=0, UINT back=0) const;
	bool		SetROI(const ImageRect& roi) const;
	UINT		NextSlice(void) const;
	UINT		PrevSlice(void) const;
	UINT		CurrentSlice(void) const;
	UINT		GotoSlice(UINT nSlice) const;
	UINT		DeleteSlice();
	ImageRect	AllSlice(void) const;
	ImageRect	FullSlice(void) const;
	ImageRect	FullROI() const;
	bool		IfFullROI() const;

	//------------------------------------------------------------------
	// @group Address operation for this Image

	// the address of the first pixel in a slice of this Image
	const PBYTE PbSliceImage(UINT z) const;
	PBYTE		PbSliceImage(UINT z);

	// Given a row containing pixels in the image, find the address
	// of the pixel in the first slice of the leftmost point in the
	// row specified.
	const PBYTE PbFirstPixelInRow(UINT y) const;
	PBYTE		PbFirstPixelInRow(UINT y);

	////////////////////////////////////////////////////////////////////////////
	//
	// @mfunc void | ZImageBase | GetOneSlice |
	//
	// Get one slice(a 2D image) from the 3D images.
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
	void		GetXYSlice(ZImageBase&, UINT z) const;
	void		GetYZSlice(ZImageBase&, UINT x) const;
	void		GetXZSlice(ZImageBase&, UINT y) const;
	void		GetYZVolume(ZImageBase& output) const;
	void		GetXZVolume(ZImageBase& output) const;

	//------------------------------------------------------------------
	// @group Other operation

	bool		FillBitmapInfo(BITMAPINFO** pbmpinfo) const;
	void		FillBitmapInfoHeader(BITMAPINFOHEADER& bmpinfoheader) const;

	virtual void SendPixelsTo(ZImageBase&) const = 0;
	virtual void CopyPixelsTo(ZImageBase&, float dstmin=0, float dstmax=0) const = 0;
#ifdef _WINDOWS
	virtual void CopyPixelsTo(ZGDImage&, bool full=true) const = 0;
#endif

	virtual void	Negative() = 0;
	virtual ZImageBase*	Reverse() const = 0;
	virtual void	operator += (float val) = 0;
	virtual void	operator -= (float val) = 0;
	virtual void	operator *= (float val) = 0;
	virtual void	operator /= (float val) = 0;

	virtual void	operator += (const ZImageBase&) = 0;
	virtual void	operator -= (const ZImageBase&) = 0;
	virtual void	operator *= (const ZImageBase&) = 0;
	virtual void	operator /= (const ZImageBase&) = 0;
	virtual void	operator &= (const ZImage<BYTE>&) = 0;
	virtual void	Diff(const ZImageBase& image) = 0;
	virtual void	Brighter(const ZImageBase& image) = 0;
	virtual void	Darker(const ZImageBase& image, bool noback=false) = 0;
	virtual void	Average(const ZImageBase& image, bool noback=false) = 0;
	virtual void	Mask(const ZImageBase& image) = 0;
	virtual void	Superimpose(const ZImageBase& image) = 0;
	virtual void	Mosaic(const ZImageBase& image, UINT size=20) = 0;

	virtual void	MinMax(float& min, float& max, bool noback=false, UINT channel = 0) const = 0;
	virtual float	Minimum(bool noback=false, UINT channel = 0) const = 0;
	virtual float	Maximum(bool noback=false, UINT channel = 0) const = 0;
	virtual UINT	Statistics(float& mean, float& stddev, bool noback=false, UINT channel=0) const = 0;
	virtual void	Histogram(UINT* hist, UINT size, bool noback=false, UINT channel=0, float min=0, float max=0) const = 0;
	virtual void	CGV(int& x, int& y, int& z) const = 0;

	virtual void	Inflate(int, int, int, int, bool bFill, float val) = 0;
	virtual void	Inflate(int, int, int, int, int, int, bool bFill, float val) = 0;
	virtual void	Inflate2D(int length, bool bFill=true, float val=0) = 0;
	virtual void	Inflate(int length, bool bFill=true, float val=0) = 0;

	virtual void	Mirror() = 0;
	virtual void	Flip() = 0;
	virtual void	Rotate(float angle, bool keepsize=true) = 0;
	virtual void	Rescale(UINT width, UINT height, int depth=-1) = 0;
	virtual void	Rescale(float scale) = 0;
	virtual void	IsotropicRescale(ZImageBase& iso) = 0;

	virtual float	GetIntensity(UINT x, UINT y, UINT z=0) const = 0;
	virtual ZRGB<float>	GetRGBIntensity(UINT x, UINT y, UINT z=0) const = 0;
	virtual void	SetIntensity(const ZRGB<float>& color, UINT x, UINT y, UINT z=0) = 0;
	virtual float	IGetIntensity(UINT x, UINT y, UINT z=0) const = 0;
	virtual ZRGB<float>	IGetRGBIntensity(UINT x, UINT y, UINT z=0) const = 0;
	virtual void	ISetIntensity(const ZRGB<float>& color, UINT x, UINT y, UINT z=0) = 0;
	virtual float	GrayBilinear(float x, float y, float z=0) const = 0;

	virtual void	RandomNoise(float min, float max) = 0;
	virtual void	GaussianNoise(float mean, float stddev) = 0;
	virtual void	RayleighNoise(float stddev) = 0;
	virtual void	MedianFilter(int radius) = 0;

	virtual void	NewDisplayRange(float contrast, float brightness) const = 0;
	virtual void	IdealDisplaySetting(float& contrast, float& brightness) const = 0;
	virtual void	AdjustContrast(float contrast, float brightness) = 0;

	virtual void	Boundary(ZImage<BYTE>& edge) const = 0;
};

// Inline functions


////////////////////////////////////////////////////////////////////////////
//
// @mfunc bool | ZImageBase | IfValid |
//
// Is this a valid image?  (Has it been initialized (not just ZImageBase*)?
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline bool ZImageBase::IfValid(void) const
{
	return m_ValidStamp == 123456789l;
}

inline void ZImageBase::CleanUp(void)
{
	delete []m_pImageBuf;
	m_epixfmt = epixfmtUnknown;
	m_pImageBuf = 0;
	m_pImageLim = 0; 
	m_nWidth = m_nHeight = m_nDepth = m_nSliceSize = m_nSize = 0;
	m_nBytesPerRow = m_nBytesPerSlice = m_cbPixel = 0;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc EPixFmt | ZImageBase | PixFmt |
//
// Returns a union of <t Pixel> values which give type
// information about the pixels used in this image.
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline EPixFmt ZImageBase::PixFmt(void) const
{
	return m_epixfmt;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc EPixFmt | ZImageBase | SetPixelDim |
// @mfunc EPixFmt | ZImageBase | GetPixelDim |
//
// Set real world measurements in mm of the pixel.
// Get real world measurements in mm of the pixel.
//
////////////////////////////////////////////////////////////////////////////
inline void	ZImageBase::SetPixelDim(float xdim, float ydim, float zdim) const
{
	m_ImageAtt.m_pixdim[0] = xdim;
	m_ImageAtt.m_pixdim[1] = ydim;
	m_ImageAtt.m_pixdim[2] = zdim;
}

inline void	ZImageBase::SetPixelDim(float dim) const
{
	m_ImageAtt.m_pixdim[0] = m_ImageAtt.m_pixdim[1] = m_ImageAtt.m_pixdim[2] = dim;
}

inline void	ZImageBase::GetPixelDim(float& xdim, float& ydim, float& zdim) const
{
	xdim = m_ImageAtt.m_pixdim[0];
	ydim = m_ImageAtt.m_pixdim[1];
	zdim = m_ImageAtt.m_pixdim[2];
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc PBYTE | ZImageBase | GetBuffer |
//
// Return the address of the image buffer
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline PBYTE ZImageBase::GetBuffer(void) const
{
	return m_pImageBuf;
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc PBYTE | ZImageBase | GetBuffer |
//
// Return the address of the last voxel;
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline PBYTE ZImageBase::GetLimit(void) const
{
	return m_pImageLim;
}


////////////////////////////////////////////////////////////////////////////
//
// Get the size of the this image.
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::ImageWidth(void) const
{
	return m_nWidth;
}

inline UINT ZImageBase::ImageHeight(void) const
{
	return m_nHeight;
}

inline UINT ZImageBase::ImageDepth(void) const
{
	return m_nDepth;
}

inline UINT ZImageBase::PixelsPerSlice(void) const
{
	return m_nSliceSize;
}

inline UINT ZImageBase::ImageSize(void) const
{
	return m_nSize;
}

inline UINT ZImageBase::MemorySize(void) const
{
	return m_nSize * m_cbPixel;
}



////////////////////////////////////////////////////////////////////////////
//
// @mfunc UINT | ZImageBase | BytesPerPixel |
//
// The size of image pixels in bytes.
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::BytesPerPixel(void) const
{
	return m_cbPixel;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc UINT | ZImageBase | BytesPerRow |
//
// The number of bytes of one image row
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::BytesPerRow(void) const
{
	return m_nWidth * m_cbPixel;
}



////////////////////////////////////////////////////////////////////////////
//
// @mfunc UINT | ZImageBase | BytesPerSlice |
//
// The number of bytes of one image slice
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::BytesPerSlice(void) const
{
	return m_nBytesPerSlice;
}


////////////////////////////////////////////////////////////////////////////
//
// Get the size of the image ROI.
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::Width(void) const
{
	return m_ImageROI.width;
}

inline UINT ZImageBase::Height(void) const
{
	return m_ImageROI.height;
}

inline UINT ZImageBase::Depth(void) const
{
	return m_ImageROI.depth;
}

inline UINT ZImageBase::SliceSize(void) const
{
	return m_ImageROI.SliceSize;
}

inline UINT ZImageBase::Size(void) const
{
	return m_ImageROI.Size;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc UINT | ZImageBase | BytesPerROIRow |
//
// The number of bytes in one row of the ROI.
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::BytesPerROIRow(void) const
{
	return m_ImageROI.BytesPerROIRow;
}



////////////////////////////////////////////////////////////////////////////
//
// @mfunc UINT | ZImageBase | BytesPerROISlice |
//
// The number of bytes in one slice of ROI
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::BytesPerROISlice(void) const
{
	return m_ImageROI.BytesPerROISlice;
}


////////////////////////////////////////////////////////////////////////////
//
// Get the coordinate of  the image ROI
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::Left(void) const
{
	return m_ImageROI.left;
}

inline UINT ZImageBase::Top(void) const
{
	return m_ImageROI.top;
}

inline UINT ZImageBase::Front(void) const
{
	return m_ImageROI.front;
}

inline UINT ZImageBase::Right(void) const
{
	return m_ImageROI.right;
}

inline UINT ZImageBase::Bottom(void) const
{
	return m_ImageROI.bottom;
}

inline UINT ZImageBase::Back(void) const
{
	return m_ImageROI.back;
}

inline ImageRect ZImageBase::GetROI() const
{
	return m_ImageROI;
}


inline bool ZImageBase::SetROI(const ImageRect& roi) const
{
	if(roi.left >= m_nWidth || roi.top >= m_nHeight || roi.front >= m_nDepth) return false;
	if(roi.right >= m_nWidth || roi.bottom >= m_nHeight || roi.back >= m_nDepth) return false;
	SetROI(roi.left, roi.right, roi.top, roi.bottom, roi.front, roi.back);
	return true;
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc bool | ZImageBase | IfInside |
//
// Is the specified point contained in the ROI rectangle?
//
// @parm UINT | x |
// The X coordinate of the point.
//
// @parm UINT | y |
// The Y coordinate of the point.
//
// @parm UINT | z |
// The Z coordinate of the point.
//
// @this const
//
////////////////////////////////////////////////////////////////////////////
inline bool ZImageBase::IfInside(int x, int y, int z) const
{
	return x >= 0 && y >= 0 && z >= 0 && x < int(Width()) && y < int(Height()) && z < int(Depth());	
}

inline bool ZImageBase::IfInsideImage(int x, int y, int z) const
{
	return x >= 0 && y >= 0 && z >= 0 && x < int(m_nWidth) && y < int(m_nHeight) && z < int(m_nDepth);
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | Origin |
//
// return the address of current original pixel.
//
// @syntax const PBYTE Origin(void) const;
// @syntax PBYTE Origin(void);
//
////////////////////////////////////////////////////////////////////////////
inline const PBYTE ZImageBase::Origin(void) const
{
	return m_ImageROI.begin;
}

inline PBYTE ZImageBase::Origin(void)
{
	return m_ImageROI.begin;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | End |
//
// return the address of current original pixel.
//
// @syntax const PBYTE Origin(void) const;
// @syntax PBYTE Origin(void);
//
////////////////////////////////////////////////////////////////////////////
inline const PBYTE ZImageBase::End(void) const
{
	return m_ImageROI.end;
}

inline PBYTE ZImageBase::End(void)
{
	return m_ImageROI.end;
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc | ZImageBase | PbPixel |
//
// Given the coordinates of a point in the image, return the reference
// of the pixel with the specified coordinates in ROI.
//
// @syntax const BYTE & PbPixel(int x, int y, int z) const;
// @syntax BYTE & PbPixel(int x, int y, int z);
//
// @parm UINT | x |
// Column of the pixel.
//
// @parm UINT | y |
// Row of the pixel.
//
// @parm UINT | z |
// Slice of the pixel.
//
////////////////////////////////////////////////////////////////////////////
inline const PBYTE ZImageBase::PbPixel(UINT x, UINT y, UINT z) const
{
	return Origin() + z * BytesPerSlice() + y * BytesPerRow() + x * BytesPerPixel();
}

inline PBYTE ZImageBase::PbPixel(UINT x, UINT y, UINT z)
{
	return Origin() + z * BytesPerSlice() + y * BytesPerRow() + x * BytesPerPixel();
}

inline const PBYTE ZImageBase::IPbPixel(UINT x, UINT y, UINT z) const
{
	return m_pImageBuf + z * BytesPerSlice() + y * BytesPerRow() + x * BytesPerPixel();
}

inline PBYTE ZImageBase::IPbPixel(UINT x, UINT y, UINT z)
{
	return m_pImageBuf + z * BytesPerSlice() + y * BytesPerRow() + x * BytesPerPixel();
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | PbSlice |
//
// Given a slice containing pixels in the ROI, find the address
// of the first pixel in the current ROI);
//
// @parm UINT | z |
// The slice of pixels.
//
////////////////////////////////////////////////////////////////////////////
inline const PBYTE ZImageBase::PbSlice(UINT z) const
{
	return Origin() + z * BytesPerSlice();
}

inline PBYTE ZImageBase::PbSlice(UINT z)
{
	return Origin() + z * BytesPerSlice();
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | PbColZeroOfRow |
//
// Given a row containing pixels in the image, find the address of
// the pixel in the first slice of the point in column zero of the
// row specified.  In the debug builds, this will assert that the
// row does contain pixels in the image, but this function
// may return an invalid pixel address if the row does not
// contain pixels in the image.
//
// @syntax const PBYTE PbColZeroOfRow(UINT y) const;
// @syntax PBYTE PbColZeroOfRow(UINT y);
//
// @parm UINT | y |
// The row of pixels.
//
////////////////////////////////////////////////////////////////////////////
inline const PBYTE ZImageBase::PbColZeroOfRow(UINT y) const
{
	return Origin() + y * BytesPerRow();
}

inline PBYTE ZImageBase::PbColZeroOfRow(UINT y)
{
	return Origin() + y * BytesPerRow();
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | PbSliceImage |
//
// Given a slice containing pixels in the image, find the address
// of the first pixel in the current slice of Image;
//
// @parm UINT | z |
// The slice of pixels.
//
////////////////////////////////////////////////////////////////////////////
inline const PBYTE ZImageBase::PbSliceImage(UINT z) const
{
	return m_pImageBuf + z * BytesPerSlice();
}

inline PBYTE ZImageBase::PbSliceImage(UINT z)
{
	return m_pImageBuf + z * BytesPerSlice();
}

////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | PbFirstPixelInRow |
//
// Given a row containing pixels in the image, find the address
// of the pixel in the first slice of the leftmost point in the
// row specified.  In the debug builds, this will assert that the
// row does contain pixels in the image, but this function
// may return an invalid pixel address if the row does not
// contain pixels in the image.
//
// @syntax const PBYTE PbFirstPixelInRow(UINT y) const;
// @syntax PBYTE PbFirstPixelInRow(UINT y);
//
// @parm UINT | y |
// The row of pixels.
//
////////////////////////////////////////////////////////////////////////////
inline const PBYTE ZImageBase::PbFirstPixelInRow(UINT y) const
{
	return m_pImageBuf + y * BytesPerRow();
}

inline PBYTE ZImageBase::PbFirstPixelInRow(UINT y)
{
	return m_pImageBuf + y * BytesPerRow();
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | NextSlice |
//
// Set the ROI to next slice; Return the current slice NO.
// 
// @syntax UINT NextSlice(void);
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::NextSlice(void) const
{
	UINT curslice = Front();
	if(m_ImageROI.front == m_nDepth-1) m_ImageROI.front = 0;
	else m_ImageROI.front++;
	SetROI(Left(), Right(), Top(), Bottom(), m_ImageROI.front, m_ImageROI.front);
	
	return curslice;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | PrevSlice |
//
// Set the ROI to previous slice; Return the current slice NO.
// 
// @syntax UINT PrevSlice(void);
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::PrevSlice(void) const
{
	UINT curslice = Front();
	if(m_ImageROI.front == 0) m_ImageROI.front = m_nDepth-1;
	else m_ImageROI.front--;
	SetROI(Left(), Right(), Top(), Bottom(), m_ImageROI.front, m_ImageROI.front);
	
	return curslice;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | CurrentSlice |
//
// Return the current slice NO.
// 
// @syntax UINT CurrentSlice()
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::CurrentSlice(void) const
{
	return Front();
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | GotoSlice |
//
// Goto slice of nslice; Return the current slice NO.
// 
// @syntax UINT GotoSlice(UINT z)
//
////////////////////////////////////////////////////////////////////////////
inline UINT ZImageBase::GotoSlice(UINT z) const
{
	UINT curslice = Front();
	if(z >= m_nDepth) z = m_nDepth - 1;
	SetROI(Left(), Right(), Top(), Bottom(), z, z);
	return curslice;
}


////////////////////////////////////////////////////////////////////////////
//
// @mfunc const PBYTE  | ZImageBase | AllSlice |
//
// Set the ROI to all slices; Return the current slice NO.
// 
// @syntax void AllSlice(void);
//
////////////////////////////////////////////////////////////////////////////
inline ImageRect ZImageBase::AllSlice(void) const
{
	ImageRect roi = m_ImageROI;
	SetROI(Left(), Right(), Top(), Bottom(), 0, m_nDepth-1);
	return roi;
}

inline ImageRect ZImageBase::FullSlice(void) const
{
	ImageRect roi = m_ImageROI;
	SetROI(0, m_nWidth-1, 0, m_nHeight-1, Front(), Back());
	return roi;
}

inline ImageRect ZImageBase::FullROI(void) const
{
	ImageRect roi = m_ImageROI;
	SetROI(0, m_nWidth-1, 0, m_nHeight-1, 0, m_nDepth-1);
	return roi;
}

inline bool	ZImageBase::IfFullROI() const
{
	return (m_nWidth == Width() && m_nHeight == Height() && m_nDepth == Depth());
}

#endif // _IMAGEBASE_H__
