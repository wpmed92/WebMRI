/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file define the image iterator class which enumerates the pixels 
// in a rectangle in an image, from the top row to the bottom row and 
// from the leftmost column to he rightmost column in each row.  

#ifndef __ZITERATOR_H__
#define __ZITERATOR_H__

#ifdef _MSC_VER
#pragma warning( disable : 4786 )  
#endif

////////////////////////////////////////////////////////////////////////////
//  
// @class
//
// This object will enumerate the pixels in (all Slices or a specified
// Slice of) a rectangle in an image, from the top
// row to the bottom row and from the leftmost column to
// the rightmost column in each row.  This class knows the type of
// the pixels being enumerated, so its methods can return pointers
// and references to objects of the pixel type.
//
// @tcarg class | value_type |
// Type of the pixels in the image being enumerated.
//
// @base public | ZIterator
//
// @xref <c ZIterator>
//
////////////////////////////////////////////////////////////////////////////

template<class image_reference, class TPixel, class pointer, class reference> 
struct ZFastIterator
{
	// Constructor.
	ZFastIterator(image_reference image, pointer offset=0)
		: m_Image(image), m_pCur(offset>0 ? offset : image.ImageOrigin()), 
		m_pImageLim(image.ImageEnd()), ref(1)
	{ }
	
	virtual ~ZFastIterator()
	{ }

	virtual void reset(void)
	{
		m_pCur = m_Image.ImageOrigin();
		m_pImageLim = m_Image.ImageEnd();
	}

	bool more(void) const
	{
		return (m_pCur < m_pImageLim);
	}

	// Advance the current pixel.  If there are no more pixels to be
	// enumerated, the current pixel will not be valid after this method
	// is called.
	virtual void operator ++() { m_pCur ++; }
	virtual void operator ++(int) { operator ++(); }

	virtual void operator +=(int off) { m_pCur += off; }

	image_reference m_Image;
	pointer	m_pCur;
	pointer	m_pImageLim;
	int ref;
};

template<class image_reference, class TPixel, class pointer, class reference> 
struct ZROIIterator	: public ZFastIterator<image_reference, TPixel, pointer, reference>
{
	typedef ZFastIterator<image_reference, TPixel, pointer, reference> fast_iterator;

	// Constructor.
	ZROIIterator(image_reference image, pointer offset)
	  : fast_iterator(image, offset), 
		m_nPixelsPerRow(image.ImageWidth()),
		m_nPixelsPerSlice(image.PixelsPerSlice()),
		m_cRowNext(image.Width()), 
		m_cSliceNext((image.Height()-1) * image.ImageWidth()),
		m_pRowLim(this->m_pCur + m_cRowNext), 
		m_pSliceLim(m_pRowLim + m_cSliceNext) 
	{ 
		if(offset > 0)
		{
			int x, y, z;
			this->m_Image.PixelOffset(PBYTE(offset), x, y, z);
			if(x >= 0 && y >= 0 && z >= 0)
			{
				this->m_pCur = offset;
				m_pRowLim = const_cast<TPixel*>(this->m_Image.PPixel(0, y, z)) + m_cRowNext;
				m_pSliceLim = m_pRowLim + m_cSliceNext;
			}
			else this->m_pCur = m_pRowLim = m_pSliceLim = this->m_pImageLim;
		}
	}

	virtual ~ZROIIterator()
	{ }

	void reset(void)
	{
		m_nPixelsPerRow = this->m_Image.ImageWidth();
		m_nPixelsPerSlice = this->m_Image.PixelsPerSlice();
		m_cRowNext = this->m_Image.Width();
		m_cSliceNext = (this->m_Image.Height()-1) * this->m_Image.ImageWidth();
		fast_iterator::m_pCur = this->m_Image.ImageOrigin();

		m_pRowLim = this->m_pCur + m_cRowNext;
		m_pSliceLim = m_pRowLim + m_cSliceNext;
		this->m_pImageLim = this->m_Image.ImageEnd();
	}

	// Advance the current pixel.  If there are no more pixels to be
	// enumerated, the current pixel will not be valid after this method
	// is called.
	void operator ++() 
	{
		if ( ++this->m_pCur < m_pRowLim) return;
		else if ((m_pRowLim += m_nPixelsPerRow) <= m_pSliceLim)
		{
			// Start at the beginning of the next row.
			this->m_pCur = m_pRowLim - m_cRowNext;
		}
		else
		{
			m_pSliceLim += m_nPixelsPerSlice;
			m_pRowLim = m_pSliceLim - m_cSliceNext;
			this->m_pCur = m_pRowLim - m_cRowNext;
		}
	}

	void operator ++(int) { operator ++(); }

	void operator +=(int off) 
	{
		pointer n = this->m_pCur + off;
		if(n < m_pRowLim) { this->m_pCur = n; return; }
		
		int dif = n - m_pRowLim;

		if(dif < m_cRowNext)
		{
			if ((m_pRowLim += m_nPixelsPerRow) <= m_pSliceLim)
			{
				// Start at the beginning of the next row.
				this->m_pCur = m_pRowLim - m_cRowNext + dif;
			}
			else
			{
				m_pSliceLim += m_nPixelsPerSlice;
				m_pRowLim = m_pSliceLim - m_cSliceNext;
				this->m_pCur = m_pRowLim - m_cRowNext;
			}
		}
		else
		{
			int x, y, z;
			this->m_Image.PixelOffset(PBYTE(n), x, y, z);
			if(x >= 0 && y >= 0 && z >= 0)
			{
				this->m_pCur = n;
				m_pRowLim = const_cast<TPixel*>(this->m_Image.PPixel(0, y, z)) + m_cRowNext;
				m_pSliceLim = m_pRowLim + m_cSliceNext;
			}
			else this->m_pCur = m_pRowLim = m_pSliceLim = this->m_pImageLim;
		}
	}

	int		m_nPixelsPerRow;
	int		m_nPixelsPerSlice;

	// The difference between the address of the first pixel enumerated
	// in a row and the address of the limit pixel in a row.
	int		m_cRowNext;
	int		m_cSliceNext;

	pointer	m_pRowLim;
	pointer	m_pSliceLim;
};

template<class image_reference, class TPixel, class pointer=TPixel*, class reference=TPixel&> 
class ZIterator
{
	typedef ZFastIterator<image_reference, TPixel, pointer, reference> fast_iterator;
	typedef ZROIIterator<image_reference, TPixel, pointer, reference> roi_iterator;
	fast_iterator*  m_pIterator;

public:
	typedef ZIterator<image_reference, TPixel, pointer, reference> image_iterator;

	ZIterator(image_reference image, pointer offset=0)
	{
		if(image.Width() == image.ImageWidth() && 
		   (image.Height() == image.ImageHeight() || image.Depth() == 1))
		   m_pIterator = new fast_iterator(image, offset);
		else
			m_pIterator = new roi_iterator(image, offset);
	}

	ZIterator(const ZIterator& x)
	{
		x.m_pIterator->ref ++;
		m_pIterator = x.m_pIterator;
	}

	~ZIterator() { if(--m_pIterator->ref==0) delete m_pIterator; }

	ZIterator& operator = (const ZIterator& x)
	{
		x.m_pIterator->ref ++;
		if(--m_pIterator->ref==0) delete m_pIterator;
		m_pIterator = x.m_pIterator;
		return *this;
	}

	void reset(void) const { m_pIterator->reset(); }
	bool more(void) const { return m_pIterator->more(); }

	// Advance the current pixel.  If there are no more pixels to be
	// enumerated, the current pixel will not be valid after this method
	// is called.
	void operator ++ () { m_pIterator->operator ++(); }
	void operator ++ (int) { m_pIterator->operator ++(); }
	void operator += (int off) { m_pIterator->operator +=(off); }

	// Get a reference to the current pixel.
	reference operator *() { return *(m_pIterator->m_pCur); }

	// Get a pointer to the current pixel.
	pointer operator ->() { return m_pIterator->m_pCur; }

	operator pointer () { return m_pIterator->m_pCur; }

	friend bool operator == (const image_iterator& x, const image_iterator& y)
	{
		return x.m_pIterator->m_pCur == y.m_pIterator->m_pCur;
	}

	friend bool operator != (const image_iterator& x, const image_iterator& y)
	{
		return x.m_pIterator->m_pCur != y.m_pIterator->m_pCur;
	}
};

#endif // __ZITERATOR_H__
