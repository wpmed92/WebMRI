/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

/*
 * Colour map (palette) class
 */

#ifndef __PALETTE_H
#define __PALETTE_H

#include "mydefine.h"
#include "bitmap.h"

#ifndef RGB
#define RGB(r,g,b)      ((UINT)(((BYTE)(r)|((WORD)((BYTE)(g))<<8))|(((UINT)(BYTE)(b))<<16)))
#endif
#define GetR(rgb)		((BYTE) (rgb)) 
#define GetG(rgb)		((BYTE) (((WORD) (rgb)) >> 8)) 
#define GetB(rgb)		((BYTE) ((rgb) >> 16))


class ZImagePalette
{
    UINT			m_nColors;
    PALETTE_ENTRY*	m_pPalette;
	bool			m_bGray;

public:
	ZImagePalette();
	ZImagePalette(UINT n, const PALETTE_ENTRY* pal=0);
	ZImagePalette(UINT n, const PBYTE red, const PBYTE green, const PBYTE blue);
	ZImagePalette(UINT n, const PBYTE rgb);		//rgb should be n rgb bytes.
	ZImagePalette(const ZImagePalette&);

	~ZImagePalette(void);

	ZImagePalette& operator = (const ZImagePalette &);
	
	const PALETTE_ENTRY& operator[] (UINT index) const 
	{ 
		if(index >= m_nColors) index = m_nColors-1;
		return m_pPalette[index];
	}
	
	PALETTE_ENTRY&	operator[] (UINT index)
	{ 
		if(index >= m_nColors) index = m_nColors-1;
		return m_pPalette[index];
	}

	bool	Create(UINT n, const PALETTE_ENTRY* pal=0);
	bool	Create(UINT n, const PBYTE red, const PBYTE green, const PBYTE blue);
	bool	Create(UINT n, const PBYTE rgb);
	bool	CreateGrayPalette();
	UINT	NBColors() const { return m_nColors; }

	UINT	GetRGB(UINT index) const;
	UINT	GetNearestPaletteIndex(BYTE red, BYTE green, BYTE blue) const;
	bool	IfGray() const { return m_bGray; }

	PALETTE_ENTRY*	GetPalette() const { return m_pPalette; }
};

#endif
