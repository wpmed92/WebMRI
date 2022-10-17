/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

/*
 * Colour map class
 *
 */

#include <memory.h>
#include "palette.h"

ZImagePalette::ZImagePalette() 
	: m_nColors(0), m_pPalette(0), m_bGray(false)
{
}

ZImagePalette::ZImagePalette(UINT n, const PALETTE_ENTRY* pal)
	: m_nColors(0), m_pPalette(0), m_bGray(false)
{
	Create(n, pal);
}

ZImagePalette::ZImagePalette(UINT n, const PBYTE red, const PBYTE green, const PBYTE blue)
	: m_nColors(0), m_pPalette(0), m_bGray(false)
{
	Create(n, red, green, blue);
}

ZImagePalette::ZImagePalette(UINT n, const PBYTE rgb)		//rgb should be n rgb bytes.
	: m_nColors(0), m_pPalette(0), m_bGray(false)
{
	Create(n, rgb);
}

ZImagePalette::ZImagePalette(const ZImagePalette& pal)
	: m_nColors(pal.m_nColors), m_pPalette(pal.m_pPalette), m_bGray(pal.m_bGray)
{
	if(m_nColors > 0)
	{
		m_pPalette = new PALETTE_ENTRY[m_nColors];
		memcpy(m_pPalette, pal.m_pPalette, sizeof(PALETTE_ENTRY) * m_nColors);
	}
}

ZImagePalette::~ZImagePalette(void)
{
	delete []m_pPalette;
}

ZImagePalette& ZImagePalette::operator= (const ZImagePalette & pal)
{
	UINT count = pal.m_nColors;
	if(count == 0) return *this;

	m_nColors = count;
	m_bGray = pal.m_bGray;
	delete []m_pPalette;
	m_pPalette = new PALETTE_ENTRY[m_nColors];
	memcpy(m_pPalette, pal.m_pPalette, sizeof(PALETTE_ENTRY) * m_nColors);

	return *this;
}

bool ZImagePalette::Create(UINT n, const PALETTE_ENTRY* pal)
{
	if(n==0) return false;

	delete []m_pPalette;

	m_pPalette = new PALETTE_ENTRY[n];
	m_nColors = n;
	m_bGray = false;

	if(pal == 0) memset(m_pPalette, 0, sizeof(PALETTE_ENTRY)*n);
	else
	{
		memcpy(m_pPalette, pal, sizeof(PALETTE_ENTRY) * m_nColors);
		m_bGray = true;
		for (UINT i=0; i<n; i++)
		{
			if(m_pPalette[i].peRed != m_pPalette[i].peGreen || m_pPalette[i].peRed != m_pPalette[i].peBlue)
			{
				m_bGray = false;
				break;
			}
		}
	}
	return true;
}

bool ZImagePalette::Create(UINT n, const PBYTE red, const PBYTE green, const PBYTE blue)
{
	if(n==0) return false;

	delete []m_pPalette;

	m_pPalette = new PALETTE_ENTRY[n];
	m_nColors = n;
	m_bGray = true;

	for (UINT i = 0; i < n; i ++)
	{
		m_pPalette[i].peRed = red[i];
		m_pPalette[i].peGreen = green[i];
		m_pPalette[i].peBlue = blue[i];
		m_pPalette[i].peFlags = 0;
		if(m_bGray)
			if(red[i]!=green[i] || red[i]!=blue[i]) m_bGray = false;
	}
	
	return true;
}

bool ZImagePalette::Create(UINT n, const PBYTE rgb)
{
	if(n==0) return false;

	delete []m_pPalette;

	m_pPalette = new PALETTE_ENTRY[n];
	m_nColors = n;
	m_bGray = true;

	for (UINT i=0; i<n; i++)
	{
		m_pPalette[i].peRed = rgb[i*3];
		m_pPalette[i].peGreen = rgb[i*3+1];
		m_pPalette[i].peBlue = rgb[i*3+2];
		m_pPalette[i].peFlags = 0;
		if(m_bGray)
			if(rgb[i*3]!=rgb[i*3+1] || rgb[i*3]!=rgb[i*3+2]) m_bGray = false;
	}
	
	return true;
}

bool ZImagePalette::CreateGrayPalette()
{
	delete []m_pPalette;

	m_nColors = 256;
	m_pPalette = new PALETTE_ENTRY[256];
	m_bGray = true;

	for (UINT i = 0; i < 256; i ++)
	{
		m_pPalette[i].peRed = (BYTE)i;
		m_pPalette[i].peGreen = (BYTE)i;
		m_pPalette[i].peBlue = (BYTE)i;
		m_pPalette[i].peFlags = 0;
	}
	return true;
}

UINT ZImagePalette::GetNearestPaletteIndex(BYTE red, BYTE green, BYTE blue) const
{
	UINT	dis1, dis = 0xFFFFFFFF; 
	UINT	index=0;

	for(UINT i=0; i<m_nColors; i++)
	{
		dis1 = (UINT)(red-m_pPalette[i].peRed)*(red-m_pPalette[i].peRed) 
			 + (UINT)(green-m_pPalette[i].peGreen)*(green-m_pPalette[i].peGreen) 
			 + (UINT)(blue-m_pPalette[i].peBlue)*(blue-m_pPalette[i].peBlue);

		if(dis>dis1) 
		{
			dis = dis1;
			index = i;
		}
		if(dis == 0) break;
	}
	return index;
}

UINT ZImagePalette::GetRGB(UINT index) const
{
	if (index > 255) return 0;

	return RGB(m_pPalette[index].peRed,m_pPalette[index].peGreen,m_pPalette[index].peBlue);
}
