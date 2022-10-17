/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "imagecore.h"
using namespace std;

ZImageBase* CreateImageFromIndex(PBYTE data, int Width, int Height, int aWidth, ZImagePalette& palette, int nbits)
{
	ZImageBase*	image = 0;
	if(palette.IfGray()) image = new ZGrayByteImage(Width, Height);
	else image = new ZRGBByteImage(Width, Height);
	if (image == 0) return 0;

	PBYTE	src = data;

	int	j;
	ZRGB<BYTE>*	dst = (ZRGB<BYTE>*)(image->GetBuffer());
	PBYTE	graydst = image->GetBuffer();

	switch(nbits)
	{
	case 1:
		for(j=0; j<Height; j++, src+=aWidth)
		{
			int num = 0;
			for(int i=0; i<aWidth; i++)
			{
				BYTE byte = src[i];
				for(int k=0; k<8; k++, num++, graydst++)
				{
					if(num >= Width) break;
					if(byte&(0x80>>k)) *graydst = 255;
					else *graydst = 0;
				}
			}
		}
		break;
		
	case 4:
		for(j=0; j<Height; j++, src+=aWidth, dst+=Width, graydst+=Width)
		{
			int num = 0;
			for(int i=0; i<(Width+1)/2; i++)
			{
				BYTE byte = src[i];
				for(int k=0; k<=4; k+=4, num++)
				{
					if(num >= Width) break;
					int index = ((0xf0 >> k) & byte) >> (4-k);
					int ii = (i<<1) +(k>>2);
					
					UINT rgb = palette.GetRGB(index);
					if(palette.IfGray()) graydst[ii] = GetR(rgb);
					else
					{
						dst[ii].r = GetR(rgb);
						dst[ii].g = GetG(rgb);
						dst[ii].b = GetB(rgb);
					}
				}
			}
		}
		break;

	case 8:
		for(j=0; j<Height; j++, src+=aWidth)
		{
			for(int i=0; i<Width; i++, dst++, graydst++)
			{
				int index = src[i];
				UINT rgb = palette.GetRGB(index);
				if(palette.IfGray()) *graydst = GetR(rgb);
				else
				{
					dst->r = GetR(rgb);
					dst->g = GetG(rgb);
					dst->b = GetB(rgb);
				}
			}
		}
		break;

	case 24:
	case 32:
		for(j=0; j<Height; j++, src+=aWidth)
		{
			int nBytesPerPixel = nbits / 8;
			for(int i=0; i<Width*nBytesPerPixel; i+=nBytesPerPixel, dst++)
			{
				dst->r = src[i];
				dst->g = src[i+1];
				dst->b = src[i+2];
			}
		}
	}

	return image;
}

ZImageBase* CreateImageFromDIB(PBYTE data)
{
	PALETTE_ENTRY*	pal=0;
	UINT	Width, Height, aWidth, numColors;
	ZImagePalette palette;

	BITMAPINFOHEADER* bmpInfoHeader = (BITMAPINFOHEADER*)data;
	Width = bmpInfoHeader->biWidth;
	Height = bmpInfoHeader->biHeight;
	aWidth = (Width * bmpInfoHeader->biBitCount + 31) /32 * 4;
	numColors = bmpInfoHeader->biBitCount > 8 ? 0 : short(1) << bmpInfoHeader->biBitCount;
	
	if(bmpInfoHeader->biBitCount <= 8) pal = (PALETTE_ENTRY*)(data + bmpInfoHeader->biSize);
	palette.Create(bmpInfoHeader->biClrUsed, pal);

	PBYTE ptr1 = data+bmpInfoHeader->biSize+numColors*sizeof(PALETTE_ENTRY);

	vector<BYTE> tmp(aWidth*Height);
	PBYTE ptr2 = &(*tmp.begin())+(Height-1)*aWidth;

	for(UINT i=0; i<Height; i++, ptr1+=aWidth, ptr2-=aWidth) memcpy(ptr2, ptr1, aWidth);

	if(numColors == 0)
	{
		PBYTE ptr = &(*tmp.begin());
		int BytesPerPixels = bmpInfoHeader->biBitCount == 24 ? 3 : 4;
		UINT width = 3 * Width;
		for(UINT j=0; j<Height; j++, ptr+=aWidth)
		for(UINT i=0; i<width; i+=BytesPerPixels) swap(ptr[i], ptr[i+2]);
	}

	return CreateImageFromIndex(&(*tmp.begin()), Width, Height, aWidth, palette, bmpInfoHeader->biBitCount);
}

bool TrueColorTo256Index(ZRGBByteImage& rgbimg, PBYTE pBuf, ZImagePalette& palette)
{
	int   index;
	int   *colors, max, mix1, mix2;

	int  dist;
	int  i, j, tmp, istart, iend;
	int  *ind;
	BYTE *rgb;

	if(pBuf == 0) return false;
	
	if ((rgb = new BYTE[4096+sizeof(int)*4096+sizeof(int)*256]) == 0) 
		return false;    // can not allocate temp buffer
	colors = (int *)(rgb+4096);
	ind = colors + 4096;

	for(i=0; i<4096; i++)  // Initialize the colors  statical table
	  colors[i]=0;
	
	ZRGBByteImage::const_iterator p_src(rgbimg);

	for(; p_src.more(); p_src++)	// Count the colors
	{
		 int tmp = ((WORD(p_src->r) & 0xf0) << 4) 
				 | (WORD(p_src->g) & 0xf8 ) 
				 | ((WORD(p_src->b) & 0xe0 ) >> 5);  // r : g : b = 4:5:3  , 4096 colors selected
		 
		 colors[tmp]++;
	}

	istart = 0;
	iend = 256;
	
	palette.Create(256);
	for(i=istart; i<iend; i++)         // Select palette colors
	{
		max=0;   index = i;
		for (j=0; j<4096; j++)
		if (colors[j]>max) 
		{
			max=colors[j];
			index=j;
		}
		ind[i] = index;
		colors[index] = -colors[index];                          // negtive it to aviod another selection
		rgb[index] = (BYTE)i;
		palette[i].peRed = BYTE((index&0x0f00)>>4);     // set it's red
		palette[i].peGreen = BYTE((index&0x0f8));         // set it's green
		palette[i].peBlue = BYTE((index&0x07)<<5);       // set it's blue
	}
	
	int   red, grn, blu;
	for (i=0; i<4096; i++)		// Map other colors to palette using distance
	{     
		tmp=0;                        // no selection done
		if(colors[i] > 0)			// if it is an selected color or an zero counted one
		{        
			dist=65535;                      // maximum distance
			for (j=istart; j<iend; j++) 
			{
				red = (palette[j].peRed>>3);
				grn = (palette[j].peGreen>>3);
				blu = (palette[j].peBlue>>3);
				red -= ((i&0x0f00)>>7);
				grn -= ((i&0x00f8)>>3);
				blu -= ((i&0x0007)<<2);
				WORD value = WORD(red*red+grn*grn+blu*blu);   // calculate the distance
				if (value<dist)
				{
					dist=value;   // an better match
					tmp=j;
				}
			}
		}
		
		if(colors[i] >= 0) rgb[i]=(BYTE)tmp;      // store the matching relation to that palette
		if(colors[i] > 0) 
		{
			mix1 = -colors[ind[tmp]];
			mix2 = colors[i];
			// mix red
			red =   (int)(( (UINT)palette[tmp].peRed*mix1
				+(UINT)((i&0x0f00)>>4)*mix2)    // mix this two colors by a weighted average
				 /(mix1+mix2));
			if(red > 255) red = 255;                // tructuate
			palette[tmp].peRed =(BYTE)red;
			// mix green
			grn =  (int)(( (UINT)palette[tmp].peGreen*mix1
			   +(UINT)(i&0x0f8)*mix2)   // mix this two colors by a weighted average
			   /(mix1+mix2));
			if(grn > 255) grn = 255;                // tructuate
			palette[tmp].peGreen =(BYTE)grn;
			// mix blue
			blu =   (int)(( (UINT)palette[tmp].peBlue*mix1
				+(UINT)((i&0x07)<<5)*mix2)      // mix this two colors by a weighted average
				/(mix1+mix2));
			if(blu > 255) blu = 255;                // tructuate
			palette[tmp].peBlue =(BYTE)blu;
			// calculate mixed numbers
			colors[ind[tmp]] -= mix2;       // mixed numbers (Negtive)
		}
	}

	// look up table for 8 bit color data

	PBYTE p_dst = pBuf;
	for(p_src.reset(); p_src.more(); p_src++, p_dst++)	// Count the colors
	{
		 tmp=((p_src->r & 0xf0)<<4) | (p_src->g & 0xf8) | ((p_src->b & 0xe0)>>5);  // r : g : b = 4:5:3  , 4096 colors selected
		 *p_dst = rgb[tmp];
	}
	
    for(i=istart; i<iend; i++) 
	{
	   palette[i].peRed = palette[i].peRed ? (palette[i].peRed |0x0f) : 0;
	   palette[i].peGreen = palette[i].peGreen ? (palette[i].peGreen |0x07) : 0;
	   palette[i].peBlue = palette[i].peBlue ? (palette[i].peBlue |0x1f) : 0;
	}
	
	delete []rgb;

	return true;
}
