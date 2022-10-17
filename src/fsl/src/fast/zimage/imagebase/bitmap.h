/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef __BITMAP_H
#define __BITMAP_H

struct PALETTE_ENTRY 
{
    BYTE        peBlue;
    BYTE        peGreen;
    BYTE        peRed;
    BYTE        peFlags;
};

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#if !defined(_WINDOWS) && !defined(_WINGDI_)

struct BITMAPINFOHEADER
{
    UINT       biSize;
    int        biWidth;
    int        biHeight;
    WORD       biPlanes;
    WORD       biBitCount;
    UINT       biCompression;
    UINT       biSizeImage;
    int        biXPelsPerMeter;
    int        biYPelsPerMeter;
    UINT       biClrUsed;
    UINT       biClrImportant;
};

struct BITMAPINFO 
{
    BITMAPINFOHEADER  bmiHeader;
    PALETTE_ENTRY	  bmiColors[1];
};

#endif

#endif
