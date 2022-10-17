/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declare the mother class for stream-based image IO classes.

#ifndef __IMAGEFORMAT_H
#define __IMAGEFORMAT_H

#include <fstream>
#include <cstdio>
#include <string>
#include "imagecore.h"
#include "imagereader.h"

class ZImageFile
{
protected:
	bool			m_bLittleEndian;
	bool			m_bReverse;

	UINT			m_nWidth, m_nHeight, m_nDepth, m_nNum;
	bool			m_rgbOrder;			//for RAW format

	std::string		m_filename;
	std::fstream	m_ImageStream;

public:
	ZImageFile();
	virtual ~ZImageFile() { Close(); }

	int				NumImages() const { return m_nNum; }

	virtual bool	Open(const char *fname, bool read);
	virtual bool	OpenText(const char *fname, bool read);
	void			SetFile(const char* fname) { m_filename = fname; }

	virtual void	Close();
	virtual bool	IsOpen() {return m_ImageStream.is_open()!=0;}
	virtual bool	Fail() {return m_ImageStream.fail()!=0; }
	virtual int		Read(void* buf, int num = 1);
	virtual int		Read2(void* buf, int num = 1);
	virtual int		Read4(void* buf, int num = 1);
	
	virtual bool	Write(char v, int num = 1);
	virtual bool	Write(short v, int num = 1);
	virtual bool	Write(int v, int num = 1);
	virtual bool	Write(float v, int num = 1);

	virtual bool	Write(void* buf, int num = 1);
	virtual bool	Write2(void* buf, int num = 1);
	virtual bool	Write4(void* buf, int num = 1);
	
	virtual bool	Seek(int offset, int origin);

	virtual ZImageBase* ReadFile() = 0;
	virtual bool	WriteFile(const ZImageBase&, bool intensityscale=false) = 0;
};

inline ZImageFile::ZImageFile()
{
	int val=1;
	char tmp, *ptr = (char *) &val;
	tmp = *ptr;
	m_bLittleEndian = (tmp==1) ? true : false;
	m_nWidth = m_nHeight = 0; 
	m_nNum = 1;
	m_rgbOrder = false;
	m_bReverse = false;
}

inline bool ZImageFile::Open(const char * fname, bool read)
{
	if(IsOpen()) return true;
	
	m_filename = std::string(fname);
	if(read) m_ImageStream.open(fname, std::ios::binary | std::ios::in);
	else m_ImageStream.open(fname, std::ios::binary | std::ios::out);

	return m_ImageStream.good();
}

inline bool ZImageFile::OpenText(const char * fname, bool read)
{
	if(IsOpen()) return true;
	
	m_filename = std::string(fname);
	if(read) m_ImageStream.open(fname, std::ios::in);
	else m_ImageStream.open(fname, std::ios::out);

	return m_ImageStream.good();
}

inline void ZImageFile::Close()
{
	if(IsOpen()) m_ImageStream.close();
}

inline int ZImageFile::Read(void* buf, int num)
{
	m_ImageStream.read((char*)buf, num);
	
	return m_ImageStream.gcount();
}

inline int ZImageFile::Read2(void* buf, int num)
{
	m_ImageStream.read((char*)buf, 2*num);
	if(m_bReverse) Reverse2(buf, num);

	return m_ImageStream.gcount() / 2;
}

inline int ZImageFile::Read4(void* buf, int num)
{
	m_ImageStream.read((char*)buf, 4*num);
	if(m_bReverse) Reverse4(buf, num);

	return m_ImageStream.gcount() / 4;
}

inline bool ZImageFile::Write(void* buf, int num)
{
	m_ImageStream.write((char*)buf, num);
	
	return m_ImageStream.good();
}

inline bool ZImageFile::Write2(void* buf, int num)
{
	if(m_bReverse) Reverse2(buf, num);
	m_ImageStream.write((char*)buf, 2*num);
	if(m_bReverse) Reverse2(buf, num);
	
	return m_ImageStream.good();
}

inline bool ZImageFile::Write4(void* buf, int num)
{
	if(m_bReverse) Reverse4(buf, num);
	m_ImageStream.write((char*)buf, 4*num);
	if(m_bReverse) Reverse4(buf, num);
	
	return m_ImageStream.good();
}

inline bool ZImageFile::Write(char v, int num)
{
	for(int i=0; i<num; i++) if(!(Write(&v))) return false;
	
	return true;
}

inline bool ZImageFile::Write(short v, int num)
{
	for(int i=0; i<num; i++) if(!(Write2((void*)(&v)))) return false;
	
	return true;
}

inline bool ZImageFile::Write(int v, int num)
{
	for(int i=0; i<num; i++) if(!(Write4((void*)(&v)))) return false;
	
	return true;
}

inline bool ZImageFile::Write(float v, int num)
{
	for(int i=0; i<num; i++) if(!(Write4((void*)(&v)))) return false;
	
	return true;
}

inline bool ZImageFile::Seek(int offset, int origin)
{
	switch(origin)
	{
	default:
	case SEEK_CUR: m_ImageStream.seekg(0, std::ios::cur); break;
	case SEEK_END: m_ImageStream.seekg(0, std::ios::end); break;
	case SEEK_SET: m_ImageStream.seekg(0, std::ios::beg); break;
	}
	
	return m_ImageStream.good();
}

ZImageBase* CreateImageFromIndex(PBYTE data, int Width, int Height, int aWidth, ZImagePalette& palette, int nbits);
bool TrueColorTo256Index(ZRGBByteImage& rgbimg, PBYTE data, ZImagePalette& palette);

#endif
