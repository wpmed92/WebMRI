/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declare the mother class of the file handle(FILE*) based image 
// IO class used for AVW format.

#ifndef __IMAGE_HANDLE_H
#define __IMAGE_HANDLE_H

#include <cstdio>
#include "format.h"

class ZImageFileHandle	:	public ZImageFile
{
protected:
	FILE*			m_ImageHandle;

public:
	ZImageFileHandle() : ZImageFile(), m_ImageHandle(0) {}
	~ZImageFileHandle() { Close(); }

	virtual bool	Open(const char *fname, bool read);
	virtual void	Close();
	virtual bool	IsOpen() {return m_ImageHandle!=0;}
	virtual bool	Fail() {return (m_ImageHandle==0) || ferror(m_ImageHandle)!=0; }
	virtual int		Read(void* buf, int num = 1);
	virtual int		Read2(void* buf, int num = 1);
	virtual int		Read4(void* buf, int num = 1);
	
	virtual bool	Write(void* buf, int num = 1);
	virtual bool	Write2(void* buf, int num = 1);
	virtual bool	Write4(void* buf, int num = 1);

	virtual bool	Seek(int offset, int origin);

	virtual ZImageBase* ReadFile() = 0;
	virtual bool	WriteFile(const ZImageBase&, bool intensityscale=false) = 0;
};

inline bool ZImageFileHandle::Open(const char * fname, bool read)
{
	if(IsOpen()) return true;

	m_filename = std::string(fname);

	if(read) m_ImageHandle = fopen(fname, "rb");
	else m_ImageHandle = fopen(fname, "wb");
	
	return m_ImageHandle != 0;
}

inline void ZImageFileHandle::Close()
{
	if(IsOpen()) fclose(m_ImageHandle);
	m_ImageHandle = 0;
}

inline int ZImageFileHandle::Read(void* buf, int num)
{
	int n = fread(buf, 1, num, m_ImageHandle);
	
	return n;
}

inline int ZImageFileHandle::Read2(void* buf, int num)
{
	int n = fread(buf, 2, num, m_ImageHandle);

	if(m_bReverse) Reverse2(buf, num);

	return n;
}

inline int ZImageFileHandle::Read4(void* buf, int num)
{
	int n = fread(buf, 4, num, m_ImageHandle);
	
	if(m_bReverse) Reverse4(buf, num);

	return n;
}

inline bool ZImageFileHandle::Write(void* buf, int num)
{
	int n = fwrite(buf, 1, num, m_ImageHandle);
	
	return n != num;
}

inline bool ZImageFileHandle::Write2(void* buf, int num)
{
	if(m_bReverse) Reverse2(buf, num);
	int n = fwrite(buf, 2, num, m_ImageHandle);
	if(m_bReverse) Reverse2(buf, num);
	
	return n != num;
}

inline bool ZImageFileHandle::Write4(void* buf, int num)
{
	if(m_bReverse) Reverse4(buf, num);
	int n = fwrite(buf, 4, num, m_ImageHandle);
	if(m_bReverse) Reverse4(buf, num);
	
	return n != num;
}

inline bool ZImageFileHandle::Seek(int offset, int origin)
{
	return fseek(m_ImageHandle, offset, origin) == 0;
}

#endif
