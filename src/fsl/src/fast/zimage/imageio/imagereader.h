/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declare the interface class for the general image IO.

#ifndef __IMAGEREADER_H
#define __IMAGEREADER_H

#include "imagecore.h"
#include <string>

enum	ImageType {UNKNOWN=0, AVW=200};

bool WriteImage(const ZImageBase& refimage, std::fstream& file, bool reverse=false);

class ZImageReader
{
	std::string		m_FileName;
	ImageType		m_nImageType;
protected:
	int				m_nAVWVolume;					// which volume to read in a 4D AVW.

public:

	ZImageReader ();
	ZImageReader (const char* fname, const char* ftype=NULL);

	static ImageType	Recognize(const char* fname, const char* ftype=NULL);
	ImageType	SetFile(const char* fname, const char* ftype=NULL);
	ImageType	SetFile(const char* fname, ImageType ntype);

	ZImageBase*	ReadFile();
	bool		ReadFile(ZImageBase& Image);
	bool		ReadFile(std::vector<ZImageBase*>& images);
	bool		WriteFile(const ZImageBase&, bool intensityscale=false);

	void		SetAVWVolume(int vol) { m_nAVWVolume = vol; }
};

ZImageBase* CreateImageFromDIB(PBYTE data);

#endif
