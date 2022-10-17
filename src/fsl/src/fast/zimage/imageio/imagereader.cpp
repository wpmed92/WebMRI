/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <cctype>
#include "imagereader.h"
using namespace std;

#include "fileavw.h"

bool WriteImage(const ZImageBase& refimage, fstream& file, bool reverse)
{
	if(!file.is_open()) return false;
	if(!refimage.IfValid() || refimage.GetBuffer() == 0) return false;

	UINT k,j;
	switch(refimage.BytesPerPixel())
	{
	default:
		for(k=0; k<refimage.Depth(); k++)
		for(j=0; j<refimage.Height(); j++)
		{
			PBYTE ptr = refimage.PbPixel(0,j,k);
			file.write((char*)ptr, refimage.BytesPerROIRow());
		}
		break;
	
	case 2:
		for(k=0; k<refimage.Depth(); k++)
		for(j=0; j<refimage.Height(); j++)
		{
			PBYTE ptr = refimage.PbPixel(0,j,k);
			if(reverse) Reverse2(ptr, refimage.Width());
			file.write((char*)ptr, refimage.BytesPerROIRow());
			if(reverse) Reverse2(ptr, refimage.Width());
		}
		break;
	
	case 4:
		for(k=0; k<refimage.Depth(); k++)
		for(j=0; j<refimage.Height(); j++)
		{
			PBYTE ptr = refimage.PbPixel(0,j,k);
			if(reverse) Reverse4(ptr, refimage.Width());
			file.write((char*)ptr, refimage.BytesPerROIRow());
			if(reverse) Reverse4(ptr, refimage.Width());
		}
		break;
	}

	return true;
}

/////////////////////////////////////////////////////////////////////////////
// ZImageReader
ZImageReader::ZImageReader() 
	:	m_nImageType(UNKNOWN), m_nAVWVolume(0) 
{
}

ZImageReader::ZImageReader(const char* fname, const char* ftype)
	:	m_nImageType(UNKNOWN), m_nAVWVolume(0)
{
	m_FileName = fname;

	m_nImageType = Recognize(fname, ftype);
}

ImageType ZImageReader::Recognize(const char* fname, const char* ftype)
{
	char	type[10];
	if(ftype != NULL) strncpy(type, ftype, 5);
	else
	{
		char*	PosOfLastDot = strrchr(fname, '.');
		if(PosOfLastDot == NULL) return UNKNOWN;
		strncpy(type, PosOfLastDot+1, 5);
	}

	for(UINT i=0; i<strlen(type); i++) type[i] = char(toupper(int(type[i])));

	if(strcmp(type, "HDR")==0) return AVW;
	if(strcmp(type, "IMG")==0) return AVW;

	ZWarning("ImageReader", "%s -- Unknown Image Format %s!", fname, type);

	return UNKNOWN;
}

ImageType  ZImageReader::SetFile(const char* fname, const char* ftype)
{
	m_FileName = fname;

	m_nImageType = Recognize(fname, ftype);
	return m_nImageType;
}

ImageType  ZImageReader::SetFile(const char* fname, ImageType ntype)
{
	if(ntype == UNKNOWN) m_nImageType = Recognize(fname);
	else m_nImageType = ntype;

	m_FileName = fname;

	return m_nImageType;
}

ZImageBase* ZImageReader::ReadFile()
{
	if(m_nImageType == UNKNOWN) return 0;

	ZImageFile*	img=0;

	switch(m_nImageType)
	{
	case AVW:
		img = new ZAVWImage;
		((ZAVWImage*)(img))->SetVolume(m_nAVWVolume);
		break;

	default:
		ZWarning("ImageReader", "%s -- Unknown Image Format!", m_FileName.c_str());
		return 0;
	}

	if(!img->Open(m_FileName.c_str(), true))
	{
		delete img;
		return 0;
	}

	ZImageBase* image = img->ReadFile();

	delete img;
	return image;
}

bool ZImageReader::ReadFile(ZImageBase& Image)
{
	if(m_nImageType == UNKNOWN)
	{
		ZWarning("ImageReader", "%s -- Unknown image format!", m_FileName.c_str());
		return false;
	}

	if(!Image.IfValid())
	{
		ZWarning("ImageReader", "%s -- Not a valid image object!", m_FileName.c_str());
		return false;
	}

	ZImageFile*	img=0;

	switch(m_nImageType)
	{
	case AVW:
		img = new ZAVWImage;
		((ZAVWImage*)(img))->SetVolume(m_nAVWVolume);
		break;

	default:
		ZWarning("ImageReader", "%s -- Unknown Image Format!", m_FileName.c_str());
		return 0;
	}
	
	if(!img->Open(m_FileName.c_str(), true))
	{
		delete img;
		return false;
	}

	ZImageBase* image = img->ReadFile();
	if(image==0) { delete img; return false; }

	delete img;

	image->CopyPixelsTo(Image);
	delete image;

	return true;
}

bool ZImageReader::ReadFile(vector<ZImageBase*>& images)
{
	if(m_nImageType != AVW)
	{
		ZWarning("ImageReader", "%s -- Only AVW format supports 4D images!", m_FileName.c_str());
		return false;
	}

	ZAVWImage* img = new ZAVWImage;
	img->Open(m_FileName.c_str(), true);
	ZImageBase* image = img->ReadFile();
	int num = img->NumImages();
	images.resize(num);
	images[0] = image;
	for(int i=1; i<num; i++) 
	{
		img->SetVolume(i);
		if((images[i] = img->ReadFile()) == 0) return false;
	}
	return true;
}

bool ZImageReader::WriteFile(const ZImageBase& image, bool intensityscale)
{
	ZImageFile*	img=0;

	switch(m_nImageType)
	{
	case AVW:
		img = new ZAVWImage;
		break;

	default:
		ZWarning("ImageReader", "%s -- Unknown Image Format!", m_FileName.c_str());
		return false;
	}

	bool ret = true;
	if(!img->Open(m_FileName.c_str(), false))
	{
		delete img;
		return false;
	}
	ret = img->WriteFile(image, intensityscale);

	delete img;
	return ret;
}

