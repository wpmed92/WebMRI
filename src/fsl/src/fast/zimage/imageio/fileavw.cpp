/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <cstdio>
#include <cctype>
#include "fileavw.h"
using namespace std;

/////////////////////////////////////////////////////////////////////
/// ZAVWImage
bool ZAVWImage::Open(const char *fname, bool read)
{
	if(IsOpen()) return true;

	m_filename = fname;
	m_nVolume = 0;
	string imgfname, hdrfname;

	char* ext=strrchr(fname, '.');
	if(ext == 0)
	{
		hdrfname = fname; hdrfname += ".hdr"; 
		imgfname = fname; imgfname += ".img"; 
	}
	else
	{
		string fext = ext+1;
		for(UINT i=0; i<fext.size(); i++) fext[i] = char(tolower(fext[i]));

		if(fext == "hdr")
		{
			hdrfname = fname; imgfname = fname; 
			imgfname.replace(imgfname.rfind("hdr"), string::npos, "img");
		}
		else if(fext == "img")
		{
			hdrfname = fname; imgfname = fname; 
			hdrfname.replace(hdrfname.rfind("img"), string::npos, "hdr");
		}
		else
		{
			hdrfname = fname; hdrfname += ".hdr"; 
			imgfname = fname; imgfname += ".img"; 
		}
	}

	if(read) m_HeaderStream.open(hdrfname.c_str(), ios::binary | ios::in);
	else m_HeaderStream.open(hdrfname.c_str(), ios::binary | ios::out);
	if(m_HeaderStream.fail())
	{
		ZWarning("ZAVWImage::Open", "%s -- Can not open header file!", m_filename.c_str());
		return false;
	}

	if(read) m_ImageStream.open(imgfname.c_str(), ios::binary | ios::in);
	else m_ImageStream.open(imgfname.c_str(), ios::binary | ios::out);
	if(m_ImageStream.fail())
	{
		ZWarning("ZAVWImage::Open", "%s -- Can not open img file!", m_filename.c_str());
		return false;
	}

	return true;
}

void ZAVWImage::Close()
{
	if(IsOpen())
	{
		m_HeaderStream.close();
		m_ImageStream.close();
	}
}

bool ZAVWImage::ReadHeader()
{
	m_HeaderStream.read((char*)(&m_Header), sizeof(struct dsr));
	if(m_HeaderStream.fail())
	{
		ZWarning("AVWReader", "%s -- Can not read header file!", m_filename.c_str());
		return false;
	}

	m_bReverse = (m_Header.hk.sizeof_hdr != 348);
	if (m_bReverse) ReverseDsr();
	if(m_Header.hk.sizeof_hdr != 348)
	{
		ZWarning("AVWReader", "%s -- Header Size is not 348. Not a valid AVW file!", m_filename.c_str());
		return false;
	}

	m_nWidth = m_Header.dime.dim[1];
	m_nHeight = m_Header.dime.dim[2];
	m_nDepth = m_Header.dime.dim[3];
	m_nNum = m_Header.dime.dim[4];

	m_nVolsize = m_nWidth * m_nHeight * m_nDepth;

	switch(m_Header.dime.datatype)
	{
	default:
	case DT_NONE: //case DT_UNKNOWN
	  m_nVoxelbytes = 0;
	  break;
	case DT_UNSIGNED_CHAR:
	  m_nVoxelbytes = sizeof(unsigned char);
	  break;
	case DT_SIGNED_SHORT:
	  m_nVoxelbytes = sizeof(short);
	  break;
	case DT_SIGNED_INT:
	  m_nVoxelbytes = sizeof(int);
	  break;
	case DT_FLOAT:
	  m_nVoxelbytes = sizeof(float);
	  break;
	case DT_DOUBLE:
		m_nVoxelbytes = sizeof(double);
		break;
	case DT_RGB:
		m_nVoxelbytes = 3;
		break;
	}
	return true;
}

void ZAVWImage::ReverseDsr()
{
	char *ptr = (char *) &(m_Header.hk);
	Reverse4(ptr);  // sizeof_hdr
	ptr += 32;
	Reverse4(ptr);  // extents
	ptr += 4;
	Reverse2(ptr);  // session_error

	ptr = (char *) &(m_Header.dime);
	Reverse2(ptr,8);  // dims
	ptr += 28;
	Reverse2(ptr,4);  // unused1, datatype, bitpix, dim_un0
	ptr += 8;
	Reverse4(ptr,18);  // pixdim, vox_offset, ... cal_min, compressed, ... glmin

	ptr = (char *) &(m_Header.hist);
	ptr += 168;
	Reverse4(ptr,8);  // views, ... smin
}


void ZAVWImage::ReverseBuffer(void *buffer, int no_elements)
{
	switch (m_nVoxelbytes) 
	{
	case 8:
		Reverse8(buffer,no_elements);
		break;
	case 4:
		Reverse4(buffer,no_elements);
		break;
	case 2:
		Reverse2(buffer,no_elements);
		break;
	default:
		break;
	}
}

ZImageBase* ZAVWImage::ReadFile()
{
	if(!IsOpen())
	{
		ZWarning("AVWReader", "%s -- File is not opened!", m_filename.c_str());
		return 0;
	}

	if(!m_bHeader)
	{
		if(ReadHeader() == false) return 0;
		m_bHeader = true;
	}
	
	ZImageBase* image=0;
	switch(m_Header.dime.datatype)
	{
	default:
	case DT_NONE:
		ZWarning("AVWReader", "%s -- Unknown Pixel Format (%d)!", m_filename.c_str(), m_Header.dime.datatype);
		return 0;
	case DT_UNSIGNED_CHAR:
		image = new ZGrayByteImage(m_nWidth, m_nHeight, m_nDepth);
		break;
	case DT_SIGNED_SHORT:
		image = new ZGrayShortImage(m_nWidth, m_nHeight, m_nDepth);
		break;
	case DT_SIGNED_INT:
		image = new ZGrayIntImage(m_nWidth, m_nHeight, m_nDepth);
		break;
	case DT_FLOAT:
		image = new ZGrayFloatImage(m_nWidth, m_nHeight, m_nDepth);
		break;
	case DT_COMPLEX:
		ZError("AVWReader", "%s -- Complex pixel format is not supported!", m_filename.c_str());
		return 0;
	case DT_DOUBLE:
		image = new ZGrayDoubleImage(m_nWidth, m_nHeight, m_nDepth);
		break;
	case DT_RGB:
		image = new ZRGBByteImage(m_nWidth, m_nHeight, m_nDepth);
		break;
	}

	SeekVolume(m_nVolume);
	if(ReadVolumes(image->GetBuffer()) == false) 
	{
		delete image;
		return 0;
	}

	image->SetPixelDim(Abs(m_Header.dime.pixdim[1]),
					  Abs(m_Header.dime.pixdim[2]),
					  Abs(m_Header.dime.pixdim[3]));

#ifdef _WINDOWS
	if(m_Header.hist.aux_file[0] != 0 && LUTexsited(m_Header.hist.aux_file)) 
	{
		float min=m_Header.dime.cal_min, max=m_Header.dime.cal_max;
		if(!strcmp(m_Header.hist.aux_file, "render3"))
		{
			image->MinMax(min, max);
			float r = Max(Abs(min), Abs(max));
			min = -r; max = r;
		}

		image->SetDisplayRange(min, max);

		ZImageBase* pimage = GrayToRGB(image, m_Header.hist.aux_file);
		delete image; image = pimage;
	}
#endif
	return image;
}

bool ZAVWImage::WriteFile(const ZImageBase& image, bool intensityscale)
{
	if(!IsOpen())
	{
		ZWarning("AVWWriter", "%s -- File has not been opend to write!", m_filename.c_str());
		return false;
	}

	if(!image.IfValid() || image.GetBuffer() == 0)
	{
		ZWarning("AVWWriter", "%s -- Image buffer not available!", m_filename.c_str());
		return false;
	}

	InitHeader();	
	
	switch(image.PixFmt())
	{
	default:
		m_Header.dime.datatype = DT_NONE;
		m_nVoxelbytes = 0;
		ZWarning("AVWReader", "%s -- AVW doesn't support %s pixel type!", m_filename.c_str(), PixelType(image));
		return false;
	case epixfmtGrayByte:
		m_Header.dime.datatype = DT_UNSIGNED_CHAR;
	    m_nVoxelbytes = sizeof(unsigned char);
		break;
	case epixfmtGrayShort:
		m_Header.dime.datatype = DT_SIGNED_SHORT;
	    m_nVoxelbytes = sizeof(short);
		break;
	case epixfmtGrayInt:
		m_Header.dime.datatype = DT_SIGNED_INT;
	    m_nVoxelbytes = sizeof(int);
		break;
	case epixfmtGrayFloat:
		m_Header.dime.datatype = DT_FLOAT;
	    m_nVoxelbytes = sizeof(float);
		break;
	case epixfmtGrayDouble:
		m_Header.dime.datatype = DT_DOUBLE;
	    m_nVoxelbytes = sizeof(double);
		break;
	case epixfmtRGBByte:
		m_Header.dime.datatype = DT_RGB;
	    m_nVoxelbytes = 3;
		break;
	}

	m_Header.dime.bitpix = short(m_nVoxelbytes * 8);

	SetDim(image.Width(), image.Height(), image.Depth());

	m_Header.dime.pixdim[0] = 0;
	image.GetPixelDim(m_Header.dime.pixdim[1], 
					  m_Header.dime.pixdim[2], 
					  m_Header.dime.pixdim[3]);

	float min = 0, max = 256; image.MinMax(min, max); 

	if(max < 1) max = 1;
	SetMinMax(int(min), int(max+0.5));

	if(WriteVolumes(image, intensityscale) == false) return false;

	m_HeaderStream.write((char*)(&m_Header), sizeof(struct dsr));

	return true;
}

void ZAVWImage::InitHeader()
{
	short *sbuf;

	m_Header.hk.sizeof_hdr=348; /* number of bytes */
	strcpy(m_Header.hk.data_type,"");
	m_Header.hk.extents=16384;
	m_Header.hk.session_error=0;
	m_Header.hk.regular='r';
	m_Header.hk.hkey_un0=' ';
	SetDim(0, 0, 0);
	SetVoxUnits("mm");
	strcpy(m_Header.dime.cal_units,"");
	m_Header.dime.unused1=0;
	m_Header.dime.dim_un0=0;
	m_Header.dime.vox_offset=0.00;
	m_Header.dime.funused1=1.00;
	m_Header.dime.funused2=0.00;
	m_Header.dime.funused3=0.00;
	m_Header.dime.cal_max=0.00;
	m_Header.dime.cal_min=0.00;
	m_Header.dime.compressed=0;
	m_Header.dime.verified=0;
	strcpy(m_Header.hist.descrip,"Header Created by Yongyue Zhang");
	strcpy(m_Header.hist.aux_file,"");
	SetOrientation(0);

	sbuf = (short *) calloc(5, sizeof(short));
	SetOriginator(sbuf);
	free(sbuf);

	strcpy(m_Header.hist.generated,"");
	strcpy(m_Header.hist.scannum,"");
	strcpy(m_Header.hist.patient_id,"");
	strcpy(m_Header.hist.exp_date,"");
	strcpy(m_Header.hist.exp_time,"");
	strcpy(m_Header.hist.hist_un0,"");
	m_Header.hist.views=0;
	m_Header.hist.vols_added=0;
	m_Header.hist.start_field=0;
	m_Header.hist.field_skip=0;
	m_Header.hist.omax=0;
	m_Header.hist.omin=0;
	m_Header.hist.smax=0;
	m_Header.hist.smin=0;  
}

void ZAVWImage::SeekVolume(UINT vols)
{
	UINT volbytes = m_nVolsize * m_nVoxelbytes;
	long pos = volbytes * vols;

	m_ImageStream.seekg(pos, std::ios::beg);
}

bool ZAVWImage::ReadVolumes(PBYTE buffer)
{
	UINT volbytes = m_nVolsize * m_nVoxelbytes;

	m_ImageStream.read((char*)buffer, volbytes);
	if(m_ImageStream.fail())
	{
		ZWarning("AVWReader", "%s -- Error occurs when reading!", m_filename.c_str());
		return false;
	}

	if(m_bReverse)		//need reverse
		ReverseBuffer(buffer, m_nVolsize);

	return true;
}

bool ZAVWImage::WriteVolumes(const ZImageBase& image, bool intensityscale)
{
	ZImageBase* pImage = const_cast<ZImageBase*>(&image);
	if(intensityscale)
	{
		pImage = TypeCopyFrom(&image);
		image.CopyPixelsTo(*pImage, 0, image.RangeMax());
	}
	if(WriteImage(*pImage, m_ImageStream) == false)
	{
		ZWarning("AVWWriter", "%s -- Error occurs when writing!", m_filename.c_str());
		return false;
	}
	if(intensityscale) delete pImage;

	return true;
}

void ZAVWImage::SetDim(int x, int y, int z)
{
	m_Header.dime.dim[0] = 4;
	m_Header.dime.dim[1] = short(x);
	m_Header.dime.dim[2] = short(y);
	m_Header.dime.dim[3] = short(z);
	m_Header.dime.dim[4] = 1;
	m_Header.dime.dim[5] = 0;
	m_Header.dime.dim[6] = 0;
	m_Header.dime.dim[7] = 0;
}

void ZAVWImage::SetVoxUnits(const char *units)
{
	strcpy(m_Header.dime.vox_units, units);
}

void ZAVWImage::SetOrientation(const char orient)
{
	m_Header.hist.orient = orient;
}


void ZAVWImage::SetOriginator(const short orig[5])
{
	memcpy((void *) m_Header.hist.originator, (void *) orig, 5*sizeof(short));
}


void ZAVWImage::SetMinMax(int min, int max)
{
	m_Header.dime.glmin = min;
	m_Header.dime.glmax = max;
}
