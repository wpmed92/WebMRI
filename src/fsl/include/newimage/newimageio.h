/*  General IO functions (images and transformation files)

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 3.3 (c) 2006, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/1112. */


#if !defined(__newimageio_h)
#define __newimageio_h

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include "newmatio.h"
#include "newimage.h"
#include "complexvolume.h"
#include "fslio/fslio.h"
#include "miscmaths/miscmaths.h"

using namespace NEWMAT;

namespace NEWIMAGE {

//class volumeinfo {
//  private:
//    FSLIO p_fslio;
//  public:
//    volumeinfo() { p_fslio = blank_vinfo(); }
//    FSLIO operator() () { return p_fslio; }
//}

typedef FSLIO volumeinfo;

volumeinfo blank_vinfo();

string fslbasename(const string& filename);
int make_basename(string& filename);
int find_pathname(string& filename);
bool fsl_imageexists(const string& filename);

  // read

template <class T>
int read_volume(volume<T>& target, const string& filename);
template <class T>
int read_volume(volume<T>& target, const string& filename, volumeinfo& vinfo);
template <class T>
int read_volume_hdr_only(volume<T>& target, const string& filename);
template <class T>
int read_volume_hdr_only(volume<T>& target, const string& filename, 
			 volumeinfo& vinfo);
int read_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename);
int read_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename, volumeinfo& vinfo);
int read_complexvolume(complexvolume& vol, const string& filename);
int read_complexvolume(complexvolume& vol, const string& filename, 
		       volumeinfo& vinfo);

template <class T>
int read_volume4D(volume4D<T>& target, const string& filename);
template <class T>
int read_volume4D(volume4D<T>& target, const string& filename, 
		  volumeinfo& vinfo);
template <class T>
int read_volume4D_hdr_only(volume4D<T>& target, const string& filename);
template <class T>
int read_volume4D_hdr_only(volume4D<T>& target, const string& filename,
			   volumeinfo& vinfo);
int read_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename);
int read_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename, volumeinfo& vinfo);
int read_complexvolume4D(complexvolume& vol, const string& filename);
int read_complexvolume4D(complexvolume &vol, const string& filename, 
			 volumeinfo& vinfo);


  // save

template <class T>
int save_volume(const volume<T>& source, const string& filename);
template <class T>
int save_volume(const volume<T>& source, const string& filename,
		const volumeinfo& vinfo);
int save_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename);
int save_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename,
		       const volumeinfo& vinfo);
int save_complexvolume(const complexvolume& vol, const string& filename);
int save_complexvolume(const complexvolume& vol, const string& filename,
		       const volumeinfo& vinfo);

template <class T>
int save_volume4D(const volume4D<T>& source, const string& filename);
template <class T>
int save_volume4D(const volume4D<T>& source, const string& filename,
		  const volumeinfo& vinfo);
int save_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename);
int save_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename,
			 const volumeinfo& vinfo);
int save_complexvolume4D(const complexvolume& vol, const string& filename);
int save_complexvolume4D(const complexvolume& vol, const string& filename,
			 const volumeinfo& vinfo);

template <class T>
int save_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype, const volumeinfo& vinfo);
template <class T>
int save_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype);
template <class T>
int save_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype, const volumeinfo& vinfo);
template <class T>
int save_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype);

template <class T>
int save_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype, const volumeinfo& vinfo);
template <class T>
int save_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype);
template <class T>
int save_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype, const volumeinfo& vinfo);
template <class T>
int save_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype);


  // Transform matrix IO
  
int read_ascii_matrix(Matrix &target, const string& filename);

  // outdated MEDx-compliant io routines
template <class T>
int read_matrix(Matrix &target, const string& filename, const volume<T>& invol);
template <class T>
int read_matrix(Matrix &target, const string& filename, const volume<T>& invol,
		const volume<T>& finalvol);
template <class T>
int write_medx_matrix(const Matrix& worldmat, const string& filename, 
		      const volume<T>& initvol, 
		      const volume<T>& finalvol, 
		      const string& mtype, const string& reffname);

// Helper functions

short dtype(const char* T);
short dtype(const short* T);
short dtype(const int* T);
short dtype(const float* T);
short dtype(const double* T);

short dtype(const volume<char>& vol);
short dtype(const volume<short>& vol);
short dtype(const volume<int>& vol);
short dtype(const volume<float>& vol);
short dtype(const volume<double>& vol);

short dtype(const volume4D<char>& vol);
short dtype(const volume4D<short>& vol);
short dtype(const volume4D<int>& vol);
short dtype(const volume4D<float>& vol);
short dtype(const volume4D<double>& vol);

short dtype(const string& filename);

volumeinfo volinfo(const string& filename);


// Boring overloads to enable different names (load and write)


// load

template <class T>
int load_volume(volume<T>& target, const string& filename);
template <class T>
int load_volume(volume<T>& target, const string& filename, volumeinfo& vinfo);
template <class T>
int load_volume_hdr_only(volume<T>& target, const string& filename);
template <class T>
int load_volume_hdr_only(volume<T>& target, const string& filename, 
			 volumeinfo& vinfo);
int load_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename);
int load_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename, volumeinfo& vinfo);
int load_complexvolume(complexvolume& vol, const string& filename);
int load_complexvolume(complexvolume& vol, const string& filename, 
		       volumeinfo& vinfo);

template <class T>
int load_volume4D(volume4D<T>& target, const string& filename);
template <class T>
int load_volume4D(volume4D<T>& target, const string& filename, 
		  volumeinfo& vinfo);
template <class T>
int load_volume4D_hdr_only(volume4D<T>& target, const string& filename);
template <class T>
int load_volume4D_hdr_only(volume4D<T>& target, const string& filename,
			   volumeinfo& vinfo);
int load_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename);
int load_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename, volumeinfo& vinfo);
// int load_complexvolume4D(complexvolume& vol, const string& filename);
// int load_complexvolume4D(complexvolume &vol, const string& filename, 
// 			 volumeinfo& vinfo);


  // write

template <class T>
int write_volume(const volume<T>& source, const string& filename);
template <class T>
int write_volume(const volume<T>& source, const string& filename,
		const volumeinfo& vinfo);
int write_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename);
int write_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename,
		       const volumeinfo& vinfo);
int write_complexvolume(const complexvolume& vol, const string& filename);
int write_complexvolume(const complexvolume& vol, const string& filename,
		       const volumeinfo& vinfo);

template <class T>
int write_volume4D(const volume4D<T>& source, const string& filename);
template <class T>
int write_volume4D(const volume4D<T>& source, const string& filename,
		  const volumeinfo& vinfo);
int write_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename);
int write_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename,
			 const volumeinfo& vinfo);
// int write_complexvolume4D(const complexvolume& vol, const string& filename);
// int write_complexvolume4D(const complexvolume& vol, const string& filename,
// 			 const volumeinfo& vinfo);

template <class T>
int write_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype, const volumeinfo& vinfo);
template <class T>
int write_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype);
template <class T>
int write_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype, const volumeinfo& vinfo);
template <class T>
int write_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype);

template <class T>
int write_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype, const volumeinfo& vinfo);
template <class T>
int write_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype);
template <class T>
int write_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype, const volumeinfo& vinfo);
template <class T>
int write_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype);




////////////////////////////////////////////////////////////////////////
///////////////////////// TEMPLATE DEFINITIONS /////////////////////////
////////////////////////////////////////////////////////////////////////

template <class T>
void FslReadBuffer(FSLIO* IP, T* tbuffer) 
{
  short sx,sy,sz,st;
  FslGetDim(IP,&sx,&sy,&sz,&st);
  size_t imagesize=sx*sy*sz;
  short type;
  float slope, intercept;
  bool doscaling = false;

  FslGetDataType(IP,&type);
  doscaling = FslGetIntensityScaling(IP,&slope,&intercept);
  if ( (dtype(tbuffer) == type) && (!doscaling) ) {
    FslReadVolumes(IP,tbuffer,1);
  } else {
    switch(type)
      {
      case DT_SIGNED_SHORT:
	{
	  short* sbuffer=new short[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_UNSIGNED_CHAR:
	{
	  unsigned char* sbuffer=new unsigned char[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_SIGNED_INT:
	{
	  int* sbuffer=new int[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_FLOAT:
	{
	  float* sbuffer=new float[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_DOUBLE:
	{
	  double* sbuffer=new double[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
	/*------------------- new codes for NIFTI ---*/
      case DT_INT8:
	{
	  signed char* sbuffer=new signed char[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_UINT16:
	{
	  unsigned short* sbuffer=new unsigned short[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_UINT32:
	{
	  unsigned int* sbuffer=new unsigned int[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_INT64:
	{
	  long signed int* sbuffer=new long signed int[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      case DT_UINT64:
	{
	  long unsigned int* sbuffer=new long unsigned int[imagesize];
	  if (sbuffer==0) { imthrow("Out of memory",99); }
	  FslReadVolumes(IP,sbuffer,1);
	  if (doscaling) convertbuffer(sbuffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer(sbuffer,tbuffer,imagesize);
	  delete[] sbuffer;
	}
	break;
      default:
	  /* includes: DT_BINARY, DT_RGB, DT_ALL, DT_FLOAT128, DT_COMPLEX's */
	ostringstream errmsg;
	errmsg << "Fslread: DT " << type <<  " not supported";
	perror(errmsg.str().c_str());
      }
  }
}


void check_filename(const string& basename);

FSLIO* NewFslOpen(const string& filename, const string& permissions, 
		  int filetype, const volumeinfo& vinfo, bool use_vinfo);

FSLIO* NewFslOpen(const string& filename, const string& permissions, 
		  const volumeinfo& vinfo, bool use_vinfo);

FSLIO* NewFslOpen(const string& filename, const string& permissions);

// External functions

// READ FUNCTIONS


template <class T>
void set_volume_properties(FSLIO* IP1, volume<T>& target)
{
  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);
  target.setdims(x,y,z);

  int sform_code, qform_code;
  mat44 smat, qmat;
  sform_code = FslGetStdXform(IP1,&smat);
  qform_code = FslGetRigidXform(IP1,&qmat);
  Matrix snewmat(4,4), qnewmat(4,4);
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      snewmat(i,j) = smat.m[i-1][j-1];
      qnewmat(i,j) = qmat.m[i-1][j-1];
    }
  }
  target.set_sform(sform_code,snewmat);
  target.set_qform(qform_code,qnewmat);

  short intent_code;
  float p1, p2, p3;
  FslGetIntent(IP1, &intent_code, &p1, &p2, &p3);
  target.set_intent(intent_code,p1,p2,p3);
}



template <class T>
int read_volume(volume<T>& target, const string& filename, 
		volumeinfo& vinfo, short& dtype, bool read_img_data)
{
  Tracer trcr("read_volume");

  FSLIO *IP1;
  IP1 = NewFslOpen(filename.c_str(), "r");
  if (IP1==0) { imthrow("Failed to read volume "+filename,22); }
  short sx,sy,sz,st;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  size_t volsize=sx*sy*sz;

  T* tbuffer;
  if (read_img_data) {
    tbuffer = new T[volsize];
    if (tbuffer==0) { imthrow("Out of memory",99); }
    FslReadBuffer(IP1,tbuffer);
  } else {
    tbuffer  = new T[volsize];  // a hack to stop reinitialize from allocating memory //originally 1
  }
  target.reinitialize(sx,sy,sz,tbuffer,true);

  FslGetDataType(IP1,&dtype);
  set_volume_properties(IP1,target);

  vinfo = blank_vinfo();
  FslCloneHeader(&vinfo,IP1);
  FslSetFileType(&vinfo,FslGetFileType(IP1));
  FslClose(IP1);
  return 0;
}


template <class T>
int read_volume(volume<T>& target, const string& filename, volumeinfo& vinfo)
{
  short dtype;
  int retval = read_volume(target,filename,vinfo,dtype,true);
  return retval;
}

template <class T>
int read_volume(volume<T>& target, const string& filename)
{
  short dtype;
  volumeinfo vinfo = blank_vinfo();
  int retval = read_volume(target,filename,vinfo,dtype,true);
  return retval;
}

template <class T>
int read_volume_hdr_only(volume<T>& target, const string& filename, 
			 volumeinfo& vinfo)
{
  short dtype;
  int retval = read_volume(target,filename,vinfo,dtype,false);
  return retval;
}

template <class T>
int read_volume_hdr_only(volume<T>& target, const string& filename)
{
  short dtype;
  volumeinfo vinfo = blank_vinfo();
  int retval = read_volume(target,filename,vinfo,dtype,false);
  return retval;
}


template <class T>
int read_volume4D(volume4D<T>& target, const string& filename, 
		  volumeinfo& vinfo, short& dtype, bool read_img_data)
{
  Tracer trcr("read_volume4D");

  target.destroy();

  FSLIO *IP1;
  IP1 = NewFslOpen(filename.c_str(), "r");
  if (IP1==0) { imthrow("Failed to read volume "+filename,22); }

  short sx,sy,sz,st;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  size_t volsize=sx*sy*sz;
  if (st<1) st=1;  // make it robust to dim4<1

  volume<T> dummyvol(sx,sy,sz);
  for (int t=0; t<st; t++) {
    target.addvolume(dummyvol);
    T* tbuffer;
    if (read_img_data) {
      tbuffer = new T[volsize];
      if (tbuffer==0) { imthrow("Out of memory",99); }
      FslReadBuffer(IP1,tbuffer);
    } else {
      tbuffer = new T[volsize];  // set 1 as a bad hack to stop reinitialize from allocating memory // 
    }
    // Note that the d_owner flag = true in the following so that the
    //  control for delete is passed to the volume class
    target[t].reinitialize(sx,sy,sz,tbuffer,true);
    set_volume_properties(IP1,target[t]);
  }

  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);
  target.setdims(x,y,z,tr);

  FslGetDataType(IP1,&dtype);

  vinfo = blank_vinfo();
  FslCloneHeader(&vinfo,IP1);
  FslSetFileType(&vinfo,FslGetFileType(IP1));
  FslClose(IP1);
  return 0;
}


template <class T>
int read_volume4D(volume4D<T>& target, const string& filename, 
		  volumeinfo& vinfo)
{
  short dtype;
  int retval = read_volume4D(target,filename,vinfo,dtype,true);
  return retval;
}

template <class T>
int read_volume4D(volume4D<T>& target, const string& filename)
{
  short dtype;
  volumeinfo vinfo = blank_vinfo();
  int retval = read_volume4D(target,filename,vinfo,dtype,true);
  return retval;
}

template <class T>
int read_volume4D_hdr_only(volume4D<T>& target, const string& filename,
			   volumeinfo& vinfo)
{
  short dtype;
  int retval = read_volume4D(target,filename,vinfo,dtype,false);
  return retval;
}

template <class T>
int read_volume4D_hdr_only(volume4D<T>& target, const string& filename)
{
  short dtype;
  volumeinfo vinfo = blank_vinfo();
  int retval = read_volume4D(target,filename,vinfo,dtype,false);
  return retval;
}


// SAVE FUNCTIONS


mat44 newmat2mat44(const Matrix& nmat);


template <class T>
int set_fsl_hdr(const volume<T>& source, FSLIO *OP, int tsize, float tdim)
{
  Tracer tr("set_fsl_hdr");
    
  FslSetDim(OP,source.xsize(),source.ysize(),source.zsize(),tsize);
  FslSetDataType(OP, dtype(source));
  FslSetVoxDim(OP,source.xdim(), source.ydim(), source.zdim(), tdim);

  FslSetStdXform(OP,source.sform_code(),newmat2mat44(source.sform_mat()));
  FslSetRigidXform(OP,source.qform_code(),newmat2mat44(source.qform_mat()));
  
  FslSetIntent(OP,source.intent_code(),source.intent_param(1),
	       source.intent_param(2),source.intent_param(3));
  
  return 0;
}



template <class T>
int save_basic_volume(const volume<T>& source, const string& filename, 
		      int filetype, const volumeinfo& vinfo, bool use_vinfo)
{
  // if filetype < 0 then it is ignored, otherwise it overrides everything
  Tracer tr("save_basic_volume");

  FSLIO *OP = NewFslOpen(filename.c_str(),"wb",filetype,vinfo,use_vinfo);
  if (OP==0) { imthrow("Failed to open volume "+filename+" for writing",23); }
  set_fsl_hdr(source,OP,1,1);
  FslWriteAllVolumes(OP,&(source(0,0,0)));
  FslClose(OP);
  return 0;
}


template <class T>
int save_basic_volume4D(const volume4D<T>& source, const string& filename,
			int filetype, const volumeinfo& vinfo, bool use_vinfo)
{
  Tracer tr("save_basic_volume4D");
  if (source.tsize()<1) return -1;

  // if filetype < 0 then it is ignored, otherwise it overrides everything
  FSLIO *OP = NewFslOpen(filename.c_str(),"wb",filetype,vinfo,use_vinfo);
  if (OP==0) { imthrow("Failed to open volume "+filename+" for writing",23); }

  set_fsl_hdr(source[0],OP,source.tsize(),source.tdim());
  if (filetype>=0) FslSetFileType(OP,filetype);
  FslWriteHeader(OP);
  if (source.nvoxels()>0) {
    for (int t=0; t<source.tsize(); t++) {
      FslWriteVolumes(OP,&(source[t](0,0,0)),1);
    }
  }
  FslClose(OP); 

  return 0;
}


template <class T>
int save_volume(const volume<T>& source, const string& filename, 
		      const volumeinfo& vinfo, bool use_vinfo)
{
  return save_basic_volume(source,fslbasename(filename),-1,vinfo,use_vinfo);
}


template <class T>
int save_volume(const volume<T>& source, const string& filename)
{
  volumeinfo vinfo = blank_vinfo();
  return save_volume(source,filename,vinfo,false);
}


template <class T>
int save_volume(const volume<T>& source, const string& filename,
		const volumeinfo& vinfo)
{
  return save_volume(source,filename,vinfo,true);
}


template <class T>
int save_volume4D(const volume4D<T>& source, const string& filename,
		  const volumeinfo& vinfo, bool use_vinfo)
{
  return save_basic_volume4D(source,fslbasename(filename),-1,vinfo,use_vinfo);
}

template <class T>
int save_volume4D(const volume4D<T>& source, const string& filename)
{
  volumeinfo vinfo = blank_vinfo();
  return save_volume4D(source,filename,vinfo,false);
}


template <class T>
int save_volume4D(const volume4D<T>& source, const string& filename,
		  const volumeinfo& vinfo)
{
  return save_volume4D(source,filename,vinfo,true);
}


template <class T>
int save_volume_dtype(const volume<T>& source, const string& filename,
		       short datatype, const volumeinfo& vinfo, bool use_vinfo)
{
  if (dtype(source) == datatype) {
    return save_volume(source,filename,vinfo,use_vinfo);
  } else {
    switch(datatype)
      {
      case DT_SIGNED_SHORT:
	{
	  volume<short> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename,vinfo,use_vinfo);
	}
	break;
      case DT_UNSIGNED_CHAR:
	{
	  volume<char> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename,vinfo,use_vinfo);
	}
	break;
      case DT_SIGNED_INT:
	{
	  volume<int> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename,vinfo,use_vinfo);
	}
	break;
      case DT_FLOAT:
	{
	  volume<float> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename,vinfo,use_vinfo);
	}
	break;
      case DT_DOUBLE:
	{
	  volume<double> svol;
	  copyconvert(source,svol);
	  return save_volume(svol,filename,vinfo,use_vinfo);
	}
	break;
      default:
	ostringstream errmsg;
	errmsg << "Fslread: DT " << datatype <<  " not supported";
	perror(errmsg.str().c_str());
      }
  }
  return -1;  // should never get here
}
  

template <class T>
int save_volume4D_dtype(const volume4D<T>& source, const string& filename,
		       short datatype, const volumeinfo& vinfo, bool use_vinfo)
{
  if (dtype(source) == datatype) {
    return save_volume4D(source,filename,vinfo,use_vinfo);
  } else {
    switch(datatype)
      {
      case DT_SIGNED_SHORT:
	{
	  volume4D<short> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename,vinfo,use_vinfo);
	}
	break;
      case DT_UNSIGNED_CHAR:
	{
	  volume4D<char> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename,vinfo,use_vinfo);
	}
	break;
      case DT_SIGNED_INT:
	{
	  volume4D<int> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename,vinfo,use_vinfo);
	}
	break;
      case DT_FLOAT:
	{
	  volume4D<float> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename,vinfo,use_vinfo);
	}
	break;
      case DT_DOUBLE:
	{
	  volume4D<double> svol;
	  copyconvert(source,svol);
	  return save_volume4D(svol,filename,vinfo,use_vinfo);
	}
	break;
      default:
	ostringstream errmsg;
	errmsg << "Fslread: DT " << datatype <<  " not supported";
	perror(errmsg.str().c_str());
      }
  }
  return -1;  // should never get here
}

template <class T>
int save_volume_datatype(const volume<T>& source, const string& filename,
		      short datatype, const volumeinfo& vinfo)
{  
  return save_volume_dtype(source,filename,datatype,vinfo,true);
}

template <class T>
int save_volume_datatype(const volume<T>& source, const string& filename,
		      short datatype)
{  
  volumeinfo vinfo = blank_vinfo();
  return save_volume_dtype(source,filename,datatype,vinfo,false);
}


template <class T>
int save_volume4D_datatype(const volume4D<T>& source, const string& filename,
			 short datatype, const volumeinfo& vinfo)
{  
  return save_volume4D_dtype(source,filename,datatype,vinfo,true);
}

template <class T>
int save_volume4D_datatype(const volume4D<T>& source, const string& filename,
			 short datatype)
{  
  volumeinfo vinfo = blank_vinfo();
  return save_volume4D_dtype(source,filename,datatype,vinfo,false);
}


// old versions call _dtype - kept for compatability
  
template <class T>
int save_volume_dtype(const volume<T>& source, const string& filename,
			 short datatype, const volumeinfo& vinfo)
{
  return save_volume_datatype(source,filename,datatype,vinfo);
}

template <class T>
int save_volume_dtype(const volume<T>& source, const string& filename,
			 short datatype)
{
  return save_volume_datatype(source,filename,datatype);
}

template <class T>
int save_volume4D_dtype(const volume4D<T>& source, const string& filename,
			   short datatype, const volumeinfo& vinfo)
{
  return save_volume4D_datatype(source,filename,datatype,vinfo);
}

template <class T>
int save_volume4D_dtype(const volume4D<T>& source, const string& filename,
			   short datatype)
{
  return save_volume4D_datatype(source,filename,datatype);
}



template <class T>
int save_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype, const volumeinfo& vinfo)
{
  return save_basic_volume(source,filename,filetype,vinfo,true);
}


template <class T>
int save_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype)
{
  volumeinfo vinfo = blank_vinfo();
  return save_basic_volume(source,filename,filetype,vinfo,false); 
}


template <class T>
int save_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype, const volumeinfo& vinfo)
{
  return save_basic_volume4D(source,filename,filetype,vinfo,true);
}


template <class T>
int save_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype)
{
  volumeinfo vinfo = blank_vinfo();
  return save_basic_volume4D(source,filename,filetype,vinfo,false); 
}



// MATRIX IO

// Non-templated helper functions (bodies in .cc)

int get_medx_small_matrix(Matrix &target, ifstream& matfile);
int get_medx_matrix(Matrix &target, ifstream& matfile);
int get_minc_matrix(Matrix &target, ifstream& matfile);
int put_medx_matrix(ofstream& matfile, const string& name, 
		    const Matrix& affmat);

// Template definitions

template <class T>
int medx2world(Matrix &target, const Matrix& medxmat, 
	       const ColumnVector& finalimsize, const Matrix& finalvoxsize,
	       const volume<T>& initvol) 
{
  Tracer tr("medx2world");
  // convert the voxel->voxel MEDx form to a world->world form
  Matrix swapy1(4,4), samp1(4,4);
  Matrix swapy2(4,4), samp2(4,4);
  Identity(swapy1);
  swapy1(2,2)=-1;
  swapy2=swapy1;
  swapy1(2,4)=initvol.ysize()-1.0;  // corrected version
  swapy2(2,4)=finalimsize(2)-1.0;   // corrected version
  Identity(samp1);
  samp1(1,1) = initvol.xdim();  
  samp1(2,2) = initvol.ydim();
  samp1(3,3) = initvol.zdim();
  samp2 = finalvoxsize;
  samp2(1,4)=0; samp2(2,4)=0; samp2(3,4)=0;  // translation can be ignored
  samp2(1,1)=fabs(samp2(1,1));    // swapping info is not used
  samp2(2,2)=fabs(samp2(2,2)); 
  samp2(3,3)=fabs(samp2(3,3));
  target = samp2 * swapy2 * medxmat * swapy1 * samp1.i(); // correct
  return 0;
}


template <class T>
int world2medx(Matrix &medxmat, const Matrix& worldmat, 
	       const ColumnVector& finalimsize, const Matrix& finalvoxsize,
	       const volume<T>& initvol) 
{
  Tracer tr("world2medx");
  // convert the world->world form to a voxel->voxel MEDx form 
  Matrix swapy1(4,4), samp1(4,4);
  Matrix swapy2(4,4), samp2(4,4);
  Identity(swapy1);
  swapy1(2,2)=-1;
  swapy2=swapy1;
  swapy1(2,4)=initvol.ysize()-1.0;  // corrected version
  swapy2(2,4)=finalimsize(2)-1.0;   // corrected version
  Identity(samp1);
  samp1(1,1) = initvol.xdim();  
  samp1(2,2) = initvol.ydim();
  samp1(3,3) = initvol.zdim();
  samp2 = finalvoxsize;
  samp2(1,4)=0; samp2(2,4)=0; samp2(3,4)=0;  // translation can be ignored
  samp2(1,1)=fabs(samp2(1,1));    // swapping info is not used
  samp2(2,2)=fabs(samp2(2,2)); 
  samp2(3,3)=fabs(samp2(3,3));
  medxmat = swapy2 * samp2.i() * worldmat * samp1 * swapy1;
  return 0;
}


template <class T>
int minc2world(Matrix &target, const Matrix& mincmat, 
	       const volume<T>& initvol, const volume<T>& finalvol)
{
  Tracer tr("minc2world");
  // convert the voxel->voxel MINC form to a world->world form
  Matrix trans1(4,4), trans2(4,4), samp1(4,4), samp2(4,4), 
         flipx1(4,4), flipx2(4,4);

  Identity(trans1);
  if (initvol.sform_code()!=NIFTI_XFORM_UNKNOWN) {
    trans1.SubMatrix(1,3,4,4) = initvol.sform_mat().i().SubMatrix(1,3,4,4);
  }

  Identity(trans2);
  if (finalvol.sform_code()!=NIFTI_XFORM_UNKNOWN) {
    trans2.SubMatrix(1,3,4,4) = finalvol.sform_mat().i().SubMatrix(1,3,4,4);
  }
  
  Identity(samp1);
  samp1(1,1) = fabs(initvol.xdim());  
  samp1(2,2) = fabs(initvol.ydim());
  samp1(3,3) = fabs(initvol.zdim());

  Identity(samp2);
  samp2(1,1) = fabs(finalvol.xdim());  
  samp2(2,2) = fabs(finalvol.ydim());
  samp2(3,3) = fabs(finalvol.zdim());
  
  Identity(flipx1);
  flipx1(1,1) = -1.0;
  flipx1(1,4) = initvol.xsize() - 1.0;

  Identity(flipx2);
  flipx2(1,1) = -1.0;
  flipx2(1,4) = finalvol.xsize() - 1.0;

  target = samp2 * flipx2 * trans2 * samp2.i() * mincmat * samp1 *
                  trans1.i() * flipx1 * samp1.i();
  return 0;
}


template <class T>
int read_medx_matrix(Matrix &target, ifstream& matfile, const volume<T>& invol)
{
  Tracer tr("read_medx_matrix");
  string curword;
  Matrix id4(4,4), gmat, imsize(1,3), voxsize(4,4), medxmat(4,4);
  Identity(id4);  // set to be the identity matrix for checking
  medxmat = id4;   // default return value
  imsize = 1.0;
  voxsize = id4;
  while (!matfile.eof()) {
    matfile >> curword;
    if ( (curword == "/outputsize") ) {
      get_medx_small_matrix(imsize,matfile);
    }
    if ( (curword == "/outputusermatrix") ) {
      get_medx_small_matrix(voxsize,matfile);
    }
    if ( (curword == "/GenericReslice")  || (curword == "/ShadowTransform")
	 || (curword == "/MotionCorrectionReslice") 
	 || (curword == "/UserTransformation") 
	 || (curword == "/IntoTalairachSpace") ) {
      if (get_medx_matrix(gmat,matfile)==0) {  // only accept good returns
	if (SumSquare(gmat - id4) > 0.001) // only accept non-identities
	  medxmat = gmat;
      }
    }
  }
  ColumnVector imsizev(3);
  imsizev(1) = imsize(1,1);
  imsizev(2) = imsize(1,2);
  imsizev(3) = imsize(1,3);
  medx2world(target,medxmat,imsizev,voxsize,invol);
  return 0;
}


template <class T>
int read_minc_matrix(Matrix &target, ifstream& matfile, const volume<T>& invol,
		     const volume<T>& finalvol)
{
  Tracer tr("read_minc_matrix");
  string curword;
  Matrix gmat, id4(4,4);
  Identity(id4);
  target = id4;  // default return value
  while (!matfile.eof()) {
    matfile >> curword;
    if ( (curword == "Linear_Transform") ) {
      get_minc_matrix(gmat,matfile);
    }
    curword="";
  }
  minc2world(target,gmat,invol,finalvol);
  return 0;
}


template <class T>
int read_matrix(Matrix &target, const string& filename, const volume<T>& invol,
		const volume<T>& finalvol)
{
  Tracer tr("read_matrix");
  if ( filename.size()<1 ) return -1;
  ifstream matfile(filename.c_str());
  if (!matfile) { 
    cerr << "Could not open matrix file " << filename << endl;
    return -1;
  }
  string firstline;
  matfile >> firstline;
  matfile.seekg(0,ios::beg);
  int rval=0;
  if (firstline == "%!VEST-Transformations")
    {
      rval = read_medx_matrix(target,matfile,invol);
    } else if (firstline == "MNI") {
      if (finalvol.xsize()>0) {
	rval = read_minc_matrix(target,matfile,invol,finalvol);
      } else {
	cerr << "Cannot read MINC transform files with finalvol" << endl;
	rval = -1;
      }
    } else {
      target = MISCMATHS::read_ascii_matrix(matfile);
      if (target.Nrows()<=0) rval = -1;
      else rval = 0;
    }
  matfile.close();
  return rval;
}

template <class T>
int read_matrix(Matrix &target, const string& filename, const volume<T>& invol)
{
  volume<T> dummy(0,0,0);
  return read_matrix(target,filename,invol,dummy);
}

int get_outputusermat(const string& filename, Matrix& oumat);

template <class T>
int write_medx_matrix(const Matrix& worldmat, const string& filename, 
		      const volume<T>& initvol, 
		      const volume<T>& finalvol, 
		      const string& mtype, const string& reffname)
{
  Tracer tr("write_medx_matrix");
  if ( (filename.size()<1) ) return -1;
  ofstream matfile(filename.c_str());
  if (!matfile) { 
    cerr << "Could not open file " << filename << " for writing" << endl;
    return -1;
  }

  // convert from world->world transform to the voxel->voxel MEDx transform
  Matrix swapy1(4,4), swapy2(4,4), samp1(4,4), samp2(4,4), medxtrans(4,4);
  Identity(samp1);  Identity(samp2);  Identity(swapy1);  Identity(swapy2);
  swapy1(2,2) = -1.0;
  swapy2(2,2) = -1.0;
  swapy1(2,4) = initvol.ysize()-1.0;  // corrected version
  swapy2(2,4) = finalvol.ysize()-1.0; // corrected version
  samp1(1,1) = initvol.xdim();  
  samp1(2,2) = initvol.ydim();
  samp1(3,3) = initvol.zdim();
  samp2(1,1) = finalvol.xdim();  
  samp2(2,2) = finalvol.ydim();
  samp2(3,3) = finalvol.zdim();
  medxtrans = swapy2 * samp2.i() * worldmat * samp1 * swapy1;


  string mattype;
  if ( (mtype.size()<1) || // default
       (mtype[0] == 'g') || (mtype[0] == 'G') )  
    { mattype="GenericReslice"; }
  if ( (mtype[0] == 's') || (mtype[0] == 'S') )  
    { mattype="ShadowTransform"; }
  if ( (mtype[0] == 't') || (mtype[0] == 'T') )  
    { mattype="IntoTalairachSpace"; }
  if ( (mtype[0] == 'u') || (mtype[0] == 'U') )  
    { mattype="UserTransformation"; }
  if ( (mtype[0] == 'm') || (mtype[0] == 'M') || 
       (mtype[0] == 'a') || (mtype[0] == 'A') )
    { mattype="MotionCorrectionReslice"; }

  matfile << "%!VEST-Transformations" << endl;
  matfile << "<<" << endl;
  matfile << "    /" << mattype << "     <<" << endl;
  put_medx_matrix(matfile,"matrix",medxtrans);
  matfile << "        /order  1" << endl;

  Matrix imsize(1,3);
  imsize(1,1) = finalvol.xsize();
  imsize(1,2) = finalvol.ysize();
  imsize(1,3) = finalvol.zsize();
  put_medx_matrix(matfile,"outputsize",imsize);

  Matrix oumat(4,4);
  get_outputusermat(reffname,oumat);
  put_medx_matrix(matfile,"outputusermatrix",oumat);

  if (mattype == "MotionCorrectionReslice") {
    matfile << "        /outputunits    (mm)" << endl;
    matfile << "        /talairachcalibrated?   false" << endl;
  }
  if (mattype == "IntoTalairachSpace") {
    matfile << "        /talairachcalibrated?   true" << endl;
  }
  matfile << "    >>" << endl;
  matfile << ">>" << endl;
  matfile.close();
  return 0;
}


////////////////////////////////////////////////////////////////////////
///// Boring overloads to enable different names (load and write) //////
////////////////////////////////////////////////////////////////////////



// load

template <class T>
int load_volume(volume<T>& target, const string& filename)
  { return read_volume(target,filename); }

template <class T>
int load_volume(volume<T>& target, const string& filename, volumeinfo& vinfo)
  { return read_volume(target,filename,vinfo); }

template <class T>
int load_volume_hdr_only(volume<T>& target, const string& filename)
  { return read_volume_hdr_only(target,filename); }

template <class T>
int load_volume_hdr_only(volume<T>& target, const string& filename, 
			 volumeinfo& vinfo)
  { return read_volume_hdr_only(target,filename,vinfo); }

template <class T>
int load_volume4D(volume4D<T>& target, const string& filename)
  { return read_volume4D(target,filename); }

template <class T>
int load_volume4D(volume4D<T>& target, const string& filename, 
		  volumeinfo& vinfo)
  { return read_volume4D(target,filename,vinfo); }

template <class T>
int load_volume4D_hdr_only(volume4D<T>& target, const string& filename)
  { return read_volume4D_hdr_only(target,filename); }

template <class T>
int load_volume4D_hdr_only(volume4D<T>& target, const string& filename,
			   volumeinfo& vinfo)
  { return read_volume4D_hdr_only(target,filename,vinfo); }

  // write

template <class T>
int write_volume(const volume<T>& source, const string& filename)
  { return save_volume(source,filename); }

template <class T>
int write_volume(const volume<T>& source, const string& filename,
		const volumeinfo& vinfo)
  { return save_volume(source,filename,vinfo); }

template <class T>
int write_volume4D(const volume4D<T>& source, const string& filename)
  { return save_volume4D(source,filename); }

template <class T>
int write_volume4D(const volume4D<T>& source, const string& filename,
		  const volumeinfo& vinfo)
  { return save_volume4D(source,filename,vinfo); }

template <class T>
int write_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype, const volumeinfo& vinfo)
  { return save_volume_datatype(source,filename,datatype,vinfo); }

template <class T>
int write_volume_datatype(const volume<T>& source, const string& filename,
			 short datatype)
  { return save_volume_datatype(source,filename,datatype); }

template <class T>
int write_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype, const volumeinfo& vinfo)
  { return save_volume4D_datatype(source,filename,datatype,vinfo); }

template <class T>
int write_volume4D_datatype(const volume4D<T>& source, const string& filename,
			   short datatype)
  { return save_volume4D_datatype(source,filename,datatype); }

template <class T>
int write_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype, const volumeinfo& vinfo)
  { return save_volume_filetype(source,filename,filetype,vinfo); }

template <class T>
int write_volume_filetype(const volume<T>& source, const string& filename,
			 int filetype)
  { return save_volume_filetype(source,filename,filetype); }

template <class T>
int write_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype, const volumeinfo& vinfo)
  { return save_volume4D_filetype(source,filename,filetype,vinfo); }

template <class T>
int write_volume4D_filetype(const volume4D<T>& source, const string& filename,
			   int filetype)
  { return save_volume4D_filetype(source,filename,filetype); }


}

#endif

