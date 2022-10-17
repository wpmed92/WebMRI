/*  General IO functions (images and transformation files)

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999 University of Oxford  */

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

#include "newimageio.h"

using namespace MISCMATHS;

namespace NEWIMAGE {



////////////////////////////////////////////////////////////////////////////

// VOLUME I/O

void WriteClonedHeader(FSLIO *dest, const FSLIO *src)
{
  FslCloneHeader(dest,src);
  FslSetIntensityScaling(dest,1.0,0.0);
}

mat44 newmat2mat44(const Matrix& nmat)
{
  mat44 ret;
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      ret.m[i-1][j-1] = nmat(i,j);
    }
  }
  return ret;
}


string fslbasename(const string& filename)
{
  // does string() copy safely and dispose of temporary c-string storage?
  return string(FslMakeBaseName(filename.c_str()));
}


int make_basename(string& filename)
{
  char *tmpname;
  tmpname = FslMakeBaseName(filename.c_str());
  if (tmpname==NULL) return -1;
  // this is now the basename
  filename = string(tmpname);
  // free(tmpname);  // is this safe to do?
  return 0;
}


int find_pathname(string& filename)
{
  Tracer tr("find_pathname");
  if (filename.size() < 1) return -1;
  string pathname = filename;
  int fsize = pathname.length(), indx;

  // working backwards, find '/' and remove everything after it

  indx = fsize-1;
  while ((pathname[indx] != '/') && (indx != 0))
    indx--;
  
  if (indx<fsize-1)
    pathname.erase(indx+1);
  
  filename = pathname;
  return 0;
}


bool fsl_imageexists(const string& filename) {
  return (FslFileExists(filename.c_str())!=0);
}


void check_filename(const string& basename)
{
  FSLIO* OP=FslOpen(basename.c_str(),"r");
  if (OP==NULL) {
    cerr << "ERROR: Cannot open volume " << basename << " for reading!\n";
    exit(1);
  }
}

FSLIO* NewFslOpen(const string& filename, const string& permissions, 
		  int filetype, const volumeinfo& vinfo, bool use_vinfo)
{
  string basename = filename;
  make_basename(basename);
  if ( basename.size()<1 ) {
    return 0;
  }

  bool writemode=false;
  if ( (permissions.find('w')!=string::npos) || 
       (permissions.find('+')!=string::npos) )  { writemode=true; }

  FSLIO* OP=FslXOpen(basename.c_str(),permissions.c_str(),filetype);
  if (OP==NULL) {
    cerr << "ERROR: Could not open image " << basename << endl;
    return NULL;
  }

  if (use_vinfo) {
    if (writemode) { 
      WriteClonedHeader(OP,&vinfo); 
    } else { 
      FSLIO* tmp = (FSLIO *) &vinfo;  
      FslCloneHeader(tmp,OP);   // cast away const
    }
  }
  return OP;
}

FSLIO* NewFslOpen(const string& filename, const string& permissions, 
		  const volumeinfo& vinfo, bool use_vinfo)
{
  return NewFslOpen(filename,permissions,-1,vinfo,use_vinfo);
}


FSLIO* NewFslOpen(const string& filename, const string& permissions)
{
  volumeinfo vinfo = blank_vinfo();
  return NewFslOpen(filename,permissions,vinfo,false);
}



short dtype(const char* T)   { return DT_UNSIGNED_CHAR; }
short dtype(const short* T)  { return DT_SIGNED_SHORT; }
short dtype(const int* T)    { return DT_SIGNED_INT; }
short dtype(const float* T)  { return DT_FLOAT; }
short dtype(const double* T) { return DT_DOUBLE; }

short dtype(const volume<char>& vol)   { return DT_UNSIGNED_CHAR; }
short dtype(const volume<short>& vol)  { return DT_SIGNED_SHORT; }
short dtype(const volume<int>& vol)    { return DT_SIGNED_INT; }
short dtype(const volume<float>& vol)  { return DT_FLOAT; }
short dtype(const volume<double>& vol) { return DT_DOUBLE; }

short dtype(const volume4D<char>& vol)   { return DT_UNSIGNED_CHAR; }
short dtype(const volume4D<short>& vol)  { return DT_SIGNED_SHORT; }
short dtype(const volume4D<int>& vol)    { return DT_SIGNED_INT; }
short dtype(const volume4D<float>& vol)  { return DT_FLOAT; }
short dtype(const volume4D<double>& vol) { return DT_DOUBLE; }

short dtype(const string& filename) 
{
  Tracer trcr("dtype");
  if ( filename.size()<1 ) return -1;
  string basename = fslbasename(filename);

  FSLIO* IP1;
  IP1 = FslOpen(basename.c_str(),"rb");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
    exit(1);
  }

  short dtype;
  FslGetDataType(IP1,&dtype);
  float slope, intercept;
  int doscaling;
  doscaling = FslGetIntensityScaling(IP1,&slope,&intercept);
  if (doscaling==1) { dtype = DT_FLOAT; }
  FslClose(IP1);
  free(IP1);

  return dtype;
}

volumeinfo blank_vinfo() 
{
  volumeinfo *vinfo;
  vinfo = FslInit();
  return *vinfo;
}

volumeinfo volinfo(const string& filename) 
{
  Tracer trcr("volinfo");
  volumeinfo vinfo = blank_vinfo();
  if ( filename.size()<1 ) return vinfo;
  string basename = filename;
  make_basename(basename);

  FSLIO* IP1 = FslOpen(basename.c_str(), "r");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
    exit(1);
  }

  FslCloneHeader(&vinfo,IP1);
  FslClose(IP1);

  return vinfo;
}

//////////////////////////////////////////////////////////////////////////

// COMPLEX IMAGE I/O

void FslReadComplexBuffer(FSLIO* IP, float* realbuffer, float* imagbuffer) 
{
  short sx,sy,sz,st;
  FslGetDim(IP,&sx,&sy,&sz,&st);
  size_t imagesize=sx*sy*sz;
  short type;
  FslGetDataType(IP,&type);
  switch(type)
    {
    case DT_COMPLEX:
      {
	float* sbuffer=new float[2*imagesize];
	if (sbuffer==0) { imthrow("Out of memory",99); }
	FslReadVolumes(IP,sbuffer,1);
	float *sptr=sbuffer, *rptr=realbuffer, *iptr=imagbuffer;
	for (size_t poff=0; poff<imagesize; poff++) {
	  *rptr++ = *sptr++;
	  *iptr++ = *sptr++;
	}
	delete[] sbuffer;
      }
      break;
    default:
      {  // read in the real part instead
	FslReadBuffer<float>(IP,realbuffer);
	float *iptr=imagbuffer;
	for (size_t poff=0; poff<imagesize; poff++) {
	  *iptr++ = 0.0;
	}
      }
    }
}
  
void FslWriteComplexVolume(FSLIO* OP, const float* realbuffer,
			    const float* imagbuffer) 
{
  short sx,sy,sz,st;
  FslGetDim(OP,&sx,&sy,&sz,&st);
  size_t imagesize=sx*sy*sz*1;
  float* sbuffer=new float[2*imagesize];
  if (sbuffer==0) { imthrow("Out of memory",99); }
  float *sptr=sbuffer;
  const float *rptr=realbuffer, *iptr=imagbuffer;
  for (size_t poff=0; poff<imagesize; poff++) {
    *sptr++ = *rptr++;
    *sptr++ = *iptr++;
  }
  FslWriteVolumes(OP,sbuffer,1);
  delete[] sbuffer;
}


int read_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename, volumeinfo& vinfo,
		       bool read_img_data)
{
  Tracer trcr("read_complexvolume");
  if ( filename.size()<1 ) return -1;
  string basename = filename;
  make_basename(basename);

  FSLIO* IP1 = FslOpen(basename.c_str(), "r");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
    exit(1);
  }
  short sx,sy,sz,st;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  size_t volsize=sx*sy*sz;

  float* realbuffer=new float[volsize];
  if (realbuffer==0) { imthrow("Out of memory",99); }
  float* imagbuffer=new float[volsize];
  if (imagbuffer==0) { imthrow("Out of memory",99); }
  if (read_img_data)  FslReadComplexBuffer(IP1,realbuffer,imagbuffer);
  realvol.reinitialize(sx,sy,sz,realbuffer,true);
  imagvol.reinitialize(sx,sy,sz,imagbuffer,true);

  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);
  realvol.setdims(x,y,z);
  imagvol.setdims(x,y,z);

  vinfo = blank_vinfo();
  FslCloneHeader(&vinfo,IP1);
  FslClose(IP1);
  return 0;
}


int read_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename)
{
  volumeinfo vinfo = blank_vinfo();
  int retval = read_complexvolume(realvol,imagvol,filename,vinfo,true);
  return retval;
}

int read_complexvolume(complexvolume& vol, const string& filename)
{
  return read_complexvolume(vol.re(),vol.im(),filename);
}

int read_complexvolume(complexvolume& vol, const string& filename, 
		       volumeinfo& vinfo)
{
  return read_complexvolume(vol.re(),vol.im(),filename,vinfo,true);
}


int read_complexvolume4D(volume4D<float>& realvols, volume4D<float>& imagvols,
			 const string& filename, volumeinfo& vinfo,
			 bool read_img_data)
{
  Tracer trcr("read_complexvolume4D");
  if ( filename.size()<1 ) return -1;
  string basename = filename;
  make_basename(basename);

  FSLIO* IP1 = FslOpen(basename.c_str(), "r");
  if (IP1==NULL) {
    cerr << "Cannot open volume " << basename << " for reading!\n";
    exit(1);
  }
  short sx,sy,sz,st;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  size_t volsize=sx*sy*sz;
  if (st<1) st=1;   //make it robust to dim4<1

  volume<float> dummyvol(sx,sy,sz);
  for (int t=0; t<st; t++) {
    realvols.addvolume(dummyvol);
    imagvols.addvolume(dummyvol);
    float* rbuffer=new float[volsize];
    if (rbuffer==0) { imthrow("Out of memory",99); }
    float* ibuffer=new float[volsize];
    if (ibuffer==0) { imthrow("Out of memory",99); }
    if (read_img_data)  FslReadComplexBuffer(IP1,rbuffer,ibuffer);
    // Note that the d_owner flag = true in the following so that the
    //  control for delete is passed to the volume class
    realvols[t].reinitialize(sx,sy,sz,rbuffer,true);
    imagvols[t].reinitialize(sx,sy,sz,ibuffer,true);
  }

  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);
  realvols.setdims(x,y,z,tr);
  imagvols.setdims(x,y,z,tr);

  vinfo = blank_vinfo();
  FslCloneHeader(&vinfo,IP1);
  FslClose(IP1);
  return 0;
}


int read_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename)
{
  volumeinfo vinfo = blank_vinfo();
  int retval = read_complexvolume4D(realvol,imagvol,filename,vinfo,true);
  return retval;
}

int read_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename, volumeinfo& vinfo)
{
  int retval = read_complexvolume4D(realvol,imagvol,filename,vinfo,true);
  return retval;
}

/*
int read_complexvolume4D(complexvolume& vol, const string& filename)
{
  return read_complexvolume4D(vol.re(),vol.im(),filename);
}

int read_complexvolume4D(complexvolume &vol, const string& filename, 
			 volumeinfo& vinfo)
{
  return read_complexvolume4D(vol.re(),vol.im(),filename,vinfo);
}
*/

int save_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename,
		       const volumeinfo& vinfo, bool use_vinfo)
{
  Tracer tr("save_complexvolume");
  string basename = filename;
  make_basename(basename);
  if ( basename.size()<1 ) return -1;

  FSLIO* OP=FslOpen(basename.c_str(),"w");
  if (OP==0) return -1;
  if (use_vinfo) WriteClonedHeader(OP,&vinfo);
    
  FslSetDim(OP,realvol.xsize(),realvol.ysize(),realvol.zsize(),1);
  FslSetDataType(OP, DT_COMPLEX);
  FslSetVoxDim(OP,realvol.xdim(), realvol.ydim(), realvol.zdim(),1.0);
  
  FslWriteHeader(OP);
  FslWriteComplexVolume(OP,&(realvol(0,0,0)),&(imagvol(0,0,0)));

  FslClose(OP); 
  return 0;
}


int save_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename)
{  
  volumeinfo vinfo = blank_vinfo();
  return save_complexvolume(realvol,imagvol,filename,vinfo,false);
}


int save_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename,
		       const volumeinfo& vinfo)
{  
  return save_complexvolume(realvol,imagvol,filename,vinfo,true);
}

int save_complexvolume(const complexvolume& vol, const string& filename)
{
  return save_complexvolume(vol.re(),vol.im(),filename);
}

int save_complexvolume(const complexvolume& vol, const string& filename,
		       const volumeinfo& vinfo)
{
  return save_complexvolume(vol.re(),vol.im(),filename,vinfo);
}

int save_complexvolume4D(const volume4D<float>& realvols, 
			 const volume4D<float>& imagvols, 
			 const string& filename,
			 const volumeinfo& vinfo, bool use_vinfo)
{
  Tracer tr("save_complexvolume4D");

  if (realvols.tsize()<=0) return -1;

  string basename = filename;
  make_basename(basename);
  if ( basename.size()<1 ) return -1;

  FSLIO* OP=FslOpen(basename.c_str(),"w");
  if (OP==0) return -1;
  if (use_vinfo) WriteClonedHeader(OP,&vinfo);
    
  FslSetDim(OP,realvols.xsize(),realvols.ysize(),realvols.zsize(),
	    realvols.tsize());
  FslSetDataType(OP, DT_COMPLEX);
  FslSetVoxDim(OP,realvols.xdim(), realvols.ydim(), realvols.zdim(), 
	       realvols.tdim());

  FslWriteHeader(OP);
  
  for (int t=0; t<realvols.tsize(); t++) {
    FslWriteComplexVolume(OP,&(realvols[t](0,0,0)),&(imagvols[t](0,0,0)));
  }

  FslClose(OP); 
  return 0;
}


int save_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, const string& filename)
{  
  volumeinfo vinfo = blank_vinfo();
  return save_complexvolume4D(realvol,imagvol,filename,vinfo,false);
}


int save_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, const string& filename,
			 const volumeinfo& vinfo)
{  
  return save_complexvolume4D(realvol,imagvol,filename,vinfo,true);
}

/*
int save_complexvolume4D(const complexvolume& vol, const string& filename)
{
  return save_complexvolume4D(vol.re(),vol.im(),filename);
}

int save_complexvolume4D(const complexvolume& vol, const string& filename,
			 const volumeinfo& vinfo)
{
  return save_complexvolume4D(vol.re(),vol.im(),filename,vinfo);
}
*/


//////////////////////////////////////////////////////////////////////////

// MATRIX I/O


int read_ascii_matrix(Matrix &target, const string& filename)
{
  Tracer tr("read_ascii_matrix");
  target=MISCMATHS::read_ascii_matrix(filename);
  if (target.Nrows()<=0) return -1;
  return 0;
}


int get_medx_small_matrix(Matrix &target, ifstream& matfile)
{
  Tracer tr("get_medx_small_matrix");
  string str1;
  matfile >> str1;
  if (str1 != "[") {
    return -1;
  }
  int i=1, j=1;
  matfile >> str1;
  while (str1 != "]") {
    target(i,j) = atof(str1.c_str());
    if (j==4) { i++; j=1; }
    else { j++; }
    matfile >> str1;
  }
  return 0;
}  


int get_medx_matrix(Matrix &target, ifstream& matfile)
{
  Tracer tr("get_medx_matrix");
  string str1, str2;
  matfile >> str1 >> str2;
  if ((str1 != "<<") || (str2 != "/matrix")) {
    return -1;
  }
  target.ReSize(4,4);
  Identity(target);  // default return value
  return get_medx_small_matrix(target,matfile);
}  


int get_minc_matrix(Matrix &target, ifstream& matfile)
{
  Tracer tr("get_minc_matrix");
  string str1;
  matfile >> str1;
  if (str1 != "=") {
    cerr << "Could not parse MINC transform file" << endl;
    return -1;
  }
  target.ReSize(4,4);
  Identity(target);  // default return value
  for (int i=1; i<=3; i++) {
    for (int j=1; j<=4; j++) {
      matfile >> str1;
      target(i,j) = atof(str1.c_str());
    }
  }
  return 0;
}  


//------------------------------------------------------------------------//

int put_medx_matrix(ofstream& matfile, const string& name, const Matrix& affmat)
{
  Tracer tr("put_medx_matrix");
  if (affmat.Nrows()<=0) { return -1; }
  matfile << "        /" << name << " [" << endl;
  for (int i=1; i<=affmat.Nrows(); i++) {
    for (int j=1; j<=affmat.Ncols(); j++) {
      matfile << "            " << affmat(i,j) << endl;
    }
  }
  matfile << "        ]" << endl;
  return 0;
}



int get_outputusermat(const string& filename, Matrix& oumat)
{
  Tracer trcr("get_outputusermat");
  if ( filename.size()<1 ) return -1;
  string basename = filename, fname;
  make_basename(basename);

  // check that the .hdr file exists
  if (!fsl_imageexists(basename)) {
    cerr << "Cannot open volume " << fname << " for reading!\n";
    exit(1);
  }

  FSLIO* IP1(FslOpen(basename.c_str(), "r"));

  float x,y,z,tr;
  FslGetVoxDim(IP1,&x,&y,&z,&tr);

  ColumnVector origin(3);
  origin=0.0;
  // read origin from sform (if set)
  mat44 sform_mat;
  short sform_code = FslGetStdXform(IP1,&sform_mat);
  if (sform_code!=NIFTI_XFORM_UNKNOWN) {
    sform_mat = mat44_inverse(sform_mat);
    origin(1) = sform_mat.m[0][3];
    origin(2) = sform_mat.m[1][3];
    origin(3) = sform_mat.m[2][3];
  }

  short sx,sy,sz,st;
  FslGetDim(IP1,&sx,&sy,&sz,&st);
  origin(2) = -sy - origin(2);   // Convert to dumb MEDx conventions

  oumat.ReSize(4,4);
  Identity(oumat);
  oumat(1,1) = x;
  oumat(2,2) = -y;
  oumat(3,3) = z;

  oumat(1,4) = origin(1)*x;
  oumat(2,4) = -origin(2)*y;
  oumat(3,4) = origin(3)*z;

  FslClose(IP1);
  return 0;
}  


//------------------------------------------------------------------------//


int load_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename)
  { return read_complexvolume(realvol,imagvol,filename); }

int load_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename, volumeinfo& vinfo)
  { return read_complexvolume(realvol,imagvol,filename,vinfo,true); }

int load_complexvolume(complexvolume& vol, const string& filename)
  { return read_complexvolume(vol,filename); }

int load_complexvolume(complexvolume& vol, const string& filename, 
		       volumeinfo& vinfo)
  { return read_complexvolume(vol,filename,vinfo); }

int load_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename)
  { return read_complexvolume4D(realvol,imagvol,filename); }

int load_complexvolume4D(volume4D<float>& realvol, volume4D<float>& imagvol,
			 const string& filename, volumeinfo& vinfo)
  { return read_complexvolume4D(realvol,imagvol,filename,vinfo,true); }

/*
int load_complexvolume4D(complexvolume& vol, const string& filename)
  { return read_complexvolume4D(vol,filename); }

int load_complexvolume4D(complexvolume &vol, const string& filename, 
			 volumeinfo& vinfo)
  { return read_complexvolume4D(vol,filename,vinfo); }
*/
int write_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename)
  { return save_complexvolume(realvol,imagvol,filename); }

int write_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename,
		       const volumeinfo& vinfo)
  { return save_complexvolume(realvol,imagvol,filename,vinfo,true); }

int write_complexvolume(const complexvolume& vol, const string& filename)
  { return save_complexvolume(vol,filename); }

int write_complexvolume(const complexvolume& vol, const string& filename,
		       const volumeinfo& vinfo)
  { return save_complexvolume(vol,filename,vinfo); }

int write_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename)
  { return save_complexvolume4D(realvol,imagvol,filename); }

int write_complexvolume4D(const volume4D<float>& realvol, 
			 const volume4D<float>& imagvol, 
			 const string& filename,
			 const volumeinfo& vinfo)
  { return save_complexvolume4D(realvol,imagvol,filename,vinfo,true); }

/*
int write_complexvolume4D(const complexvolume& vol, const string& filename)
  { return save_complexvolume4D(vol,filename); }

int write_complexvolume4D(const complexvolume& vol, const string& filename,
			 const volumeinfo& vinfo)
  { return save_complexvolume4D(vol,filename,vinfo); }
*/

}








