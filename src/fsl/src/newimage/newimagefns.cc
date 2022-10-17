/*  newimagefns.cc

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

// General image processing functions


#include "newimagefns.h"



using namespace MISCMATHS;

namespace NEWIMAGE {

  ///////////////////////////////////////////////////////////////////////////

  // BASIC IMAGE SUPPORT FUNCTIONS




volume4D<float> sqrt(const volume4D<char>& vol4)
  {
    volume4D<float> retvol;
    retvol=sqrt_float(vol4);
    return retvol;
  }

 volume4D<float> sqrt(const volume4D<short>& vol4)
  {
    volume4D<float> retvol;
    retvol=sqrt_float(vol4);
    return retvol;
  }


 volume4D<float> sqrt(const volume4D<int>& vol4)
  {
    volume4D<float> retvol;
    retvol=sqrt_float(vol4);
    return retvol;
  }

  volume4D<float> sqrt(const volume4D<float>& vol4)
  {
    volume4D<float> retvol;
    retvol=sqrt_float(vol4);
    return retvol;
  }

 volume4D<double> sqrt(const volume4D<double>& vol4)
  {
    if (vol4.mint()<0) { volume4D<double> newvol; return newvol; }
    volume4D<double> retvol;
    copyconvert(vol4,retvol);
    for (int t=vol4.mint(); t<=vol4.maxt(); t++) {
      for (int z=vol4.minz(); z<=vol4.maxz(); z++) {
	for (int y=vol4.miny(); y<=vol4.maxy(); y++) {
	  for (int x=vol4.minx(); x<=vol4.maxx(); x++) {
	    if (vol4(x,y,z,t)>0) {
	      retvol(x,y,z,t) = sqrt((double) vol4(x,y,z,t));
	    } else {
	      retvol(x,y,z,t) = 0;
	    }
	  }
	}
      }
    }
    return retvol;
  }


 volume<float> sqrt(const volume<char>& vol)
  {
    volume<float> retvol;
    retvol=sqrt_float(vol);
    return retvol;
  }

 volume<float> sqrt(const volume<short>& vol)
  {
    volume<float> retvol;
    retvol=sqrt_float(vol);
    return retvol;
  }


 volume<float> sqrt(const volume<int>& vol)
  {
    volume<float> retvol;
    retvol=sqrt_float(vol);
    return retvol;
  }

  volume<float> sqrt(const volume<float>& vol)
  {
    volume<float> retvol;
    retvol=sqrt_float(vol);
    return retvol;
  }

  volume<double> sqrt(const volume<double>& vol)
  {
    volume<double> retvol;
    copyconvert(vol,retvol);
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (vol(x,y,z)>0) {
	    retvol(x,y,z) = sqrt((double) vol(x,y,z));
	  } else {
	    retvol(x,y,z) = 0;
	  }
	}
      }
    }
    return retvol;
  }

  float length(float x, float y) { return sqrt(x*x + y*y); }

  volume<float> abs(const volume<float>& realvol, const volume<float>& imagvol)
    {
      volume<float> absmap;
      absmap = realvol;
      for (int z=realvol.minz(); z<=realvol.maxz(); z++) {
	for (int y=realvol.miny(); y<=realvol.maxy(); y++) {
	  for (int x=realvol.minx(); x<=realvol.maxx(); x++) {
	    absmap(x,y,z) = length(imagvol(x,y,z),realvol(x,y,z));
	  }
	}
      }
      return absmap;
    }
  
  volume<float> phase(const volume<float>& realvol, 
		      const volume<float>& imagvol)
    {
      volume<float> phasemap;
      phasemap = realvol;
      for (int z=realvol.minz(); z<=realvol.maxz(); z++) {
	for (int y=realvol.miny(); y<=realvol.maxy(); y++) {
	  for (int x=realvol.minx(); x<=realvol.maxx(); x++) {
	    phasemap(x,y,z) = atan2(imagvol(x,y,z),realvol(x,y,z));
	  }
	}
      }
      return phasemap;
    }


  volume<float> real(const volume<float>& absvol, 
		      const volume<float>& phasevol)
    {
      volume<float> realmap;
      realmap = absvol;
      for (int z=absvol.minz(); z<=absvol.maxz(); z++) {
	for (int y=absvol.miny(); y<=absvol.maxy(); y++) {
	  for (int x=absvol.minx(); x<=absvol.maxx(); x++) {
	    realmap(x,y,z) = absvol(x,y,z) * cos(phasevol(x,y,z));
	  }
	}
      }
      return realmap;
    }


  volume<float> imag(const volume<float>& absvol, 
		      const volume<float>& phasevol)
    {
      volume<float> imagmap;
      imagmap = absvol;
      for (int z=absvol.minz(); z<=absvol.maxz(); z++) {
	for (int y=absvol.miny(); y<=absvol.maxy(); y++) {
	  for (int x=absvol.minx(); x<=absvol.maxx(); x++) {
	    imagmap(x,y,z) = absvol(x,y,z) * sin(phasevol(x,y,z));
	  }
	}
      }
      return imagmap;
    }

  ///////////////////////////////////////////////////////////////////////////

  // IMAGE PROCESSING ROUTINES

  void make_grad_masks(volume<float>& maskx, volume<float>& masky, 
		       volume<float>& maskz)
    {
      maskx.reinitialize(3,3,3);
      masky.reinitialize(3,3,3);
      maskz.reinitialize(3,3,3);
      for (int z=0; z<3; z++) {
	for (int y=0; y<3; y++) {
	  for (int x=0; x<3; x++) {
	    maskx(x,y,z)=(x-1.0)*pow(3.0,1.0-fabs(y-1.0)-fabs(z-1.0));
	    masky(x,y,z)=(y-1.0)*pow(3.0,1.0-fabs(x-1.0)-fabs(z-1.0));
	    maskz(x,y,z)=(z-1.0)*pow(3.0,1.0-fabs(x-1.0)-fabs(y-1.0));
	  }
	}
      }
      return;
    }

  void make_blur_mask(ColumnVector& bmask, const float final_vox_dim, 
		     const float init_vox_dim)
    {
      // construct the default output
      bmask.ReSize(1);
      bmask = 1.0;
      if (fabs(init_vox_dim)<1e-8) { return; }

      float sampling_ratio = final_vox_dim / init_vox_dim;
      if (sampling_ratio < 1.1) { return; }

      float sigma = 0.85*(sampling_ratio/2.0);
      if (sigma<0.5) { return; }

      int n=((int) (sigma-0.001))*2 + 3;
      int midn = n/2 + 1;
      bmask.ReSize(n);
      for (int x=1; x<=n; x++) {
	bmask(x) = exp(-((float) Sqr(x-midn))/( Sqr(sigma) * 4.0));
      }
      bmask = bmask / Sum(bmask);
      return;
    }


  ColumnVector gaussian_kernel1D(float sigma, int radius)
    {
      ColumnVector kern(2*radius+1);
      float sum=0.0, val=0.0;
      
      for(int j=-radius; j<=radius; j++) {
	if (sigma>1e-6) {
	  val = exp(-(j*j)/(2.0*sigma*sigma));
	} else {
	  if (j==0) { val=1; } else { val=0; }
	}
	kern(j+radius+1) = val;
	sum += val;
      }
      
      kern *= (1.0/sum);
      return kern;
    }


  volume<float> gaussian_kernel2D(float sigma, int radius)
    {
      volume<float> new_kernel((2*radius+1),(2*radius+1),1); 
      float sum=0.0, val=0.0;
      
      for(int i=-radius; i<=radius; i++) {
	for(int j=-radius; j<=radius; j++) {
	  if (sigma>1e-6) {
	    val = exp(-(i*i+j*j)/(2.0*sigma*sigma));
	  } else {
	    if ((i*i + j*j)==0) { val=1; } else { val=0; }
	  }
	  new_kernel((j+radius),(i+radius),0) = val;
	  sum += val;
	}
      }
      
      new_kernel *= (1.0/sum);
      return new_kernel;
    }


  volume<float> gaussian_kernel3D(float sigma, int radius)
    {
      volume<float> new_kernel((2*radius+1),(2*radius+1),(2*radius+1)); 
      float sum=0.0, sum2=0.0, val=0.0;
      
      for(int i=-radius; i<=radius; i++) {
	for(int j=-radius; j<=radius; j++) {
	  for(int k=-radius; k<=radius; k++) {
	    if (sigma>1e-6) {
	      val = exp(-(i*i+j*j+k*k)/(2.0*sigma*sigma));
	    } else {
	      if ((i*i + j*j + k*k)==0) { val=1; } else { val=0; }
	    }
	    new_kernel((j+radius),(i+radius),(k+radius)) = val;
	    sum += val;
	  }
	}
	sum2 += sum; sum=0.0;
      }
      
      new_kernel *= (1.0/sum2);
      return new_kernel;
    }




  volume<float> gaussian_kernel3D(float sigma, float xdim, float ydim, float zdim) {
  int sx = ((int) ceil(sigma*4.0/xdim))*2 + 1;
  int sy = ((int) ceil(sigma*4.0/ydim))*2 + 1;
  int sz = ((int) ceil(sigma*4.0/zdim))*2 + 1;
  volume<float> vker(sx,sy,sz);
  float dx2=Sqr(xdim);
  float dy2=Sqr(ydim);
  float dz2=Sqr(zdim);
  for (int z=-sz/2; z<=sz/2; z++) {
    for (int y=-sy/2; y<=sy/2; y++) {
      for (int x=-sx/2; x<=sx/2; x++) {
	vker(x+sx/2,y+sy/2,z+sz/2)=exp(-(x*x*dx2+y*y*dy2+z*z*dz2)/(2*sigma*sigma));
      }
    }
  }
  return vker;
  }


  volume<float> spherical_kernel(float radius, float xdim, float ydim, float zdim)
  {
  int sx = MISCMATHS::round(radius/xdim)*2 + 1;
  int sy = MISCMATHS::round(radius/ydim)*2 + 1;
  int sz = MISCMATHS::round(radius/zdim)*2 + 1;
  volume<float> vker(sx,sy,sz);
  vker = 0.0;
  float dx2=Sqr(xdim);
  float dy2=Sqr(ydim);
  float dz2=Sqr(zdim);
  for (int z=-sz/2; z<=sz/2; z++) {
    for (int y=-sy/2; y<=sy/2; y++) {
      for (int x=-sx/2; x<=sx/2; x++) {
	if ((x*x*dx2+y*y*dy2+z*z*dz2)<=Sqr(radius)) { 
	  vker(x+sx/2,y+sy/2,z+sz/2)=1.0; 
	}
      }
    }
  }
  return vker;
  }

  volume<float> box_kernel(float length, float xdim,float ydim,float zdim)  //mm dimensions
  {
      int x = ((int) floor(length/xdim/2))*2 + 1;
      int y = ((int) floor(length/ydim/2))*2 + 1;
      int z = ((int) floor(length/zdim/2))*2 + 1;
      volume<float> new_kernel(x,y,z);
      new_kernel=1.0;
      return new_kernel;          
  }




  volume<float> box_kernel(int x,int y, int z)  //voxel dimensions
  {
      volume<float> new_kernel(x,y,z);
      new_kernel=1.0;
      return new_kernel;          
  }








float fsllog2(float x)
{
  // a cygwin annoyance!
  return log(x)/log(2);
}




  ///////////////////////////////////////////////////////////////////////////

  // support functions for connected components

  int find_first_nonzero(const Matrix& mat)
    {
      Tracer tr("first");
      for (int idx=1; idx<=mat.Nrows(); idx++) {
	if (mat(idx,1)!=0.0) return idx;
      }
      return -1;  // Failed
    }

  void addpair2set(int x, int y, std::vector<int>& sx, std::vector<int>& sy)
    {
      sx.push_back(x);
      sy.push_back(y);
    }


  inline void get_parent_label(int& idx, const Matrix& idxmap) 
  {
    while (idxmap(idx,1)>0.0) { idx = MISCMATHS::round(float(idxmap(idx,1))); }
  }


  void relabel_components_uniquely(volume<int>& labelvol, 
				   const std::vector<int>& equivlista,
				   const std::vector<int>& equivlistb, ColumnVector& clustersizes) 
  {
    int labelnum = labelvol.max();
    Matrix emap(labelnum,1);
    emap = -0.2;

    int n1, n2;
    for (unsigned int n=0; n<equivlista.size(); n++) {
      n1 = equivlista[n];
      get_parent_label(n1,emap);
      n2 = equivlistb[n];
      get_parent_label(n2,emap);
      if (n1!=n2) emap(Max(n1,n2),1) = Min(n1,n2);
    }

    // re-parse emap to assign sequential, unique numbers
    int newlabel=1;
    for (int n=1; n<=labelnum; n++) {
      int n1 = n;
      get_parent_label(n1,emap);
      if (n1<n) {  // it points to another label
	emap(n,1) = emap(n1,1);
      } else {  // it is a newly found label
	emap(n,1) = -newlabel;
	newlabel++;
      }
    }
    
    int numclusts=newlabel-1;
    clustersizes.ReSize(numclusts);
    clustersizes=0;
    
    // Change the old labels to new ones
    
    for (int z=labelvol.minz(); z<=labelvol.maxz(); z++) {
      for (int y=labelvol.miny(); y<=labelvol.maxy(); y++) {
	for (int x=labelvol.minx(); x<=labelvol.maxx(); x++) {
	  if (labelvol(x,y,z)>0) {

	    int tmp = MISCMATHS::round(-float(emap(labelvol(x,y,z),1)));
	    labelvol(x,y,z)=tmp;
	    clustersizes(tmp)+=1;
	    
	  }
	}
      }
    }
  }

  void relabel_components_uniquely(volume<int>& labelvol, 
				   const std::vector<int>& equivlista,
				   const std::vector<int>& equivlistb){
    ColumnVector clustersize;
    relabel_components_uniquely(labelvol, equivlista, equivlistb,clustersize); 
    
    
  }
















  ///////////////////////////////////////////////////////////////////////////

}

