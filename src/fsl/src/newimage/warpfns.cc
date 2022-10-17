/*  warpfns.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

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


// *** NOTE *** //
// MUST DEAL WITH THE PROBLEM OF EXTRAPOLATING VALUES FOR WARPS AS IT IS
//  NOT GUARANTEED THAT THE FOV IS SUFFICIENT

#include "warpfns.h"

namespace NEWIMAGE {

////////////////////////////////////////////////////////////////////////////

int affine2warp(const Matrix& affmat, volume4D<float>& warpvol,
		const volume<float>& outvol)
{
  if (outvol.nvoxels() <= 0) {
    cerr << "Cannot do affine2warp as outvol has no size" << endl;
    return -1;
  }
  warpvol.reinitialize(outvol.xsize(),outvol.ysize(),outvol.zsize(),3);
  warpvol[0] = outvol;
  warpvol[1] = outvol;
  warpvol[2] = outvol;

  ColumnVector xin(4), xout(4);
  xin(4) = 1.0;  xout(4)=1.0;

  for (int z=outvol.minz(); z<=outvol.maxz(); z++) {
    for (int y=outvol.miny(); y<=outvol.maxy(); y++) {
      for (int x=outvol.minx(); x<=outvol.maxx(); x++) {
	//   convert x,y,z to mm coords (xout)
	xout(1) = x;  xout(2) = y;  xout(3) = z;
	xout = outvol.sampling_mat() * xout;
	xin = affmat.i() * xout;
	// use the mm coordinates to store the results
	warpvol[0](x,y,z) = xin(1);
	warpvol[1](x,y,z) = xin(2);
	warpvol[2](x,y,z) = xin(3);
      }
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int calc_dir(const string& shiftdir, int& dir, int& sign)
{
  // dir is 0,1,2 for x,y,z  :  sign is +/- 1 for x/x-  etc..
  if (shiftdir=="x") {
    dir=0;  sign=1;
  } else if (shiftdir=="y") {
    dir=1;  sign=1;
  } else if (shiftdir=="z") {
    dir=2;  sign=1;
  } else if (shiftdir=="x-") {
    dir=0;  sign=-1;
  } else if (shiftdir=="y-") {
    dir=1;  sign=-1;
  } else if (shiftdir=="z-") {
    dir=2;  sign=-1;
  } else {
    cerr << "Cannot interpret shift direction = " << shiftdir << endl;
    return -1;
  }
  return 0;
}


int shift2warp(const volume<float>& shiftmap, 
	       volume4D<float>& warp, const string& shiftdir)
{
  affine2warp(Identity(4),warp,shiftmap);  // use shiftmap as refvol (set size)
  int dir, sign;
  calc_dir(shiftdir,dir,sign);
  float voxdim = shiftmap.sampling_mat()(dir+1,dir+1);

  for (int z=shiftmap.minz(); z<=shiftmap.maxz(); z++) {
    for (int y=shiftmap.miny(); y<=shiftmap.maxy(); y++) {
      for (int x=shiftmap.minx(); x<=shiftmap.maxx(); x++) {
	// get amount of shift in mm
	float shift = shiftmap(x,y,z) * voxdim * sign;
	warp[dir](x,y,z) += shift;
      }
    }
  }
  return 0;
}


////////////////////////////////////////////////////////////////////////////

int convertwarp_rel2abs(volume4D<float>& warpvol)
{
  // conversion is: w(x) = x + u(x)  (all in mm)
  for (int z=0; z<warpvol.zsize(); z++) {
    for (int y=0; y<warpvol.ysize(); y++) {
      for (int x=0; x<warpvol.xsize(); x++) {
	warpvol(x,y,z,0) += x*warpvol.xdim();
	warpvol(x,y,z,1) += y*warpvol.ydim();
	warpvol(x,y,z,2) += z*warpvol.zdim();
      }
    }
  }
  return 0;
}

int convertwarp_abs2rel(volume4D<float>& warpvol)
{
  // conversion is: w(x) = x + u(x)  (all in mm)
  for (int z=0; z<warpvol.zsize(); z++) {
    for (int y=0; y<warpvol.ysize(); y++) {
      for (int x=0; x<warpvol.xsize(); x++) {
	warpvol(x,y,z,0) -= x*warpvol.xdim();
	warpvol(x,y,z,1) -= y*warpvol.ydim();
	warpvol(x,y,z,2) -= z*warpvol.zdim();
      }
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int concat_warps(const volume4D<float>& prewarp, 
                 const volume4D<float>& postwarp,
		 volume4D<float>& totalwarp)
{
  totalwarp = postwarp;  // set size
  totalwarp = 0.0;
  ColumnVector xmid(4), xpre(4);
  xmid(4) = 1.0;  xpre(4)=1.0;
  for (int z=postwarp.minz(); z<=postwarp.maxz(); z++) {
    for (int y=postwarp.miny(); y<=postwarp.maxy(); y++) {
      for (int x=postwarp.minx(); x<=postwarp.maxx(); x++) {
	xmid(1) = postwarp[0](x,y,z);
	xmid(2) = postwarp[1](x,y,z);
	xmid(3) = postwarp[2](x,y,z);
	// convert xmid from mm to voxels (of prewarp image)
	xmid = prewarp[0].sampling_mat().i() * xmid;
	// look up the coordinates in prewarp
	xpre(1) = prewarp[0].interpolate(xmid(1),xmid(2),xmid(3));
	xpre(2) = prewarp[1].interpolate(xmid(1),xmid(2),xmid(3));
	xpre(3) = prewarp[2].interpolate(xmid(1),xmid(2),xmid(3));
	// set these mm coordinates as the result
	totalwarp[0](x,y,z) = xpre(1);
	totalwarp[1](x,y,z) = xpre(2);
	totalwarp[2](x,y,z) = xpre(3);
      }
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int fast_apply_warp(const volume<float>& invol, volume<float>& outvol,
		    const volume4D<float>& warpvol)
{
  ColumnVector xin(4);
  xin(4)=1.0;
  float I_in;
  // assumes that warpvol has same number and size of voxels as invol
  for (int z=outvol.minz(); z<=outvol.maxz(); z++) {
    for (int y=outvol.miny(); y<=outvol.maxy(); y++) {
      for (int x=outvol.minx(); x<=outvol.maxx(); x++) {
	//   assume outvol and invol have same voxel size
	//   look up the warp dest coordinate (in mm)
	if (warpvol[0].in_bounds(x,y,z)) {
	  xin(1) = warpvol[0](x,y,z);
	  xin(2) = warpvol[1](x,y,z);
	  xin(3) = warpvol[2](x,y,z);
	  //   convert xin from mm to voxel coords
	  I_in = invol.interpolate(xin(1)/invol.xdim(),xin(2)/invol.ydim(),
				   xin(3)/invol.zdim());
	} else {
	  I_in = invol.getpadvalue();
	}
	outvol(x,y,z) = I_in;
      }
    }
  }
  return 0;
}

int raw_apply_warp(const volume<float>& invol, volume<float>& outvol,
		   const volume4D<float>& warpvol, 
		   const Matrix& premat, const Matrix& postmat)
{
  if (outvol.nvoxels() <= 0) {
    cerr << "Cannot apply warp to outvol as it has no size" << endl;
    return -1;
  }
  if ( (samesize(invol,warpvol[0])) && (invol.xdim() == warpvol.xdim() )
       && (invol.ydim() == warpvol.ydim() )
       && (invol.zdim() == warpvol.zdim() ) 
       && ( (premat - Identity(4)).MaximumAbsoluteValue() < 1e-3)
       && ( (postmat - Identity(4)).MaximumAbsoluteValue() < 1e-3) )
 {
    return fast_apply_warp(invol,outvol,warpvol);
  }
  ColumnVector xin(4), xout(4);
  xin(4) = 1.0;  xout(4)=1.0;
  float I_in;
  for (int z=outvol.minz(); z<=outvol.maxz(); z++) {
    for (int y=outvol.miny(); y<=outvol.maxy(); y++) {
      for (int x=outvol.minx(); x<=outvol.maxx(); x++) {
	//   convert x,y,z to mm coords (xout)
	xout(1) = x;  xout(2) = y;  xout(3) = z;
	xout = outvol.sampling_mat() * xout;
	//   apply inverse of postmat
	xout = postmat.i() * xout;
	//   convert xout to warpvol voxel coords
	xout = warpvol[0].sampling_mat().i() * xout;
	//   look up the warp dest coordinate (in mm)
	if (warpvol[0].in_bounds(MISCMATHS::round(xout(1)),
				 MISCMATHS::round(xout(2)),
				 MISCMATHS::round(xout(3)))) {
	  xin(1) = warpvol[0].interpolate(xout(1),xout(2),xout(3));
	  xin(2) = warpvol[1].interpolate(xout(1),xout(2),xout(3));
	  xin(3) = warpvol[2].interpolate(xout(1),xout(2),xout(3));
	  //   apply inverse of premat
	  xin = premat.i() * xin;
	  //   convert xin from mm to voxel coords
	  xin = invol.sampling_mat().i() * xin;
	  I_in = invol.interpolate(xin(1),xin(2),xin(3));
	} else {
	  I_in = invol.getpadvalue();
	}
	outvol(x,y,z) = I_in;
      }
    }
  }
  return 0;
}


int apply_warp(const volume<float>& invol, volume<float>& outvol,
	       const volume4D<float>& warpvol, 
	       const Matrix& premat, const Matrix& postmat)
{
  // set the desired extrapolation settings
  extrapolation oldin = invol.getextrapolationmethod();
  extrapolation oldwarp = warpvol.getextrapolationmethod();
  warpvol.setextrapolationmethod(extraslice);
  invol.setextrapolationmethod(extraslice);
  float oldpad = invol.getpadvalue();
  invol.setpadvalue(invol.backgroundval());

  int retval = raw_apply_warp(invol,outvol,warpvol,premat,postmat);

  // restore extrapolation settings
  warpvol.setextrapolationmethod(oldwarp);
  invol.setextrapolationmethod(oldin);
  invol.setpadvalue(oldpad);
  
  return retval;
}



int apply_warp(const volume<float>& invol, volume<float>& outvol,
	       const volume4D<float>& warpvol)
{
  Matrix ident(4,4);
  ident = Identity(4);
  return apply_warp(invol,outvol,warpvol,ident,ident);
}


}

