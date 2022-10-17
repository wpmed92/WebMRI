/*  avwconv.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2003-2004 University of Oxford  */

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

// General convolution and filtering utility

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;


string title="avwconv (Version 1.0)\nCopyright(c) 2003, University of Oxford (Mark Jenkinson)";
string examples="avwconv [options] -i <input image> -k <kernel image> -o <output image>";


Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> forcesum(string("--forcesum"), false,
		  string("forces convolution to be done by summation"),
		  false, no_argument);
Option<bool> forcefft(string("--forcefft"), false,
		  string("forces convolution to be done by fft"),
		  false, no_argument);
Option<bool> dilate(string("--dilate"), false,
		  string("perform gray-value dilation (max filter) using kernel as mask"),
		  false, no_argument);
Option<bool> erode(string("--erode"), false,
		  string("perform gray-value erosion (min filter) using kernel as mask"),
		  false, no_argument);
Option<bool> medianfilt(string("--median"), false,
		  string("perform median filtering using kernel as mask"),
		  false, no_argument);
Option<float> boxsize(string("-b,--box"), 0,
		  string("size of box kernel (in mm)"),
		  false, requires_argument);
Option<float> sphere(string("-s,--sphere"), 0,
		  string("radius of spherical kernel (in mm)"),
		  false, requires_argument);
Option<float> gaussian(string("-g,--gaussian"), 0,
		  string("size of Gaussian kernel (sigma in mm)"),
		  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input filename for base image"),
		  true, requires_argument);
Option<string> kername(string("-k,--kernel"), string(""),
		  string("input filename for kernel"),
		  false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename for convolved image"),
		  true, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Convolution Kernel Calculations

volume<float> box_kernel(float boxsize, const volume<float>& vin) {
  float b = boxsize;
  int sx = ((int) floor(b/vin.xdim()/2))*2 + 1;
  int sy = ((int) floor(b/vin.ydim()/2))*2 + 1;
  int sz = ((int) floor(b/vin.zdim()/2))*2 + 1;
  volume<float> vker(sx,sy,sz);
  vker.copyproperties(vin);
  vker = 1.0;
  return vker;
}


volume<float> spherical_kernel(float radius, const volume<float>& vin) {
  int sx = MISCMATHS::round(radius/vin.xdim())*2 + 1;
  int sy = MISCMATHS::round(radius/vin.ydim())*2 + 1;
  int sz = MISCMATHS::round(radius/vin.zdim())*2 + 1;
  volume<float> vker(sx,sy,sz);
  vker.copyproperties(vin);
  vker = 0.0;
  float dx2=Sqr(vin.xdim());
  float dy2=Sqr(vin.ydim());
  float dz2=Sqr(vin.zdim());
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


volume<float> gaussian_kernel(float sigma, const volume<float>& vin) {
  int sx = ((int) ceil(sigma*3/vin.xdim()))*2 + 1;
  int sy = ((int) ceil(sigma*3/vin.ydim()))*2 + 1;
  int sz = ((int) ceil(sigma*3/vin.zdim()))*2 + 1;
  volume<float> vker(sx,sy,sz);
  vker.copyproperties(vin);
  float dx2=Sqr(vin.xdim());
  float dy2=Sqr(vin.ydim());
  float dz2=Sqr(vin.zdim());
  for (int z=-sz/2; z<=sz/2; z++) {
    for (int y=-sy/2; y<=sy/2; y++) {
      for (int x=-sx/2; x<=sx/2; x++) {
	vker(x+sx/2,y+sy/2,z+sz/2)=exp(-(x*x*dx2+y*y*dy2+z*z*dz2)/(2*sigma*sigma));
      }
    }
  }
  return vker;
}

////////////////////////////////////////////////////////////////////////////

// Convolution code

// template <class T, class S>
// int insertpart(volume<T>& v1, const volume<S>& v2) 
// {
//   for (int z=v2.minz(); z<=v2.maxz(); z++) {
//     for (int y=v2.miny(); y<=v2.maxy(); y++) {
//       for (int x=v2.minx(); x<=v2.maxx(); x++) {
// 	v1(x,y,z)=(T) v2(x,y,z);
//       }
//     }
//   }
//   return 0;
// }


// template <class T, class S>
// volume<S> extractpart(const volume<T>& v1, const volume<S>& v2, 
// 		      const volume<S>& kernel) 
// {
//   volume<S> vout=v2;
//   vout = (S) 0.0;
//   int kxoff = (kernel.xsize()-1)/2;
//   int kyoff = (kernel.ysize()-1)/2;
//   int kzoff = (kernel.zsize()-1)/2;
//   for (int z=v2.minz(); z<=v2.maxz(); z++) {
//     for (int y=v2.miny(); y<=v2.maxy(); y++) {
//       for (int x=v2.minx(); x<=v2.maxx(); x++) {
// 	vout(x,y,z)=(S) v1(x+kxoff,y+kyoff,z+kzoff);
//       }
//     }
//   }
//   return vout;
// }


// float fsllog2(float x)
// {
//   // a cygwin annoyance!
//   return log(x)/log(2);
// }


volume<float> efficient_convolve(const volume<float>& vin,
				 const volume<float>& vker)
{
  bool usefft=true;
  if (forcesum.set()) { usefft=false; }
  else if (forcefft.set()) { usefft=true; }
  else {
    // estimate calculation time for the two methods and pick the best
    float offt = 2 * vin.nvoxels() * fsllog2(2 * vin.nvoxels());
    float osum = vin.nvoxels() * vker.nvoxels();
    float scalefactor = 44;  // relative unit operation cost for fft vs sum
    usefft = (osum > offt * scalefactor);
  }
  if (usefft) {
    int sx = Max(vin.xsize(),vker.xsize())*2;
    int sy = Max(vin.ysize(),vker.ysize())*2;
    int sz = Max(vin.zsize(),vker.zsize())*2;
    complexvolume vif, vkf;
    vif.re().reinitialize(sx,sy,sz);
    //vif.re().copyproperties(vin);
    vif.re() = 0.0;
    vif.im() = vif.re();
    vkf = vif;
    insertpart(vif.re(),vin);
    insertpart(vkf.re(),vker);
    fft3(vif);
    fft3(vkf);
    vif *= vkf;
    ifft3(vif);
    return extractpart(vif.re(),vin,vker);
  } else {
    // Direct intensity-based summation method
    return convolve(vin,vker);
  }
}

////////////////////////////////////////////////////////////////////////////

// Main execution stuff

int do_work(int argc, char* argv[]) 
{
  volume<float> vin, vout, vker;
  read_volume(vin,inname.value());

  if (boxsize.set()) { vker = box_kernel(boxsize.value(),vin); }

  if (sphere.set()) { vker = spherical_kernel(sphere.value(),vin); }

  if (gaussian.set()) { vker = gaussian_kernel(gaussian.value(),vin); }

  if (kername.set()) { read_volume(vker,kername.value()); }

  if (verbose.value()) { print_volume_info(vker,"kernel"); }

  if (dilate.value()) {
    vout = morphfilter(vin,vker,"dilate");
  } else if (erode.value()) {
    vout = morphfilter(vin,vker,"erode");
  } else if (medianfilt.value()) {
    vout = morphfilter(vin,vker,"median");
  } else {
    if (kername.unset()) { vker /= vker.sum(); }
    vout = efficient_convolve(vin,vker);
  }

  save_volume(vout,outname.value());
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(inname);
    options.add(outname);
    options.add(kername);
    options.add(gaussian);
    options.add(boxsize);
    options.add(sphere);
    options.add(dilate);
    options.add(erode);
    options.add(medianfilt);
    options.add(forcesum);
    options.add(forcefft);
    options.add(verbose);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ||
	 (kername.unset() && boxsize.unset() && sphere.unset() 
	  && gaussian.unset()) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions

  return do_work(argc,argv);
}

