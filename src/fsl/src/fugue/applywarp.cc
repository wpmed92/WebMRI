/*  applywarp.cc

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

#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "newimage/warpfns.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

////////////////////////////////////////////////////////////////////////////

// COMMAND LINE OPTIONS

string title="applywarp (Version 1.2)\nCopyright(c) 2001, University of Oxford (Mark Jenkinson)";
string examples="applywarp -i invol -o outvol -r refvol -w warpvol";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> abswarp(string("--abs"), false,
		  string("treat warp field as absolute: x' = w(x)"),
		  false, no_argument);
Option<bool> relwarp(string("--rel"), false,
		  string("treat warp field as relative: x' = x + w(x)"),
		  false, no_argument);
Option<string> interp(string("--interp"), string(""),
		   string("interpolation method {nn,trilinear,sinc}"),
		   false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		       string("filename of input image (to be warped)"),
		       true, requires_argument);
Option<string> refname(string("-r,--ref"), string(""),
		       string("filename for reference image"),
		       false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("filename for output (warped) image"),
		       true, requires_argument);
Option<string> warpname(string("-w,--warp"), string(""),
			string("filename for warp transform (volume)"),
			true, requires_argument);
Option<string> maskname(string("-m,--mask"), string(""),
		       string("filename for mask image (in reference space)"),
		       false, requires_argument);
Option<string> prematname(string("--premat"), string(""),
		       string("filename for pre-transform (affine matrix)"),
		       false, requires_argument);
Option<string> postmatname(string("--postmat"), string(""),
		       string("filename for post-transform (affine matrix)"),
		       false, requires_argument);


bool abs_warp=false;

////////////////////////////////////////////////////////////////////////////

int applywarp()
{

  // read in pre/post transforms
  Matrix premat, postmat;
  premat = Identity(4);
  postmat = Identity(4);
  if (prematname.set()) {
    read_ascii_matrix(premat,prematname.value());
  }
  if (postmatname.set()) {
    read_ascii_matrix(postmat,postmatname.value());
  }
  
  // read in images
  volume4D<float> invol, outvol;
  volume<float> refvol, mask;
  volumeinfo vinfo;
  volume4D<float> warpvol;

  read_volume4D(invol,inname.value(),vinfo);
  read_volume4D(warpvol,warpname.value());
  if (warpvol.tsize()<3) {
    cerr << "Warp volume does not have the correct dimensions (min. 3 volumes)"
	 << endl;
    exit(1);
  }

  if (abswarp.set()) { abs_warp = true; }
  else if (relwarp.set()) { abs_warp = false; }
  else {
    if (verbose.value()) { 
      cout << "Automatically determining relative/absolute warp conventions" 
	   << endl; 
    }
    // try to determine this automatically
    float stddev0 = warpvol[0].stddev()+warpvol[1].stddev()+warpvol[2].stddev();
    convertwarp_abs2rel(warpvol);
    float stddev1 = warpvol[0].stddev()+warpvol[1].stddev()+warpvol[2].stddev();
    // restore to the original form
    convertwarp_rel2abs(warpvol);
    // assume that relative warp always has less stddev
    if (stddev0>stddev1) {
      // the initial one (greater stddev) was absolute
      if (verbose.value()) {cout << "Warp convention = absolute" << endl;}
      abs_warp=true;
    } else {
      // the initial one was relative
      if (verbose.value()) {cout << "Warp convention = relative" << endl;}
      abs_warp=false;
    }
  }

  if (!abs_warp) {
    // warpvol needs to be in absolute convention from here on : x' = w(x)
    convertwarp_rel2abs(warpvol);
  }


  if (maskname.set()) { read_volume(mask,maskname.value()); }

  if (refname.set()) {
    read_volume(refvol,refname.value());
  } else {
    refvol = warpvol[0];
  }
  outvol = invol;

  volume<float> tmpvol;
  for (int t=0; t<invol.tsize(); t++) {
    invol[t].setpadvalue(invol[t].backgroundval());
    invol[t].setextrapolationmethod(extraslice);
    
    // set interpolation method
    if (interp.value() == "nn" ) {
      invol.setinterpolationmethod(nearestneighbour);
    } else if (interp.value() == "trilinear") {
      invol.setinterpolationmethod(trilinear);
    } else if (interp.value() == "sinc") {
      invol.setinterpolationmethod(sinc);
    }
    
    // do the deed
    tmpvol = refvol;
    apply_warp(invol[t],tmpvol,warpvol,premat,postmat);
    if (maskname.set()) {
      outvol[t] = tmpvol * mask;
    } else {
      outvol[t] = tmpvol;
    }
  }

  // save the results
  save_volume4D_dtype(outvol,outname.value(),dtype(inname.value()),vinfo);
  
  return 0;
}




int main(int argc, char *argv[])
{

  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(inname);
    options.add(refname);
    options.add(warpname);
    options.add(abswarp);
    options.add(relwarp);
    options.add(outname);
    options.add(prematname);
    options.add(postmatname);
    options.add(maskname);
    options.add(interp);
    options.add(verbose);
    options.add(help);
    
    options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
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

  return applywarp();
}

