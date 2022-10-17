/*  avwmeants.cc

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


// Creates a mean time series (ignoring zeros) from the input 4D volume
// Saves the results as a column in a text file

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;


using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="avwmeants (Version 1.1)\nCopyright(c) 2004, University of Oxford (Mark Jenkinson)\nPrints average timeseries (intensities) to the screen (or saves to a file).\nThe average is taken over all voxels in the mask (or all voxels in the image if no mask is specified).\n";
string examples="avwmeants -i filtered_func_data -o meants.txt -m my_mask\navwmeants -i filtered_func_data -m my_mask\navwmeants -i filtered_func_data -c 24 19 10";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> usemm(string("--usemm"), false,
		  string("\tuse mm instead of voxel coordinates (for -c option)"),
		  false, no_argument);
Option<bool> showall(string("--showall"), false,
		  string("show all voxel time series (within mask) instead of averaging"),
		  false, no_argument);
Option<string> inname(string("-i"), string(""),
		  string("~<filename>\tinput 4D image"),
		  true, requires_argument);
Option<string> maskname(string("-m"), string(""),
		  string("~<filename>\tinput 3D mask"),
		  false, requires_argument);
Option<string> outmat(string("-o"), string(""),
		  string("~<filename>\toutput text matrix"),
		  false, requires_argument);
Option<float> coordval(string("-c"), 0.0,
		       string("~<x y z>\trequested spatial coordinate (instead of mask)"), 
		       false, requires_3_arguments);

int nonoptarg;



int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  options.add(inname);
  options.add(outmat);
  options.add(maskname);
  options.add(coordval);
  options.add(usemm);
  options.add(showall);
  options.add(verbose);
  options.add(help);
  
  nonoptarg = options.parse_command_line(argc, argv);
  
  // line below stops the program if the help was requested or 
  //  a compulsory option was not set
  if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
    {
      options.usage();
      exit(EXIT_FAILURE);
    }
  

  // OK, now do the job ...

  volume4D<float> vin;
  read_volume4D(vin,inname.value());

  volume<float> mask;
  if (maskname.set()) {
    read_volume(mask,maskname.value());
  } else {
    mask = vin[0];
    mask = 1.0;
  }


  // override the mask if a coordinate is specified
  if (coordval.set()) {
    mask = vin[0];
    mask = 0.0;
    float x, y, z;
    x = coordval.value(0);
    y = coordval.value(1);
    z = coordval.value(2);
    if (usemm.value()) {
      // convert from mm to voxels
      x /= vin[0].xdim();
      y /= vin[0].ydim();
      z /= vin[0].zdim();
    }
    mask(MISCMATHS::round(x),MISCMATHS::round(y),MISCMATHS::round(z)) = 1.0;
  }

  if (!samesize(vin[0],mask)) {
    cerr << "ERROR: Mask and Input volumes have different (x,y,z) size." 
	 << endl;
    return -1;
  }

  mask.binarise(1e-8);  // arbitrary "0" threshold

  int nt = vin.tsize();
  int nvox = 1;
  if (showall.value()) {
    nvox = MISCMATHS::round(mask.sum()+0.5);
    nt += 3;
  }
  Matrix meants(nt,nvox);
  meants = 0;
  long int num=0;

  for (int z=mask.minz(); z<=mask.maxz(); z++) {
    for (int y=mask.miny(); y<=mask.maxy(); y++) {
      for (int x=mask.minx(); x<=mask.maxx(); x++) {
	if (mask(x,y,z)>0.5) {
	  num++;
	  if (showall.value()) {
	    meants(1,num) = x;
	    meants(2,num) = y;
	    meants(3,num) = z;
	    meants.SubMatrix(4,nt,num,num) = vin.voxelts(x,y,z);
	  } else {
	    meants += vin.voxelts(x,y,z);
	  }
	}
      }
    }
  }

  if (verbose.value()) {
    cout << "Number of voxels used = " << num << endl;
  }

  // normalise for number of valid entries if averaging
  if (!showall.value()) {
    if (num>0) meants /= (float) num;
  }

  // save the result
  if (outmat.set()) {
    write_ascii_matrix(meants,outmat.value());
  } else {
    cout << meants << endl;
  }

  return 0;
}

