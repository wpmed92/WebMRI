/*  surface_norm.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

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

// Calculates the surface normals for a mask, using a smoothed
//  gradient calculation (all non-surface points get zero ouput)

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="surface_norm (Version 1.0)\nCopyright(c) 2004, University of Oxford (Mark Jenkinson)";
string examples="surface_norm [options] -i <input mask image> -s <smoothing in mm> -o <output surface normal image>";

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
Option<float>  smoothing(string("-s"),3.0,
			 string("smoothing in mm (default is 3.0)"),
			 false, requires_argument);
Option<string> inname(string("-i"), string(""),
		      string("input mask filename"),
		      true, requires_argument);
Option<string> outname(string("-o"), string(""),
		       string("output surface normal filename"),
		       true, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions

int do_work(int argc, char* argv[]) 
{
  volume<float> vin;
  volume4D<float> snorm;
  read_volume(vin,inname.value());

  // force input to be a binary mask
  vin.binarise(0.5);
  if (verbose.value()) print_info(vin,"vin");

  // smooth the mask volume and calculate the gradient
  if (verbose.value()) { cerr << "Performing smoothing" << endl; }
  { 
    volume<float> smoothv;
    smoothv = smooth(vin,smoothing.value());
    if (verbose.value()) print_info(smoothv,"smoothv");
    if (verbose.value()) { cerr << "Calculating the gradient" << endl; }
    gradient(smoothv,snorm);
  }

  if (verbose.value()) print_info(snorm,"snorm");

  bool atedge;
  if (verbose.value()) { cerr << "Normalising gradients" << endl; }
  // zero all non-surface values and normalise vectors at the surface
  for (int z=vin.minz(); z<=vin.maxz(); z++) {
    for (int y=vin.miny(); y<=vin.maxy(); y++) {
      for (int x=vin.minx(); x<=vin.maxx(); x++) {
	atedge = false;
	if ( (vin(x,y,z)>0.5) ) {
	  if (vin(x,y,z-1)<0.5) atedge=true;
	  else { 
	    if (vin(x,y-1,z)<0.5) atedge=true;
	    else { 
	      if (vin(x-1,y,z)<0.5) atedge=true;
	      else {
		if (vin(x+1,y,z)<0.5) atedge=true;
		else {
		  if (vin(x,y+1,z)<0.5) atedge=true;
		  else {
		    if (vin(x,y,z+1)<0.5) atedge=true;
		  }
		}
	      }
	    }
	  }
	}
	if (atedge) {
	  float norm;
	  norm = 1.0 / sqrt(Sqr(snorm(x,y,z,0)) + Sqr(snorm(x,y,z,1)) 
			    + Sqr(snorm(x,y,z,2)) );
	  snorm(x,y,z,0) *= norm;
	  snorm(x,y,z,1) *= norm;
	  snorm(x,y,z,2) *= norm;
	} else {
	  snorm(x,y,z,0) = 0.0;
	  snorm(x,y,z,1) = 0.0;
	  snorm(x,y,z,2) = 0.0;
	}
      }
    }
    if (verbose.value()) { cerr << "."; }
  }
  if (verbose.value()) { cerr << endl; }

  // save the results
  save_volume4D(snorm,outname.value());

  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(smoothing);
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

