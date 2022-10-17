/*  cylsamp.cc

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

// Performs cylindrical sampling over a surface.  That is, averaging
//  all values within a cylindrical region centred at each surface
//  point on a mask.


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

string title="cylsamp (Version 1.0)\nCopyright(c) 2004, University of Oxford (Mark Jenkinson)";
string examples="cylsamp [options] -i <input mask image> -s <smoothing in mm> -o <output surface normal image>";

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
Option<float>  radius(string("-r"),5.0,
		      string("radius of cylinder in mm (default is 5.0)"),
		      false, requires_argument);
Option<float>  height(string("-h"),10.0,
		      string("half-height of cylinder in mm (default is 10.0)"),
		      false, requires_argument);
Option<string> normname(string("-n"), string(""),
		      string("input surface normal image filename"),
		      true, requires_argument);
Option<string> edgemaskname(string("-e"), string(""),
		      string("input edge mask image filename"),
		      true, requires_argument);
Option<string> flowname(string("-f"), string(""),
		      string("input flow image filename"),
		      true, requires_argument);
Option<string> outname(string("-o"), string(""),
		       string("output surface normal filename"),
		       true, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions

void zeronans(volume<float>& vol) 
{
  for (int z=vol.minz(); z<=vol.maxz(); z++) {
    for (int y=vol.miny(); y<=vol.maxy(); y++) {
      for (int x=vol.minx(); x<=vol.maxx(); x++) {
	if (isnan(vol.value(x,y,z))) {
	  vol(x,y,z)=0;
	}
      }
    }
  }
}


void zeronans(volume4D<float>& vol)
{
  for (int t=vol.mint(); t<=vol.maxt(); t++) {
    zeronans(vol[t]);
  }
}


int do_work(int argc, char* argv[]) 
{
  volume<float> vflow, vedgemask, vsamp;
  volume4D<float> snorm;
  read_volume(vflow,flowname.value());
  read_volume(vedgemask,edgemaskname.value());
  read_volume4D(snorm,normname.value());
  zeronans(vflow);
  zeronans(vedgemask);
  zeronans(snorm);

  // set up output image
  if (verbose.value()) print_info(vflow,"vflow");
  vsamp = vflow;

  float r=radius.value();
  float r2 = r*r;
  float h=height.value();
  int len=(int) (sqrt(r*r + h*h)) + 1;
  int lenx = (int) ceil(len/vsamp.xdim());
  int leny = (int) ceil(len/vsamp.ydim());
  int lenz = (int) ceil(len/vsamp.zdim());

  if (verbose.value()) { cerr << "Performing Cylindrical Sampling" << endl; }
  for (int z=vsamp.minz(); z<=vsamp.maxz(); z++) {
    for (int y=vsamp.miny(); y<=vsamp.maxy(); y++) {
      for (int x=vsamp.minx(); x<=vsamp.maxx(); x++) {
	// check to see if it is an edge point (where snorm != 0)
	float vx, vy, vz;
	vx = snorm(x,y,z,0);
	vy = snorm(x,y,z,1);
	vz = snorm(x,y,z,2);
	if ( (fabs(vx)>1e-5) || (fabs(vy)>1e-5) || (fabs(vz)>1e-5) ) {
	  // OK, now do the cylindrical sampling
	  float tot=0.0;
	  int num=0;
	  for (int z1=Max(z-lenz,vsamp.minz()); z1<=Min(z+lenz,vsamp.maxz()); z1++) {
	    for (int y1=Max(y-leny,vsamp.miny()); y1<=Min(y+leny,vsamp.maxy()); y1++) {
	      for (int x1=Max(x-lenx,vsamp.minx()); x1<=Min(x+lenx,vsamp.maxx()); x1++) {
		if (vedgemask(x1,y1,z1)>0.5) {
		  // yy is the projection of vector onto major axis of cylinder
		  float yy = vx * (x1 - x) + vy * (y1 - y) + vz * (z1 - z);
		  if (fabs(yy)<=h) {
		    // xx2 is the square radius of the planar component of the vector
		    float xx2 = ( Sqr(x1-x) + Sqr(y1-y) + Sqr(z1-z) ) - Sqr(yy);
		    if (xx2 <= r2) {
		      // inside cylinder
		      num++;
		      tot += vflow(x1,y1,z1);
		    }
		  }
		}
	      }
	    }
	  }
	  vsamp(x,y,z) = tot / Max(num,1);  // prevent divide by zero
	} else {
	  vsamp(x,y,z) = 0.0;
	}
      }
    }
    if (verbose.value()) { cerr << "."; }
  }
  if (verbose.value()) { cerr << endl; }

  // save the results
  save_volume(vsamp,outname.value());

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
    options.add(flowname);
    options.add(edgemaskname);
    options.add(normname);
    options.add(outname);
    options.add(radius);
    options.add(height);
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

