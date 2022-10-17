/*  avcheck.cc - Utility to check validity of avw header file 
    
    Peter Bannister, FMRIB Image Analysis Group
    
    Copyright (C) 1999-2001 University of Oxford  */

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


#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <sstream>

#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "utils/options.h"
#include "newimage/fmribmain.h"

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace Utilities;

string title="avwcheck (Version 1.0)\nChecks Header file info for voxel dimensions\nCopyright(c) 2001, University of Oxford (Peter R Bannister)";
string examples="avwcheck -i <header_file>\n";

Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<string> inputname(string("-i,--in"), string(""),
			 string("filename of input timeseries"),
			 true, requires_argument);
const float tolerance = 0.00000001;

template <class T>
int fmrib_main(int argc, char* argv[])
{
  volume4D<T> timeseries;
  volumeinfo vinfo;
  float x_dim, y_dim, z_dim, t_dim;
  short datatype;

  if (inputname.set()) {
    if (verbose. value()) { cout << "Reading header file" << endl; }
    read_volume4D_hdr_only(timeseries,inputname. value(),vinfo);
  }

  if (fabs(timeseries. xdim()) <= tolerance) {
    do {
      cerr << "Invalid voxel x-dimension (" << timeseries. xdim() << " mm): Please enter a new value" << endl;
      cin >> x_dim;
    } while (fabs(x_dim) <= tolerance);
    timeseries. setxdim(x_dim);
  }

  if (fabs(timeseries. ydim()) <= tolerance) {
    do {
      cerr << "Invalid voxel y-dimension (" << timeseries. ydim() << " mm): Please enter a new value" << endl;
      cin >> y_dim;
    } while (fabs(y_dim) <= tolerance);
    timeseries. setydim(y_dim);
  }
  
  if (fabs(timeseries. zdim()) <= tolerance) {
    do {
      cerr << "Invalid voxel z-dimension (" << timeseries. zdim() << " mm): Please enter a new value" << endl;
      cin >> z_dim;
    } while (fabs(z_dim) <= tolerance);
    timeseries. setzdim(z_dim);
  }
  
  if (fabs(timeseries. tdim()) <= tolerance) {
    do {
      cerr << "Invalid TR (" << timeseries. tdim() << " secs): Please enter a new value" << endl;
      cin >> t_dim;
    } while (fabs(t_dim) <= tolerance);
    timeseries. settdim(t_dim);
  }
  
  // set bitpix correctly (via FslSetDataType) - NB: type T = datatype
  FslGetDataType(&vinfo,&datatype);
  FslSetDataType(&vinfo,datatype);
  // read in whole file(!) and save it again
  volume4D<T> tmp;
  read_volume4D(tmp, inputname.value());
  save_volume4D_filetype(tmp, inputname.value(), FslGetFileType(&vinfo), vinfo);
  return 0;
}

int main (int argc,char** argv)
{
  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(inputname);
    options.add(help);
    options.add(verbose);

    options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
        
    if ( inputname.unset()) 
      {
	options.usage();
	cerr << endl 
	     << "--in or -i MUST be used." 
	     << endl;
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  int retval=call_fmrib_main(dtype(inputname.value()), argc, argv);
  if (retval!=0) {
    cerr << "Failed to correctly read file, please check the .hdr using avwhd" << endl;
  } else return retval;
}

