/*  midtrans.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2003 University of Oxford  */

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

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"
#include <vector>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace std;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="midtrans (Version 1.0)\nCopyright(c) 2004, University of Oxford (Mark Jenkinson)";
string examples="midtrans [options] transform1 transform2 ... transformN";

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
Option<float> outlier(string("-r,--reject"), 3.0,
		  string("threshold for outlier rejection (in stddev): default=3.0"),
		  false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename for matrix"),
		  false, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions


//DUMMY TEST STUFF - MUST REMOVE!!!
typedef Matrix complexmatrix;
complexmatrix logm(const Matrix& mat) { return mat; }
Matrix expm(const complexmatrix& mat) { return mat; }
//END OF DUMMY TEST STUFF (MJ)

int do_work(int argc, char* argv[]) 
{
  int ntrans = argc - nonoptarg;

  // read in matrices
  Matrix midtrans(4,4);
  vector<Matrix> matlist;
  for (int n=nonoptarg; n<argc; n++) {
    Matrix trans(4,4);
    trans = read_ascii_matrix(string(argv[n]));
    if (fabs(trans.Determinant()) < 1e-6) {
      cerr << "Could not read matrix " << argv[n] << endl;
      exit(EXIT_FAILURE);
    } 
    matlist.push_back(trans);
  }

  // find the geometric mean via logm and expm
  complexmatrix logsum(4,4);
  logsum=0.0;
  for (int n=0; n<ntrans; n++) {
    logsum += logm(matlist[n]);
  }
  midtrans = expm(logsum / (double) ntrans);

  // perform outlier reject (if necessary and possible)
  if (outlier.set() && (ntrans < 5)) {
     cerr << "Warning: outlier reject is not done when there are less than 5 transforms"
	  << endl;
  }
  if ((ntrans>=5) && (outlier.value()>0)) {
    /* use rmsdiff to measure differences from midtrans with and
	without including the 2 max values (assume no more than
	two outliers!)  
       should really correct for standard deviation estimation from
	ordered statistics (excluding top 2 values) but probably 
	doesn't make much difference - maybe later */
    float sum=0.0, sum2=0.0, max1=0.0, max2=0.0;
    float rmax=80;  // maximum radius in mm (for rmsdiff)
    ColumnVector centre(3);
    centre=0.0;
    int n1=-1, n2=-1;
    for (int n=0; n<ntrans; n++) {
	float rmsd=rms_deviation(matlist[n],midtrans,centre,rmax);
	if (rmsd>max2) {
	  if (rmsd>max1) { max2=max1; max1=rmsd; n2=n1; n1=n; }
	  else { max2=rmsd; n2=n; }
        }
	sum+=rmsd;
	sum2+=Sqr(rmsd);
    }
    // adjust sums to remove the effect of max1 and max2
    sum -= (max1 + max2);
    sum2 -= (Sqr(max1) + Sqr(max2));
    int n = ntrans-2;
    float stddev = (sum2 - Sqr(sum)/n)/(n-1);

    // remove max 1 or 2 if they are over outlier rejection threshold
    bool rm1=false, rm2=false;
    rm1 = (max1/stddev) > outlier.value();
    rm2 = (max2/stddev) > outlier.value();
    if (rm1 || rm2) {
	// now re-estimate the mid-trans
        logsum=0.0;
        for (int n=0; n<ntrans; n++) {
          if ( (rm2 && (n==n2)) || (rm1 && (n==n1)) ) {
	    // do nothing 
	  } else {
	    logsum += logm(matlist[n]);
          }
        }
        int nval=ntrans;
	if (rm1) nval--;
	if (rm2) nval--;
        midtrans = expm(logsum / (double) nval);
    }
  }
  

  // show/save output matrix
  if (outname.set()) {
    write_ascii_matrix(midtrans,outname.value());
  }
  if (verbose.value() || outname.unset()) {
    cout << midtrans << endl;
  }

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
    options.add(outname);
    options.add(outlier);
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
    
    if ((argc - nonoptarg)<2) {
      cerr << "Must specify at least 2 transforms to find the mid-transform"
  	<< endl;
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

