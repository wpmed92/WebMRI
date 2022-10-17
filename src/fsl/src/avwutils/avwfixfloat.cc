/*  avwfixfloat.cc

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

#include <math.h>
#include <iostream>
#include <string>
#include "newimage/newimageall.h"

using namespace NEWIMAGE;

void dump_float(float *fptr) {
  char *cptr;
  cptr = (char *) fptr;
  cerr << (int) *cptr << " , " << (int) *(cptr+1) << " , " 
       << (int) *(cptr+2) << " , " << (int) *(cptr+3) << " ";
}


//  void check_volume(const volume4D<float>& vol1) {
//    long int count=0;
//    float val, sum=0.0;
//    for (int t=0; t<vol1.tsize(); t++) {
//      for (int z=0; z<vol1.zsize(); z++) {
//        for (int y=0; y<vol1.ysize(); y++) {
//  	for (int x=0; x<vol1.xsize(); x++) {
//  	  count++;
//  	  val = vol1(x,y,z,t);
//  	  sum += val;
//  	}
//        }
//      }
//    }
//  }


int main(int argc, char *argv[]) 
{
 
  if (argc<2) { 
    cerr << "Usage: " << argv[0] << " <input_volume> [output_volume] [-v]" 
	 << endl; 
    return -1;
  }
  
  string outname="";
  bool verbose = false;
  if (argc>=3) {
    outname = argv[2];
    if (outname == "-v" ) {
      verbose = true;
      outname = "";
    }
  }

  string optarg;
  if (argc>=4) {
    optarg = argv[3];
    if (optarg=="-v") {
      verbose = true;
    }
  }

  volume4D<float> vol1;
  volumeinfo vinfo;
  read_volume4D(vol1,argv[1],vinfo);
  long int count=0, badcount=0, nancount=0;
  float val, sum=0;
  for (int t=0; t<vol1.tsize(); t++) {
    for (int z=0; z<vol1.zsize(); z++) {
      for (int y=0; y<vol1.ysize(); y++) {
	for (int x=0; x<vol1.xsize(); x++) {
	  count++;
	  if (verbose) {
	    cerr << "x,y,z,t = " << x <<","<<y<<","<<z<<","<<t<<" : ";
	    dump_float(&vol1(x,y,z,t));
	  }
	  char *cptr;
	  cptr = (char *) &(vol1(x,y,z,t));
	  int ival0 = (int) *cptr;
	  int ival1 = (int) *(cptr+1);
	  int ival2 = (int) *(cptr+2);
	  int ival3 = (int) *(cptr+3);

	  if ( ( (ival3==0) && ( (ival0!=0) || (ival1!=0) || (ival2!=0) ) ) ||
	       (ival3<-125)  )
	    {
	      badcount++;
	      if (verbose) {
		cerr << "BAD VALUE DETECTED - fixing it" << endl;
		dump_float((float *) cptr);  cerr << endl;
	      }
	      *(cptr+0) = 0;
	      *(cptr+1) = 0;
	      *(cptr+2) = 0;
	      *(cptr+3) = 0;
	    }

	  if (! finite(vol1(x,y,z,t))) {
	    nancount++;
	    vol1(x,y,z,t) = 0.0;
	  }

	  val = vol1(x,y,z,t);
	  sum += val;
	  if (verbose) cerr << "  ::  " << vol1(x,y,z,t) << "  ::  " << endl;
	}
	if (verbose) cerr << y << ",";
      }
      if (verbose) cerr << endl << "Z = " << z << endl;
    }
    if (verbose) cerr << endl << "T = " << t << endl;
  }
  cout << "Successfully parsed volume" << endl;
  cout << "Found " << badcount << " invalid elements out of " << count  << endl;
  cout << "Found " << nancount << " non-finite elements out of " << count  
       << endl;
  if (verbose) print_volume_info(vol1,"Volume");
  if (outname.size()>0) {
    cout << "Saving volume ... " << endl;
    save_volume4D(vol1,outname,vinfo);
    cout << "done" << endl;
  }
}
