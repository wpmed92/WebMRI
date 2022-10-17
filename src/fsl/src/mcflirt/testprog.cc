/*  test.cc Evaluates RMSdiff between matricies
    
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


#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>

#include "newmatap.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"

#include "Globaloptions.h"
#include "Log.h"

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace UTILS;
using namespace MISCMATHS;

int main (int argc,char** argv){
  
  volume4D<float> timeseries;
  volume<float> vout;
  volumeinfo vinfo;
  Matrix displacement[200], correction[200];
  ColumnVector cvec(3);
  cvec = 0.0;
  int no_volumes;
  float rmax=80.0;
  string inputfname;
  string inputbasename = "/usr/people/prb/medx/motion/";
  string directory;
  ColumnVector vtemp(6);
  ColumnVector centre(3);
  centre = 0;
  Log& logger = Log::getInstance();

  if (argc < 3){
    cerr << "Usage: testprog <reference_mats_dir> <test_mats_dir>" << endl;
    exit(1);
  }

  inputfname = "/usr/people/prb/medx/motion/validation/fmri_3.2/fmri_3.2.hdr";
  
  cerr << "Reading stationary time series" << inputfname << " ... " << endl;
  read_volume4D(timeseries,inputfname,vinfo);
  no_volumes = timeseries.tsize();

  cerr << no_volumes << endl;
  
  vout = timeseries[0];
  
  centre(1) = 0.5*(vout.xsize() - 1.0)*vout.xdim();
  centre(2) = 0.5*(vout.ysize() - 1.0)*vout.ydim();
  centre(3) = 0.5*(vout.zsize() - 1.0)*vout.zdim();

  directory = argv[1];
  logger.setDir(directory);

  cerr << directory << endl;

  for (int i = 0; i < no_volumes; i++) {
    
    ostringstream osc;
    osc << "MAT_"<< setw(4) << setfill('0') << i;
    cerr << "Reading matrix: " << osc.str() << endl;
    logger.in(osc.str().c_str(), displacement[i],4,4);
    cerr << "read in matrix" << endl;
    cerr << displacement[i] << endl;
  }

  /*
   logger.setDir("/usr/people/prb/medx/motion/roberto/group2/subj9/mats");
   for (int i = 0; i < no_volumes; i++) {
    
     char strc[50];
     ostringstream osc(strc,50);
     if (i < 10)
       osc << "MAT_000" << i << '\0';
     else if (i < 100)
       osc << "MAT_00" << i << '\0';
     else if (i < 1000)
       osc << "MAT_0" << i << '\0';
     cerr << "Reading matrix: " << strc << endl;
     logger.in(strc, displacement[i],4,4);
     cerr << "read in matrix" << endl;
     cerr << displacement[i] << endl;
   }
  */

  /*
   logger.setDir("/usr/people/prb/medx/motion/sinctest/spm/");
  
   for (int i = 0; i < no_volumes; i++) {
     char strc[50];
     ostringstream osc(strc,50);
     osc << "spm_mat_corrected_"<< i << '\0';
     cerr << "Reading matrix: " << strc << endl;
     logger.in(strc, correction[i],4,4);
     cerr << "read in matrix" << endl;
     cerr << correction[i] << endl;
   } 
  */

  /*
  logger.setDir("/usr/people/prb/medx/motion/validation/reference/air/");
  
  
  for (int i = 0; i < no_volumes; i++) {
     char strc[50];
     ostringstream osc(strc,50);
     osc << "newairmat_"<< i+1  << '\0';
     cerr << "Reading matrix: " << strc << endl;
     logger.in(strc, correction[i], 4, 4);
     cerr << "read in matrix" << endl;
     cerr << correction[i] << endl;
  } 
  */
 
 
  directory = argv[2];
  logger.setDir(directory);
  
  /*
  for (int i = 0; i < no_volumes; i++) {
    char strc[50];
    ostringstream osc(strc,50);
    osc << "spm_mat_corrected_"<< setw(4) << setfill('0') << i << '\0';
    cerr << "Reading matrix: " << strc << endl;
    logger.in(strc, correction[i], 4, 4);
    cerr << "read in matrix" << endl;
    cerr << correction[i] << endl;
  }
  */
  
  /*
  for (int i = 0; i < 10; i++) {
    char strc[50];
    ostringstream osc(strc,50);
    osc << "mat_000"<< i << '\0';
    cerr << "Reading matrix: " << strc << endl;
    logger.in(strc, correction[i], 4, 4);
    cerr << "read in matrix" << endl;
    cerr << correction[i] << endl;
  }

  for (int i = 10; i < no_volumes; i++) {
    char strc[50];
    ostringstream osc(strc,50);
    osc << "mat_00"<< i << '\0';
    cerr << "Reading matrix: " << strc << endl;
    logger.in(strc, correction[i], 4, 4);
    cerr << "read in matrix" << endl;
    cerr << correction[i] << endl;
  }
  */

   for (int i = 0; i < no_volumes; i++) {    
     char strc[50];
     ostringstream osc(strc,50);
     osc << "MAT_"<< setw(4) << setfill('0') << i << '\0';
     cerr << "Reading matrix: " << strc << endl;
     logger.in(strc, correction[i], 4, 4);
     cerr << "read in matrix" << endl;
     cerr << correction[i] << endl;
   }


    for (int i = 0; i < no_volumes; i++) {
      //float rms = rms_deviation(Identity(4),displacement[no_volumes/2]*displacement[i].i(),centre,rmax);// ref
      float rms = rms_deviation(correction[i],displacement[no_volumes/2]*displacement[i].i(),centre,rmax);// middle ref 
      //float rms = rms_deviation(correction[i], displacement[0]*displacement[i].i(),centre,rmax); // first ref
      cout << rms << endl;
  } 
 
}























