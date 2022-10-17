/*  vmapper.cc
    
    Peter Bannister, FMRIB Image Analysis Group
    
    Copyright (C) 1999-2000 University of Oxford  */

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

#include "newimage/newimageall.h"
#include "newmatap.h"

#include "Globaloptions.h"

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace NEWIMAGE;

//------------------------------------------------------------------------//


int main (int argc,char** argv)
{
  volume4D<float> timeseries;
  volume<float> refvol, anisorefvol, variancevol, meanvol, sigmavol;
  int vmax;
  Globaloptions& globalopts = Globaloptions::getInstance();
  

  cerr << endl << "vmapper v 0.1.1" << endl <<"This software is currently alpha." << endl << "Please send any comments/ queries to: prb@fmrib.ox.ac.uk" << endl << endl;
  globalopts.parse_command_line(argc, argv);
  cerr << "Reading time series... " << endl;
  volumeinfo vinfo;
  read_volume4D(timeseries, globalopts.inputfname, vinfo);
  vmax = timeseries.tsize();
  
  meanvol = timeseries[0];
  variancevol = timeseries[0];
  sigmavol = timeseries[0];
  
  meanvol = 0.0;
  variancevol = 0.0;
  sigmavol = 0.0;

  cerr << "This version has been hacked for AIR/ SPM corrected data." << endl  
       << "The end slices are ignored when calculating mean and variance statistics" << endl;

  // calculate the mean value and variance at each voxel
  for (int x=0; x< timeseries[0].xsize(); x++) {
    for (int y=0; y< timeseries[0].ysize(); y++) {
      for (int z=1; z< (timeseries[0].zsize()-1); z++) {
	for (int i=0; i< vmax; i++) 
	  meanvol(x,y,z) += timeseries[i](x,y,z);
	meanvol(x,y,z) = meanvol(x,y,z)/(float)vmax;
	}
    }
  }

  cerr << "generated mean image" << endl;

  // change limits on z index for end slice exclusion
  for (int x=0; x< timeseries[0].xsize(); x++) {
    for (int y=0; y< timeseries[0].ysize(); y++) {
      for (int z=1; z< (timeseries[0].zsize()-1); z++) {
	for (int i=0; i< vmax; i++)
	  variancevol(x,y,z) += (timeseries[i](x,y,z) - meanvol(x,y,z))*(timeseries[i](x,y,z) - meanvol(x,y,z));
	variancevol(x,y,z) = variancevol(x,y,z)/((float)(vmax - 1));
	sigmavol(x,y,z) = sqrt(variancevol(x,y,z));
      }
    }
  }
  
  cerr << "Saving mean volume... " << endl;
  save_volume_dtype(meanvol, "mean_"+globalopts.inputfname, dtype(globalopts.inputfname), vinfo);

  cerr << "Saving variance volume... " << endl;
  save_volume_dtype(variancevol, "variance_"+globalopts.inputfname, dtype(globalopts.inputfname), vinfo);

  cerr << "Saving standard deviation volume... " << endl;
  save_volume_dtype(sigmavol, "sigma_"+globalopts.inputfname, dtype(globalopts.inputfname), vinfo);
}








