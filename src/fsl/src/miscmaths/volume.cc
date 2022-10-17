/*  volume.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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

#include <iostream>
#include <cstdlib>
#include "fslio/fslio.h"

#include "newmatap.h"
#include "newmatio.h"
#include "volume.h"
#include "utils/time_tracer.h"

using namespace NEWMAT;
using namespace Utilities;

namespace MISCMATHS {
  
  void Volume::unthreshold()
     {  
       Time_Tracer ts("Volume::unthreshold");
       
       int sizeVolUnthresholded = getUnthresholdSize();
       int sizeVolThresholded = getVolumeSize();

       // Increase this to required size and set to 0
       Release(); 
       ColumnVector X=*this; 
       ReSize(sizeVolUnthresholded);
       ColumnVector::operator=(0);
       
       // Loop through volume restoring non thresholded values:
       for(int i = 1; i <= sizeVolThresholded; i++)
	 {
	   (*this)(int(preThresholdPositions(i)))=X(i);
	 }
     }
  
  void Volume::unthreshold(const VolumeInfo& pvolinfo,const ColumnVector& in)
  {
    volinfo = pvolinfo;
    preThresholdPositions = in;
    unthreshold();
  }
  
  void Volume::threshold()
  {  
    Time_Tracer ts("Volume::threshold");
    
    int sizeVol = preThresholdPositions.Nrows();
    
    ColumnVector X(sizeVol); 
    
    // Loop through pulling out non thresholded values
    for(int i = 1; i <= sizeVol; i++)
      {
	X(i)=(*this)(int(preThresholdPositions(i)));
      }
    
    (*this) = X;
  }
  
  void Volume::threshold(float thresh)
  {
    Time_Tracer ts("Volume::threshold");
    
    int size = getVolumeSize();
    int j = 0;
    
    preThresholdPositions.ReSize(size);
       
    float m = 0;
    for(int i = 1; i <= size; i++)
      {
	m = (*this)(i);
	   
	if(m > thresh)
	  {
	    j++;	
	    preThresholdPositions(j) = i;
	     
	    (*this)(j) = (*this)(i);
	  }
      }

    // Discard rest:
    *this = Rows(1,j);

    preThresholdPositions = preThresholdPositions.Rows(1,j);
  }



  void Volume::writeAsFloat(const string& fname)
  {    
    Time_Tracer ts(string("Volume::writeAsFloat" + fname).c_str());
       
    const ColumnVector& outputvol = *this;
    FSLIO* OP = FslOpen(fname.c_str(), "wb");

    FslCloneHeader(OP,volinfo.miscinfo);
   
    FslSetDim(OP,volinfo.x, volinfo.y, volinfo.z, 1);
    FslSetVoxDim(OP,volinfo.vx, volinfo.vy, volinfo.vz, 0);
    FslSetDataType(OP, DT_FLOAT);
    FslSetIntent(OP, volinfo.intent_code, volinfo.intent_p1, volinfo.intent_p2, 
		 volinfo.intent_p3);

    int sizeVol = outputvol.Nrows();

    float *qv = new float[sizeVol];

    for (int i=1; i<=sizeVol; i++) 
      { 
	qv[i-1] = outputvol(i);
      }

    FslWriteHeader(OP);
    FslWriteVolumes(OP, qv, 1);

    delete [] qv;

    FslClose(OP);
  }

  void Volume::writeAsInt(const string& fname)
  {    
    Time_Tracer ts("Volume::writeAsInt");

    const ColumnVector& outputvol = *this;
    FSLIO* OP = FslOpen(fname.c_str(), "wb");
 
    FslCloneHeader(OP,volinfo.miscinfo);
   
    FslSetDim(OP, volinfo.x, volinfo.y, volinfo.z, 1);  
    FslSetVoxDim(OP, volinfo.vx, volinfo.vy, volinfo.vz, 0);  
    FslSetDataType(OP, DT_SIGNED_SHORT);
    FslSetIntent(OP, volinfo.intent_code, volinfo.intent_p1, volinfo.intent_p2, 
		 volinfo.intent_p3);

    int sizeVol = outputvol.Nrows();

    short *qv = new short[sizeVol];
  

    for (int i=1; i<=sizeVol; i++) 
      { 
	qv[i-1] = (short)outputvol(i);
      }

    FslWriteHeader(OP);
    FslWriteVolumes(OP, qv, 1);

    delete [] qv;

    FslClose(OP);
  }
 
  void Volume::read(const string& fname)
  {
    Time_Tracer ts(string("Volume::read" + fname).c_str());
      
    FSLIO* IP = FslOpen(fname.c_str(), "rb");
    ColumnVector& output = *this;
       
    short x,y,z,v,type;
    float vx,vy,vz,tr;
    float slope, intercept;
    int doscaling;


    FslGetDim(IP,&x,&y,&z,&v);
    FslGetVoxDim(IP,&vx,&vy,&vz,&tr);
    FslGetIntent(IP, &(volinfo.intent_code), &(volinfo.intent_p1), &(volinfo.intent_p2), 
		 &(volinfo.intent_p3));
    doscaling = FslGetIntensityScaling(IP,&slope,&intercept);

    volinfo.x = x; volinfo.y = y; volinfo.z = z; volinfo.v = v;
    volinfo.vx = vx; volinfo.vy = vy; volinfo.vz = vz; volinfo.tr = tr;

    volinfo.miscinfo = FslInit();
    FslCloneHeader(volinfo.miscinfo,IP);

    size_t imagesize=x*y*z;
    FslGetDataType(IP,&type);
  
    output.ReSize(x*y*z);
  
    switch(type)
      {
      case DT_SIGNED_SHORT:
	{
	  short* sbuffer=new short[imagesize];
	  FslReadVolumes(IP,sbuffer,v);
	     
	  for(size_t j = 1; j<=(size_t)x*y*z; j++)
	    {
	      if (doscaling==0) { output(j)=sbuffer[j-1]; }
	      else { output(j)=(slope * sbuffer[j-1]) + intercept; }
	    }
	  
	  delete[] sbuffer;
	}
	break;
      case DT_FLOAT:
	{
	  float* fbuffer=new float[imagesize];
	  FslReadVolumes(IP,fbuffer,v);

	  for(size_t j = 1; j<=(size_t)x*y*z; j++)
	    {
	      if (doscaling==0) { output(j)=fbuffer[j-1]; }
	      else { output(j)=(slope * fbuffer[j-1]) + intercept; }
	    }
	      
	  delete[] fbuffer;
	}
	break;
      case DT_UNSIGNED_CHAR:
	{
	  unsigned char* cbuffer=new unsigned char[imagesize];
	  FslReadVolumes(IP,cbuffer,v);

	  for(size_t j = 1; j<=(size_t)x*y*z; j++)
	    {
	      if (doscaling==0) { output(j)=cbuffer[j-1]; }
	      else { output(j)=(slope * cbuffer[j-1]) + intercept; }
	    }
	      
	  delete[] cbuffer;
	}
	break;
      default:
	perror("FslRead: DT not supported");
      }

    FslClose(IP);
    
    return;
  }
}



   









