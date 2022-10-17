/*  volumeseries.cc

    Mark Woolrich - FMRIB Image Analysis Group

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
#include "volumeseries.h"
#include "miscmaths.h"
#include "utils/time_tracer.h"

using namespace NEWMAT;
using namespace Utilities;

namespace MISCMATHS {

  void VolumeSeries::replaceMeans()
  {
    Time_Tracer ts("VolumeSeries::replaceMeans");
    for(int i = 1; i <= getNumSeries(); i++)
      {
	Column(i) = Column(i) + means(i);
      }
  }

  void VolumeSeries::unthresholdSeries()
    {
       Time_Tracer ts("VolumeSeries::unthresholdSeries");

       int numUnThresholdedSeries = getUnthresholdNumSeries();
       int numThresholdedSeries = getNumSeries();
       int numVolumes = getNumVolumes();
      
       // Increase this to required size and set to 0
       Release();        
       Matrix X=*this; 
       ReSize(numVolumes, numUnThresholdedSeries);
       Matrix::operator=(0);
       
       // Loop through restoring non thresholded series: 
       for(int i = 1; i <= numThresholdedSeries; i++)
	 {
	   Column(int(preThresholdPositions(i)))=X.Column(i);
	 }
    }

  void VolumeSeries::thresholdSeries()
  {
    Time_Tracer ts("VolumeSeries::thresholdSeries");

       int numSeries = preThresholdPositions.Nrows();
       int numVolumes = getNumVolumes();

       Matrix X(numVolumes, numSeries); 
      
       for(int i = 1; i <= numSeries; i++)
	 {
	   X.Column(i)=Column(int(preThresholdPositions(i)));
	 }

       *this = X;
    }
  
   void VolumeSeries::thresholdSeries(float thresh, bool removeMean)
     {
       Time_Tracer ts("VolumeSeries::thresholdSeries");

       int numSeries = getNumSeries();
       int j = 0;

       if(removeMean)
	 {
	   means.ReSize(numSeries);
	   means = 0;
	 }

       preThresholdPositions.ReSize(numSeries);
       
       float m = 0,s = 0;
       for(int i = 1; i <= numSeries; i++)
	 {
	   m = MISCMATHS::mean(ColumnVector(getSeries(i))).AsScalar();
	   s = MISCMATHS::var(ColumnVector(getSeries(i))).AsScalar();
	   /*
	   if(m > thresh && s == 0.0)
	     {
	       cerr << "m = " << m << endl;
	       cerr << "s = " << s << endl;
	       cerr << "i = " << i << endl;
	       cerr << "j = " << j+1 << endl;
	     }
	   */

	   if(m > thresh && s > 1e-10)
	     {
	       j++;	
	       preThresholdPositions(j) = i;
	     
	       if(removeMean)
		 {
		   Column(j) = getSeries(i) - m;
		   means(i) = m;
		 }
	       else
		 Column(j) = getSeries(i);
	     }
	 }

       // Discard rest:
       *this = Columns(1,j);

       preThresholdPositions = preThresholdPositions.Rows(1,j);
     }

  void VolumeSeries::removeSeriesMeans()
    {
      int numSeries = getNumSeries();
    
      for(int i = 1; i <= numSeries; i++)
	{
	  float m = MISCMATHS::mean(ColumnVector(getSeries(i))).AsScalar();
	  Column(i) = getSeries(i) - m;
	}
    }

  void VolumeSeries::read(const string& fname)
    {
      Time_Tracer ts(string("VolumeSeries::read-" + fname).c_str());

      FSLIO* IP = FslOpen(fname.c_str(), "rb");
      Matrix& output = *this;

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
  
      output.ReSize(v,x*y*z);
  
      for(size_t i=1; i<=(size_t)v; i++)
	{
	  switch(type)
	    {
	    case DT_SIGNED_SHORT:
	      {
		short* sbuffer=new short[imagesize];
		FslReadVolumes(IP, sbuffer,1);
		for(size_t j = 1; j<=imagesize; j++)
		  {
		    if (doscaling==0) { output(i,j)=sbuffer[j-1]; }
		    else { output(i,j)=(slope * sbuffer[j-1]) + intercept;}
		  }
		delete[] sbuffer;
	      }
	      break;
	    case DT_FLOAT:
	      {
		float* fbuffer=new float[imagesize];
		FslReadVolumes(IP,fbuffer,1);
		for(size_t j = 1; j<=imagesize; j++)
		  {
		    if (doscaling==0) { output(i,j)=fbuffer[j-1]; }
		    else { output(i,j)=(slope * fbuffer[j-1]) + intercept;}
		  }
		delete[] fbuffer;
	      }
	      break;
	    case DT_UNSIGNED_CHAR:
	      {
		unsigned char* cbuffer=new unsigned char[imagesize];
		FslReadVolumes(IP,cbuffer,1);
		for(size_t j = 1; j<=imagesize; j++)
		  {
		    if (doscaling==0) { output(i,j)=cbuffer[j-1]; }
		    else { output(i,j)=(slope * cbuffer[j-1]) + intercept;}
		  }
		delete[] cbuffer;
	      }
	      break;
	    default:
	      perror("FslRead: DT not supported");
	    }
	}

      FslClose(IP);

      return;
    }
  
  void VolumeSeries::unthresholdSeries(const VolumeInfo& pvolinfo,const ColumnVector& in)
  {
    volinfo = pvolinfo;
    preThresholdPositions = in;
    unthresholdSeries();
  }

  void VolumeSeries::writeAsFloat(const string& fname)
    {
       Time_Tracer ts(string("VolumeSeries::writeAsFloat" + fname).c_str());

       FSLIO* OP = FslOpen(fname.c_str(), "wb");

       FslCloneHeader(OP,volinfo.miscinfo);
       
       FslSetDim(OP,volinfo.x, volinfo.y, volinfo.z, volinfo.v);
       FslSetVoxDim(OP,volinfo.vx, volinfo.vy, volinfo.vz, volinfo.tr);
       FslSetDataType(OP, DT_FLOAT);
       FslSetIntent(OP, volinfo.intent_code, volinfo.intent_p1, volinfo.intent_p2, 
		    volinfo.intent_p3);

       int volSize = getNumSeries();
       int volNum = getNumVolumes();

       FslWriteHeader(OP);

       float *qv = new float[volSize];
     
       for(int i = 1; i<= volNum; i++)
	 {
	    for(int j = 1; j <= volSize; j++)
	      qv[j-1] = (*this)(i,j);
	    FslWriteVolumes(OP, qv, 1);
	 }

       delete [] qv;

       FslClose(OP);
       
    }

  void VolumeSeries::writeThresholdedSeriesAsFloat(const VolumeInfo& pvolinfo,const ColumnVector& in,const string& fname)
    {
       volinfo = pvolinfo;
       preThresholdPositions = in;

       Time_Tracer ts(string("VolumeSeries::writeThresholdedSeriesAsFloat" + fname).c_str());

       FSLIO* OP = FslOpen(fname.c_str(), "wb");

       FslCloneHeader(OP,volinfo.miscinfo);
       
       FslSetDim(OP,volinfo.x, volinfo.y, volinfo.z, volinfo.v);
       FslSetVoxDim(OP,volinfo.vx, volinfo.vy, volinfo.vz, volinfo.tr);
       FslSetDataType(OP, DT_FLOAT);
       FslSetIntent(OP, volinfo.intent_code, volinfo.intent_p1, volinfo.intent_p2, volinfo.intent_p3);

       int volSize = getUnthresholdNumSeries();
       int numThresholdedSeries = getNumSeries();
       int volNum = getNumVolumes();

       FslWriteHeader(OP);

       float *qv = new float[volSize];
      
       for(int i = 1; i<= volNum; i++)
	 {
	   for(int j = 1; j <= volSize; j++)
	     qv[j-1]=0;
	   for(int j = 1; j <= numThresholdedSeries; j++)
	     qv[int(preThresholdPositions(j))-1] = (*this)(i,j);
	   FslWriteVolumes(OP, qv, 1);
	 }

       delete [] qv;

       FslClose(OP);
       
    }

}



