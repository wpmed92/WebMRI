/*  glim.cc

    Mark Woolrich, FMRIB Image Analysis Group

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

#include <sstream>

#include "glim.h"
#include "miscmaths/miscmaths.h"
#include "utils/log.h"
#include "miscmaths/volume.h"

#ifndef NO_NAMESPACE
using namespace MISCMATHS;
using namespace Utilities;

namespace FILM {
#endif

  Glim::Glim(VolumeSeries& p_y, const Matrix& p_x):
    y(p_y),
    x(p_x),
    numTS(p_y.Ncols()),
    sizeTS(p_y.Nrows()),
    numParams(p_x.Ncols()),
    r(sizeTS, numTS, p_y.getInfo(), p_y.getPreThresholdPositions()),
    pinv_x(p_x.Ncols(), sizeTS),
    V(sizeTS,sizeTS),
    RV(sizeTS,sizeTS),
    RMat(sizeTS,sizeTS),
    batch_size(BATCHSIZE),
    corrections(numTS, numParams*numParams),
    b(numTS, numParams),
    sigmaSquareds(numTS),
    dof(sizeTS - p_x.Ncols())
    { 
    }

  void Glim::Save()
    {
      // Need to save b, sigmaSquareds, corrections and dof 
      Log& logger = LogSingleton::getInstance();
     
      // b:
      Volume peVol;
      for(int i = 1; i <= numParams; i++)
	{
	  peVol = b.Row(i).AsColumn();
	  peVol.setInfo(y.getInfo());
	  peVol.setPreThresholdPositions(y.getPreThresholdPositions());
	  peVol.unthreshold();	
	  
	  // Add param number to "pe" to create filename:
	  ostringstream osc;
	  osc << i;
	  
	  peVol.writeAsFloat(logger.getDir() + "/pe" + osc.str().c_str());
	}

      // sigmaSquareds:
      sigmaSquareds.setInfo(y.getInfo());
      sigmaSquareds.setPreThresholdPositions(y.getPreThresholdPositions());
      sigmaSquareds.unthreshold();	
      sigmaSquareds.writeAsFloat(logger.getDir() + "/sigmasquareds");

      // dof:
      ColumnVector dofVec(1);
      dofVec = dof;
      write_ascii_matrix(logger.appendDir("dof"), dofVec);

      // corrections:
      write_ascii_matrix(logger.appendDir("corrections"), corrections);
    }

  // Called on entire data set:
  VolumeSeries& Glim::ComputeResids()
    {
      Tracer ts("ComputeResids");

      int batch_pos = 1;

      // pinv(x) = inv(x'x)x'
      //pinv_x = (x.t()*x).i()*x.t();
      pinv_x = pinv(x);

      // R = I - x*pinv(x)
      Matrix I(sizeTS, sizeTS);
      Identity(I);

      RMat = I - x*pinv_x;
      
      while(batch_pos <= numTS)
	{
	  if(batch_pos+batch_size - 1 > numTS)
	    r.Columns(batch_pos, numTS) = RMat*y.Columns(batch_pos, numTS);
	  else
	    r.Columns(batch_pos, batch_pos+batch_size-1) = RMat*y.Columns(batch_pos, batch_pos+batch_size-1);
	
	  batch_pos += batch_size;
	}
      
      return r;
    }

  // Called on entire data set:
  void Glim::ComputePes()
    { 
      Tracer ts("ComputePe");
      
      b = pinv_x*y;
    }

  void Glim::ConstructV(const ColumnVector& p_vrow)
    {
      Tracer ts("ConstructV");
      V = 0;

      for (int i = 1; i <= sizeTS; i++)
	{
	  V.SubMatrix(i,i,i,sizeTS) = p_vrow.Rows(1,sizeTS-i+1).t();
	  V.SubMatrix(i,i,1,i) = p_vrow.Rows(1,i).Reverse().t();
	}
    }

  void Glim::SetVrow(const ColumnVector& p_vrow, const int ind)
    {
      Tracer ts("SetVrow");
      
      ConstructV(p_vrow);

      // var/e = inv(x'x)x'*V*x*inv(x'x)
      Matrix corr = pinv_x*V*pinv_x.t();
      SetCorrection(corr, ind);
    }

  void Glim::SetGlobalVrow(const ColumnVector& p_vrow)
    {
      Tracer ts("Glim::SetGlobalVrow");

      ConstructV(p_vrow);
      RV = RMat*V;

      dof = Trace(RV)*Trace(RV)/Trace(RV*RV);
    }


  void Glim::UseGlobalVrow()
    {
      Tracer ts("Glim::UseGlobalVrow");

      // var/e = inv(x'x)x'*V*x*inv(x'x)
      Matrix corr = pinv_x*V*pinv_x.t();

      for(int i = 1; i <= numTS; i++)
	{
	  SetCorrection(corr, i);
	}
    }

  void Glim::ComputeSigmaSquared(const int ind)
    {
      Tracer ts("Glim::ComputeSigmaSquared");

      sigmaSquareds(ind) = (r.Column(ind).t()*r.Column(ind)/Trace(RV)).AsScalar();
    }

  void Glim::SetCorrection(const Matrix& corr, const int ind)
    {
      Tracer ts("SetCorrection");

      // puts Matrix corr which is p*p into Matrix correction
      // as a p*p length row:
      int p = corr.Nrows();
      
      for (int i = 1; i <= p; i++)
	{
	  for (int j = 1; j <= p; j++)
	    {
	      corrections(ind, (i-1)*p + j) = corr(i,j); 
	    }
	}
    }

#ifndef NO_NAMESPACE
}
#endif







