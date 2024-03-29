/*  Gam.h

    Mark Woolrich, Tim Behrens - FMRIB Image Analysis Group

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

#if !defined(Gam_h)
#define Gam_h

#include <iostream>
#include <math.h>
//extern "C"
//{
#include "libprob.h"
  //};
#include "utils/tracer_plus.h"

using namespace Utilities;

namespace Gs {
  
  class Gam
    {
    public:
      static Gam& getInstance();
      ~Gam() 
	{ delete gam; 
	delete gamma; 
	}
      
      void setParams(const double& pa, const double& pb) {
	b = pb; 
	if(pa!=a) 
	  {
	    a=pa; 
	    delete gamma; gamma=new Gamma(a);
	  }
      }
	
      double pdf(const double& x);
      double neglog_pdf(const double& x);
      double cdf(const double& x);
      double rnd();      

      double pdf(const double& pa, const double& pb, const double& x) {a = pa; b = pb; return pdf(x);}
      double neglog_pdf(const double& pa, const double& pb, const double& x) {a = pa; b = pb; return neglog_pdf(x);}
      double cdf(const double& pa, const double& pb, const double& x) {a = pa; b = pb; return cdf(x);}
      double rnd(const double& pa, const double& pb) {	
       	setParams(pa,pb); 
       	return rnd();
      }

    private:
      Gam():a(0.0),b(0.0),gamma(NULL) 
	{}
      
      const Gam& operator=(Gam&);
      Gam(Gam&);
      
      static Gam* gam;

      double a; 
      double b;
      Gamma* gamma;
    };

  inline Gam& Gam::getInstance(){
    if(gam == NULL)	
      gam = new Gam();
      
    return *gam;
  }

  inline double Gam::rnd()
    {
      return gamma->Next()/b;
    }

  inline double Gam::pdf(const double& x)
    {
      if(x<=0)
	return 0;
      else
	return exp(a*log(b) + (a-1)*log(x) - b*x - lgam(a));
    } 

  inline double Gam::neglog_pdf(const double& x)
    {
      if(x<=0)
	return 1e32;
      else
	return -(a*log(b) + (a-1)*log(x) - b*x - lgam(a));
    }      
}

#endif
