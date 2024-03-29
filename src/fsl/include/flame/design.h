/*  design.h

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


#if !defined(design_h)
#define design_h
  
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "newmatap.h"
#include "newmatio.h"

using namespace NEWMAT;

namespace Gs{

  class Design
    {
    public:

      // constructor
      Design() :	
	nevs(0),
	ntpts(0),
	ngs(0)
	{ 
	}

      ~Design() {}

      // load design matrix in from file and set it up
      void setup(bool loadcontrasts = true);

      // getters
      const int getnevs() const { return nevs; }
      const int getntpts() const { return ntpts; }
      const int getngs() const { return ngs; }

      const Matrix& getdm() const { return dm; }
      const Matrix& getcovsplit() const { return covsplit; }

      ReturnMatrix gettcontrast(const int p_num) const { RowVector ret = tcontrasts.Row(p_num); ret.Release(); return ret; }

      const Matrix& getfcontrast(const int p_num) const { return fc[p_num-1]; }

      int getgroup(int t) const { return int(gind(t)); } 
      int getntptsingroup(int g) const { return int(ntptsing(g)); }

      int getnumfcontrasts() const { return numFcontrasts; }
      int getnumtcontrasts() const { return numTcontrasts; }
      //const Matrix& getfreduceddm(const int p_num) const { return reduceddms[p_num-1]; }
    private:

      const Design& operator=(Design& par);     
      Design(Design& des) { operator=(des); }
      
      void setupfcontrasts();

      int nevs;
      int ntpts;
      int ngs;

      Matrix dm;
      Matrix covsplit;

      Matrix tcontrasts;
      Matrix fcontrasts;
      vector<Matrix> fc;

      //      vector<Matrix> reduceddms;

      ColumnVector gind;
      ColumnVector ntptsing;

      int numFcontrasts;
      int numTcontrasts;
    };
} 
#endif


