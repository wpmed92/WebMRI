/*  paradigm.cc

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

#include "paradigm.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include <fstream>

using namespace FILM;

namespace FILM {

  void Paradigm::read_vest_waveform(string p_fname, Matrix& p_mat)
  {
    ifstream in;
    in.open(p_fname.c_str(), ios::in);
    
    if(!in)
      {
	cerr << "Unable to open " << p_fname << endl;
	throw;
      }
    
    int numWaves = 0;
    int numPoints = 0;
  
    string str;

    while(true)
      {
	in >> str;
	if(str == "/Matrix")
	  break;
	else if(str == "/NumWaves")
	  {
	    in >> numWaves;
	  }
	else if(str == "/NumPoints" || str == "/NumContrasts")
	  {
	    in >> numPoints;
	  }
      }

    if(numWaves != 0)
      {
	p_mat.ReSize(numPoints, numWaves);
      }
    else
      {
	numWaves = p_mat.Ncols();
	numPoints = p_mat.Nrows();
      }
      
    for(int i = 1; i <= numPoints; i++)
      {
	for(int j = 1; j <= numWaves; j++)    
	  {
	    in >> p_mat(i,j);
	  }
      }

    in.close();
  }
  
  void Paradigm::load(const string& p_paradfname, const string& p_tcontrastfname, const string& p_fcontrastfname, bool p_blockdesign, int p_sizets)
    {
      if(p_paradfname != "")
	read_vest_waveform(p_paradfname, designMatrix);

      if(p_tcontrastfname != "")
	read_vest_waveform(p_tcontrastfname, tcontrasts);       

      if(p_fcontrastfname != "")
	read_vest_waveform(p_fcontrastfname, fcontrasts);   
       
      // Check time series match up with design
      if(p_paradfname != "" && designMatrix.Nrows() != p_sizets)
	{
	  cerr << "Num scans = "<< p_sizets << ", design matrix rows = " << designMatrix.Nrows() << endl;
	  throw Exception("size of design matrix does not match number of scans\n");
	}
      
      // Check contrasts match
      if(p_tcontrastfname != "" &&  p_fcontrastfname != "" && tcontrasts.Nrows() != fcontrasts.Ncols())
	{ 
	  cerr << "tcontrasts.Nrows()  = " << tcontrasts.Nrows() << endl;
	  cerr << "fcontrasts.Ncols()  = " << fcontrasts.Ncols() << endl;
	  throw Exception("size of tcontrasts does not match fcontrasts\n");
	}

      // Check contrast matches design matrix
      if(p_paradfname != "" &&  p_tcontrastfname != "" && designMatrix.Ncols() != tcontrasts.Ncols())
	 { 
	  cerr << "Num tcontrast cols  = " << tcontrasts.Ncols() << ", design matrix cols = " 
	    << designMatrix.Ncols() << endl;
	  throw Exception("size of design matrix does not match t contrasts\n");
	 }

      if(p_blockdesign)
	{
	  // Need to adjust amplitude from (-1,1) to (0,1)
	  designMatrix = (designMatrix + 1)/2;
	}
    }
  
}







