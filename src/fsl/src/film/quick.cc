/*  quick.cc

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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include <string>
#include <math.h>
#include "glmrand.h"
#include "miscmaths/volume.h"
#include "utils/log.h"
#include "histogram.h"
#include "miscmaths/t2z.h"
#include "miscmaths/f2z.h"
#include "miscmaths/miscmaths.h"
#include "libprob.h"
#include <time.h>

using namespace Utilities;
using namespace NEWMAT;
using namespace FILM;
using namespace MISCMATHS;

int main(int argc, char *argv[])
{
  try{
    rand();

//     Log::getInstance().setDir(".");

//     //    cerr << "z = " << T2z::getInstance().convert(9.02,958) << endl;
//     int X = 20;
//     int Y = 20;
//     int Z = 1;
//     int T = 200;
    
//     int N = X*Y*Z*T;
//     //Matrix mat(N,N);
//     SymmetricBandMatrix mat(N,X*Y*Z);
//     mat = 0;

//     cerr << "X*Y*Z*T=" << N << endl;
//     cerr << "X*Y*Z=" << X*Y*Z << endl;
//     int mba = 0;

//     for(float j=0;j<=0.8;j+=0.4)
//       {
// 	for(float i=-0.8;i<=0.8;i+=0.8)
// 	  {
	
//       for(int t = 1; t <= T; t++)
//       for(int z = 1; z <= Z; z++)
// 	for(int y = 1; y <= Y; y++)
// 	  for(int x = 1; x <= X; x++)
// 	    {
// 	      int a = x+(y-1)*X+(z-1)*X*Y+(t-1)*Y*X*Z;

// 	      for(int t2 = Max(t-1,1); t2 <= Min(t+1,T); t2++)
// 		for(int z2 = Max(z-1,1); z2 <= Min(z+1,Z); z2++)
// 		  for(int y2 = Max(y-1,1); y2 <= Min(y+1,Y); y2++)
// 		    for(int x2 = Max(x-1,1); x2 <= Min(x+1,X); x2++)
// 		      {		    
// 			int b = x2+(y2-1)*X+(z2-1)*X*Y+(t2-1)*Y*X*Z;						
// 			if(b-a > mba)
// 			  mba = b-a;

// // 			cerr <<"t="<<t;
// // 			cerr <<"z="<<z;
// // 			cerr <<"y="<<z;
// // 			cerr <<"x="<<x;
// // 			cerr <<"t2="<<t2;
// // 			cerr <<"z2="<<z2;
// // 			cerr <<"y2="<<z2;
// // 			cerr <<"x2="<<x2<<endl;

// 			if (x==x2 && y==y2 && z==z2 && t==t2)
// 			  {			    
// 			    mat(a,b) = 1;
// 			  }
// 			else
// 			  {
// 			    if(t==t2)
// 			      mat(a,b) = -j;
// 			    else if(x==x2 && y==y2 && z==z2)
// 			      mat(a,b) = -i;
// 			  }			
// 		      }
// 	    }
//       //cerr << "mba=" << mba <<endl;
//       //write_ascii_matrix(mat,"mat");
      
//       clock_t start = clock();
	        
//       //cerr << "det(mat)=" << mat.Determinant() << endl;
      
//       cerr << "logdet(mat)=" << mat.LogDeterminant().LogValue() << endl;
    //	  }
    //      }
      

    //      cerr << "time taken=" << (clock()-start)/float(CLOCKS_PER_SEC) << endl;
  }
  catch(Exception p_excp) 
    {
      cerr << p_excp.what() << endl;
    }
  return 0;
}












