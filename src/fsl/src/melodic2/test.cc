/*  test.cc
    
    Copyright (C) 1999-2004 University of Oxford */

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

#include <iomanip>
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"
#include <iostream>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include <string>
#include "utils/log.h"
#include "melgmix.h"
#include "meloptions.h"
#include <time.h>
#include "miscmaths/miscprob.h"


using namespace std;
using namespace Utilities;
using namespace NEWIMAGE;
using namespace MISCPLOT;
using namespace MISCPIC;
using namespace Melodic;


Matrix calcFFT(const Matrix& Mat)
  {
    Matrix res;
    for(int ctr=1; ctr <= Mat.Ncols(); ctr++)
      {
	ColumnVector tmpCol;
	tmpCol=Mat.Column(ctr);
	ColumnVector FtmpCol_real;
	ColumnVector FtmpCol_imag;
	ColumnVector tmpPow;
	if(tmpCol.Nrows()%2 != 0){
	  Matrix empty(1,1); empty=0;
	  tmpCol &= empty;}
	RealFFT(tmpCol,FtmpCol_real,FtmpCol_imag);
	tmpPow = pow(FtmpCol_real,2)+pow(FtmpCol_imag,2);
	tmpPow = tmpPow.Rows(2,tmpPow.Nrows());
	//if(opts.logPower.value()) tmpPow = log(tmpPow);
	if(res.Storage()==0){res= tmpPow;}else{res|=tmpPow;}
      }
    return res;
  } //Matrix calc_FFT()

string float2str(float f, int width, int prec, int scientif)
	{
	    ostringstream os;
	    int redw = int(std::abs(std::log10(std::abs(f))))+1;
	    if(width>0)
	      os.width(width);
	    if(scientif>0)
	      os.setf(ios::scientific);
	    os.precision(redw+std::abs(prec));
	    os.setf(ios::internal, ios::adjustfield);
	    os << f';
	    return os.str();
	}

void usage(void)
{
  cout << "Usage: test infile ICfile [mixfile]" << endl;
  exit(1);
}

ReturnMatrix normpdf3(const RowVector& vals, const RowVector& mu, const RowVector& var)
{
  Matrix res(mu.Ncols(),vals.Ncols());
  for (int mc=1; mc<=res.Ncols(); mc++){
    for (int mr=1; mr<=res.Nrows(); mr++){
      res(mr,mc) = std::exp(-0.5*(std::pow(vals(mc)-mu(mr),2)/var(mr)))*std::pow(2*M_PI*var(mr),-0.5);
    }
  }

  res.Release();
  return res;
}




int main(int argc, char *argv[])
{

  RowVector mutmp(1);
  RowVector vartmp(1);

  RowVector valtmp(1);

  mutmp(1) =0.0;
  vartmp(1) =1.0;

  float mintmp=-1;
  float inttmp=0.1;

  for(int ctr=1; ctr<=5; ctr++){
    valtmp(1) = mintmp + ctr*inttmp;
    cerr << valtmp << "    " << normpdf(valtmp,mutmp,vartmp) << endl;
    // cerr << valtmp << "    " << normpdf(valtmp,mutmp,vartmp) << endl;
  }

  //exit(1);

  if (argc<3)
    usage();

  Matrix ICs;
  Matrix mixMatrix;
  Matrix fmixMatrix;
  volumeinfo ICvolInfo;
  volume<float> Mask;
  volume<float> Mean;

  string RAWfname;
  RAWfname = string(argv[1]);
  string ICfname;
  ICfname = string(argv[2]);

  string MIXfname;
  if (argc>3)
    MIXfname = string(argv[3]);

  string MASKfname;
  if (argc>4)
    MASKfname = string(argv[4]);



  cerr << argc << "  " << RAWfname << " " <<ICfname << " " << MIXfname << MASKfname << endl;
 
  /* {
    volume4D<float> RawData;
    cout << " Reading orig. data " << RAWfname << " ... ";
    read_volume4D(RawData,RAWfname,ICvolInfo);
    Mean = meanvol(RawData);

    float howmuch = 0.5*(std::min(std::min(std::abs(Mean.xdim()),std::abs(Mean.ydim())),std::abs(Mean.zdim())));

    cout << " Smoothing by " << howmuch << endl;
    volume<float> tmpvol = smooth(Mean,howmuch);
 
    miscpic newpic;
    char instr[10000];

    sprintf(instr," ");
    strcat(instr,"-s 2");
    strcat(instr," -A 950 ");
    strcat(instr,string("./res/m1.png").c_str());      
    newpic.slicer(Mean, instr, &ICvolInfo); 

    char instr2[10000];

    sprintf(instr2," ");
    strcat(instr2,"-s 2");
    strcat(instr2," -A 950 ");
    strcat(instr2,string("./res/m2.png").c_str());      
    newpic.slicer(tmpvol, instr2, &ICvolInfo); 



    cout << " done! " << endl;
    }*/

   {
    volume4D<float> RawIC;
    cout << "Reading components " << ICfname << "  ... ";
    read_volume4D(RawIC,ICfname);
    cout << " done" << endl;

    cout << "Creating mask   ... ";

    read_volume(Mask,MASKfname);

    //  Mask = binarise(RawIC[0],float(RawIC[0].min()),float(RawIC[0].max()));

    ICs = RawIC.matrix(Mask);
    /*if(ICs.Nrows()>1){
      Matrix DStDev=stdev(ICs);
      
      volume4D<float> tmpMask;
      tmpMask.setmatrix(DStDev,Mask);
      
      float tMmax;
      volume<float> tmpMask2;
      tmpMask2 = tmpMask[0];
      tMmax = tmpMask2.max();
      double st_mean = DStDev.Sum()/DStDev.Ncols();
      double st_std  = stdev(DStDev.t()).AsScalar();
      
      Mask = binarise(tmpMask2,(float) max((float) st_mean-3*st_std,
					   (float) 0.01*st_mean),tMmax);  
      ICs = RawIC.matrix(Mask);
    }
    else{
      Mask = binarise(RawIC[0],float(0.001),float(RawIC[0].max())) 
	+ binarise(RawIC[0],float(RawIC[0].min()),float(-0.001));
      ICs = RawIC.matrix(Mask);
      }*/

    cerr << "ICs : " << ICs.Ncols() << ICs.Nrows() << endl;
    cout << " done" << endl;
  }

  if(MIXfname.length()>0){
    cout << "Reading mixing matrix " << MIXfname << " ... ";
    mixMatrix = read_ascii_matrix(MIXfname);
    if (mixMatrix.Storage()<=0) {
      cerr <<" Please specify the mixing matrix correctly" << endl;
      exit(2);
      cout << " done " << endl;
    }
  }

  cout << " ICs: " << ICs.Nrows() << " x " << ICs.Ncols() << endl;
for(int ctr=1; ctr <= ICs.Nrows(); ctr++){ 

  cout << " Plotting histogram for map " << ctr << endl;
  
    miscplot newplot;
    // newplot.add_label("legend1");
    //newplot.add_label("legend2");
    //newplot.add_ylabel(string("yll1"));
    
    //newplot.add_xlabel(string("xl1"));
    //    newplot.set_xysize(600,600);

    Matrix mu(1,3);
    Matrix pi(1,3);
    Matrix std(1,3);
    mu(1,1)=0;mu(1,2)=1.10;mu(1,3)=-1.1;
    pi(1,1)=0.96;pi(1,2)=0.02;pi(1,3)=0.02;
    std(1,1)=1;std(1,2)=0.1;std(1,3)=0.1;

    //   Matrix returnM;
    // returnM = normpdf(ICs.Row(ctr)*10,10,11);

    newplot.gmmfit(ICs.Row(ctr),mu,std,pi,string("./res/g"+num2str(ctr)+".png"),string("Title"));
    //newplot.histogram(ICs.Row(ctr),string("./res/g"+num2str(ctr)+".png"),string("Title"));

  }

  cout << endl << endl;

  return 0;
}
