/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    meldata.h - data container class

    Christian F. Beckmann, FMRIB Image Analysis Group
    
    Copyright (C) 1999-2002 University of Oxford */

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


#ifndef __MELODICDATA_h
#define __MELODICDATA_h

#include "newimage/newimageall.h"
#include "utils/log.h"
#include "meloptions.h"

using namespace Utilities;
using namespace NEWIMAGE;

namespace Melodic{
  
  class MelodicData
    {
    public:

      //constructor
      MelodicData(MelodicOptions &popts, Log &plogger):  
	opts(popts),
	logger(plogger)
	{	
	  after_mm = false;
	}  
 
      void save();

      inline void save4D(Matrix what, string fname){
	 volume4D<float> tempVol;
	 tempVol.setmatrix(what,Mask);
	 save_volume4D(tempVol,logger.appendDir(fname),tempInfo);
      }
      
      inline void saveascii(Matrix what, string fname){
	 write_ascii_matrix(logger.appendDir(fname),what);      }

      int  remove_components();
      void setup();
      void status(const string &txt);

      inline Matrix& get_pcaE() {return pcaE;}
      inline void set_pcaE(Matrix& Arg) {pcaE = Arg;}

      inline DiagonalMatrix& get_pcaD() {return pcaD;}
      inline void set_pcaD(DiagonalMatrix& Arg) {pcaD = Arg;}

      inline Matrix& get_data() {return Data;}
      inline void set_data(Matrix& Arg) {Data = Arg;}

      inline Matrix& get_IC() {return IC;}
      inline void set_IC(Matrix& Arg) {IC = Arg;}
      inline void set_IC(int ctr, Matrix& Arg) {IC.Row(ctr) = Arg;}
      
      inline Matrix& get_white() {return whiteMatrix;}
      inline void set_white(Matrix& Arg) {whiteMatrix = Arg;}
      
      inline Matrix& get_dewhite() {return dewhiteMatrix;}
      inline void set_dewhite(Matrix& Arg) {dewhiteMatrix = Arg;}
      
      inline Matrix& get_meanC() {return meanC;}
      inline Matrix& get_meanR() {return meanR;}

      inline Matrix& get_stdDevi() {return stdDevi;}
      inline void set_stdDevi(Matrix& Arg) {stdDevi = Arg;}
  
      inline Matrix& get_mix() {return mixMatrix;}
      inline void set_mix(Matrix& Arg) 
	{mixMatrix = Arg; unmixMatrix=pinv(mixMatrix);}
      inline void set_mix(int ctr, Matrix& Arg) 
	{mixMatrix.Column(ctr) = Arg; unmixMatrix=pinv(mixMatrix);}
   
      inline Matrix& get_fmix() {return mixFFT;}
      inline void set_fmix(Matrix& Arg) {mixFFT = Arg;}

      inline Matrix& get_unmix() {return unmixMatrix;}
      inline void set_unmix(Matrix& Arg) {unmixMatrix = Arg;}

      inline volume<float>& get_mask() {return Mask;}
      inline void set_mask(volume<float>& Arg) {Mask = Arg;}
  
      inline volume<float>& get_mean() {return Mean;}
      inline void set_mean(volume<float>& Arg) {Mean = Arg;}
   
      inline Matrix& get_Data() {return Data;}
      inline void set_Data(Matrix& Arg) {Data = Arg;}
    
      inline Matrix& get_DataVN() {return DataVN;}
      inline void set_DataVN(Matrix& Arg) {DataVN = Arg;}

      inline Matrix& get_RXweight() {return RXweight;}
      inline void set_RXweight(Matrix& Arg) {RXweight = Arg;}
 
      inline Matrix& get_ICstats() {return ICstats;}
      inline void set_ICstats(Matrix& Arg) {ICstats = Arg;}
     
      inline Matrix& get_EVP() {return EVP;}
      inline void set_EVP(Matrix& Arg) {EVP = Arg;}
      
      inline Matrix& get_EV() {return EV;}
      inline void set_EV(Matrix& Arg) {EV = Arg;}

      inline Matrix& get_PPCA() {return PPCA;}
      inline void set_PPCA(Matrix& Arg) {PPCA = Arg;}

      inline Matrix& get_stdNoisei() {return stdNoisei;}
      inline void set_stdNoisei(Matrix& Arg) {stdNoisei = Arg;}

      inline int data_dim() {return Data.Nrows();}
      inline int data_samples() {return Data.Ncols();}
     
      inline float get_resels() {return Resels;}
      inline void set_resels(float& Arg) {Resels = Arg;}

      inline void set_after_mm(bool val) {after_mm = val;}

      Matrix smoothColumns(const Matrix &inp);

      inline void flipres(int num){
	IC.Row(num) = -IC.Row(num);
	mixMatrix.Column(num) = -mixMatrix.Column(num);
	mixFFT=calc_FFT(mixMatrix);
	unmixMatrix = pinv(mixMatrix);

        if(ICstats.Storage()>0&&ICstats.Ncols()>3){
	  double tmp;
	  tmp = ICstats(num,3);
	  ICstats(num,3) = -1.0*ICstats(num,4);
	  ICstats(num,4) = -1.0*tmp;
	}
      }

      
      inline void varnorm(){
	Data = SP(Data,ones(Data.Nrows(),1)*stdDevi);
      }

      void sort();

      volumeinfo tempInfo;
      
      Matrix calc_FFT(const Matrix& Mat);

    private:
      MelodicOptions &opts;     
      Log &logger;       
      //      static MelodicData* melodat;

      Matrix pcaE;
      DiagonalMatrix pcaD;
      Matrix whiteMatrix;
      Matrix dewhiteMatrix;
      Matrix meanC;
      Matrix meanR;
      Matrix stdDev;
      Matrix stdDevi;
      Matrix RXweight;
      Matrix mixMatrix;
      Matrix unmixMatrix;
      Matrix mixFFT;
      Matrix IC;
      Matrix ICstats;
      Matrix EVP;
      Matrix EV;
      Matrix stdNoisei;
      volume<float> Mask;
      volume<float> Mean;
      Matrix Data;
      Matrix DataVN;
      Matrix PPCA;
      
      bool after_mm;

      float Resels;

      char Mean_fname[1000];

      void create_mask(volume4D<float> &theData, volume<float> &theMask);
      //      void remove_comp(Matrix& Mat, const int numIC);
      void create_RXweight();
         
    };
}

#endif
