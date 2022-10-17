/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    meldata.cc - data handler / container class

    Christian F. Beckmann, FMRIB Image Analysis Group
    
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

#include "newimage/newimageall.h"
#include "meloptions.h"
#include "meldata.h"
#include "melodic.h"
#include "utils/log.h"
#include <time.h>
#include "miscmaths/miscprob.h"
#include <sstream>

using namespace Utilities;
using namespace NEWIMAGE;

namespace Melodic{

  void MelodicData::setup()
  {
    {
      volume4D<float> RawData;
      message("Reading data file " << opts.inputfname.value() << "  ... ");
      read_volume4D(RawData,opts.inputfname.value(),tempInfo);
      message(" done" << endl);
      for(int ctr=1; ctr<=opts.dummy.value(); ctr++){
	RawData.deletevolume(ctr);
      }    
    
      // calculate a Mean image and save it
      Mean = meanvol(RawData);
      tmpnam(Mean_fname); // generate a tmp name
      save_volume(Mean,Mean_fname);
      create_mask(RawData, Mask);
      if(Mask.xsize()==RawData.xsize() &&
	 Mask.ysize()==RawData.ysize() &&
	 Mask.zsize()==RawData.zsize())
	 Data = RawData.matrix(Mask);
      else{
	cerr << "ERROR:: mask and data have different dimensions  \n\n";
	exit(2);
      }
	
           
      // clean /tmp
      ostringstream osc;
      osc  << "rm " << string(Mean_fname) <<"*  ";
      system(osc.str().c_str());
         
      //mask out constant voxels
      message("Excluding voxels with constant value " << endl);
      Matrix DStDev=stdev(Data);
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
      Data = RawData.matrix(Mask);
    }

    message("Data size : " << Data.Nrows() << " x " << Data.Ncols() <<endl);

    {// remove mean volume
      //    if((opts.remove_meanvol.value()||opts.varnorm.value())){
      message(string("Removing mean image ... "));
      meanR=mean(Data);
      Data=remmean(Data); 
      message("done" << endl);
      //}else{
      // meanR=zeros(1,Data.Ncols());
      //}
    }
    {// remove mean time course
      meanC=mean(Data,2);
      Data=remmean(Data,2); 
    }

    {//switch dimension in case temporal ICA is required
      if(opts.temporal.value()){
	message(string("Switching dimensions for temporal ICA") << endl);
	Data = Data.t();
	Matrix tmp;
	tmp = meanC;
	meanC = meanR.t();
	meanR = tmp.t();
	message("Data size : " << Data.Nrows() << " x " << Data.Ncols() <<endl);
      }
    }


    {// variance-normalize the data
      DiagonalMatrix tmpD;Matrix tmpE;
      message(string("Estimating data covariance ... "));
      SymmetricMatrix Corr;
      Corr = cov(Data.t());
      EigenValues(Corr,tmpD,tmpE);
      
      Matrix RE; DiagonalMatrix RD;
      //RE = tmpE.Columns(2,Corr.Ncols());
      RE = tmpE;
      //RD << abs(tmpD.SymSubMatrix(2,Corr.Ncols()));    
      RD << abs(tmpD);
      Matrix tmpWhite;Matrix tmpDeWhite;
      tmpWhite = sqrt(abs(RD.i()))*RE.t();
      tmpDeWhite = RE*sqrt(RD);
      message("done"<<endl);
      if(opts.varnorm.value())
	message(string("Perform variance-normalisation ... "));
      Matrix WS;
      WS = tmpWhite * Data;
      for(int ctr1 =1; ctr1<=WS.Nrows(); ctr1++){
	for(int ctr2 =1; ctr2<=WS.Ncols(); ctr2++){
	  if(abs(WS(ctr1,ctr2))<3.1)
	    WS(ctr1,ctr2)=0.0;
	}
      }
      
      stdDev = stdev(Data - tmpDeWhite*WS);   
      stdDevi = pow(stdDev,-1);   
      Data = Data + meanC*ones(1,Data.Ncols());
    }
    
    DataVN = SP(Data,ones(Data.Nrows(),1)*stdDevi);

    if(opts.output_all.value()){
      volume4D<float> tempVol; 
      tempVol.setmatrix(stdDevi,Mask);
      save_volume4D(tempVol,logger.appendDir(opts.outputfname.value() 
					     + "_vn_stdev"),tempInfo);
      tempVol.setmatrix(DataVN,Mask);
      save_volume4D(tempVol,logger.appendDir(opts.outputfname.value() 
					     + "_vn"),tempInfo);
    }

    if(opts.varnorm.value()){
      Data = DataVN;
      message("done"<<endl);
    }
    {//remove row mean
      if(opts.temporal.value()){
	message(string("Removing mean image ... "));
      }else{
	message(string("Removing mean time course ... "));
      }
      meanC=mean(Data,2);
      Data=remmean(Data,2); 
      message("done"<<endl);
    }
    
    if(opts.segment.value().length()>0){
       create_RXweight();
    }
 
    //save the mask
    save_volume(Mask,logger.appendDir("mask"));

    //seed the random number generator
    double tmptime = time(NULL);
    // cerr << tmptime << endl << endl;;
    srand((unsigned int) tmptime);
  } // void setup()

  void MelodicData::save()
  {   

    //check for temporal ICA
    if(opts.temporal.value()){
      message(string("temporal ICA: transform back the data ... "));
      Matrix tmpIC = mixMatrix.t();
      mixMatrix=IC.t();
      IC=tmpIC;

      unmixMatrix=pinv(mixMatrix);
      Data = Data.t();
      tmpIC = meanC;
      meanC = meanR.t();
      meanR = tmpIC.t();
      //  whiteMatrix = whiteMatrix.t;
      //  dewhiteMatrix = dewhiteMatrix.t();
      message(string("done") << endl);
      opts.temporal.set_T(false); // Do not switch again!
    }
 
    message("Writing results to : " << endl);

    //Output IC	
    if((IC.Storage()>0)&&(opts.output_origIC.value())&&(after_mm==false)){
      volume4D<float> tempVol; 
      tempVol.setmatrix(IC,Mask);
      //strncpy(tempInfo.header.hist.aux_file,"render3",24);
      save_volume4D(tempVol,logger.appendDir(opts.outputfname.value() 
					     + "_oIC"),tempInfo);
      message("  " << logger.appendDir(opts.outputfname.value() + "_oIC") <<endl);
    }

    //Output IC -- adjusted for noise	
      if(IC.Storage()>0){
      volume4D<float> tempVol;	
	//      volumeinfo tempInfo;
	//  read_volume4D(tempVol,opts.inputfname.value(),tempInfo); 

      //Matrix ICadjust = IC;
      
      Matrix ICadjust;
      if(after_mm)
	ICadjust = IC;
      else{
	stdNoisei = pow(stdev(Data - mixMatrix * IC)*std::sqrt((float)(Data.Nrows()-1))/
			std::sqrt((float)(Data.Nrows()-IC.Nrows())),-1);

	ColumnVector diagvals;
	diagvals=pow(diag(unmixMatrix*unmixMatrix.t()),-0.5);
	ICadjust = SP(IC,diagvals*stdNoisei);
      }

      tempVol.setmatrix(ICadjust,Mask);
      //strncpy(tempInfo.header.hist.aux_file,"render3",24);
      save_volume4D(tempVol,logger.appendDir(opts.outputfname.value() 
					     + "_IC"),tempInfo);
      message("  " << logger.appendDir(opts.outputfname.value() + "_IC") <<endl);

      if(opts.output_origIC.value()){
	tempVol.setmatrix(stdNoisei,Mask);
	save_volume4D(tempVol,logger.appendDir(string("Noise_stddev_inv")),tempInfo);
	message("  " << logger.appendDir(string("Noise_stddev_inv")) <<endl);
      }
      
    }

    //Output mixMatrix
    if(mixMatrix.Storage()>0){
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_mix"),
			 mixMatrix); 
      mixFFT=calc_FFT(mixMatrix);
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_FTmix"),
			 mixFFT);      
      message("  "<<
	      logger.appendDir(opts.outputfname.value() + "_mix") <<endl);
      message("  "<<
	      logger.appendDir(opts.outputfname.value() + "_FTmix") <<endl);
    }

    //Output ICstats
    if(ICstats.Storage()>0){
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_ICstats"),
			 ICstats); 
      message("  "<<
	      logger.appendDir(opts.outputfname.value() + "_ICstats") <<endl);
    }

    //Output unmixMatrix
    if(opts.output_unmix.value() && unmixMatrix.Storage()>0){
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_unmix"),unmixMatrix);
      message("  "<<
	  logger.appendDir(opts.outputfname.value() + "_unmix") <<endl);
    }

    //Output Mask
    message("  "<< logger.appendDir("mask") <<endl);

    //Output mean
    if(opts.output_mean.value() && meanC.Storage()>0 && meanR.Storage()>0){
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_meanR"),
			 meanR);
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_meanC"),
			 meanC);
      message("  "<<
	  logger.appendDir(opts.outputfname.value() + "_meanR") <<endl);
      message("  "<<
	  logger.appendDir(opts.outputfname.value() + "_meanC") <<endl);
    }

    //Output white
    if(opts.output_white.value() && whiteMatrix.Storage()>0&&
       dewhiteMatrix.Storage()>0){
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_white"),
			 whiteMatrix);
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_dewhite"),dewhiteMatrix);
      Matrix tmp;
      tmp=calc_FFT(dewhiteMatrix);
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_FTdewhite"),tmp);
      message("  "<<
	  logger.appendDir(opts.outputfname.value() + "_white") <<endl);
      message("  "<<
	  logger.appendDir(opts.outputfname.value() + "_dewhite") <<endl);
    }

    //Output PCA
    if(opts.output_pca.value() && pcaD.Storage()>0&&pcaE.Storage()>0){
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_pcaE"),
			 pcaE);
      message("  "<<
	      logger.appendDir(opts.outputfname.value() + "_pcaE") <<endl);
      write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_pcaD"),
			 (Matrix) diag(pcaD));
      message("  "<<
	      logger.appendDir(opts.outputfname.value() + "_pcaD") <<endl);
      
      volume4D<float> tempVol;	

      Matrix PCAmaps;

      if(whiteMatrix.Ncols()==Data.Ncols()){
	PCAmaps = dewhiteMatrix.t();
      }else
	PCAmaps = whiteMatrix * Data;
     

      tempVol.setmatrix(PCAmaps,Mask);
      //strncpy(tempInfo.header.hist.aux_file,"render3",24);
      save_volume4D(tempVol,logger.appendDir(opts.outputfname.value() 
					     + "_pca"),tempInfo);
      message("  " << 
	  logger.appendDir(opts.outputfname.value() + "_pca") <<endl);
      
          }
  } //void save()
  
  int MelodicData::remove_components()
  {  
    message("Reading " << opts.filtermix.value() << endl) 
    mixMatrix = read_ascii_matrix(opts.filtermix.value());
    if (mixMatrix.Storage()<=0) {
      cerr <<" Please specify the mixing matrix correctly" << endl;
      exit(2);
    }
    
    unmixMatrix = pinv(mixMatrix);
    IC = unmixMatrix * Data;

    string tmpstr;
    tmpstr = opts.filter.value();

    Matrix noiseMix;
    Matrix noiseIC;

    int ctr=0;    
    char *p;
    char t[1024];
    const char *discard = ", [];{(})abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*_-=+|\':><./?";
    
    message("Filtering the data...");
    strcpy(t, tmpstr.c_str());
    p=strtok(t,discard);
    ctr = atoi(p);
    if(ctr>0 && ctr<=mixMatrix.Ncols()){
      message(" "<< ctr );
      noiseMix = mixMatrix.Column(ctr);
      noiseIC  = IC.Row(ctr).t();    
    }else{
      cerr << endl<< "component number "<<ctr<<" does not exist" << endl;
    }
    
    do{
      p=strtok(NULL,discard);
      if(p){
	ctr = atoi(p);
	
        if(ctr>0 && ctr<=mixMatrix.Ncols()){
	  message(" "<<ctr);
	  noiseMix |= mixMatrix.Column(ctr);
	  noiseIC  |= IC.Row(ctr).t();
	}
	else{
	  cerr << endl<< "component number "<<ctr<<" does not exist" << endl;
	}
      }
    }while(p);
    message(endl);
    Matrix newData;
    newData = Data - noiseMix * noiseIC.t();

    //cerr << newData.Nrows() << " x " << newData.Ncols() << endl;
    //cerr << meanC.Nrows() << " x " << meanC.Ncols() << endl;
    //cerr << meanR.Nrows() << " x " << meanR.Ncols() << endl;
    newData = newData + meanC*ones(1,newData.Ncols());
    newData = newData + ones(newData.Nrows(),1)*meanR;
    
    volume4D<float> tmp;
    read_volume4D(tmp,opts.inputfname.value()); 
    tmp.setmatrix(newData,Mask);
    save_volume4D(tmp,logger.appendDir(opts.outputfname.value() + "_ICAfiltered")); 
   
    return 0;
  } // int remove_components()

  void MelodicData::create_RXweight()
  {
    message("Reading the weights for the covariance R_X from file "<< opts.segment.value() << endl);
  
    volume4D<float> tmpRX;
    read_volume4D(tmpRX,opts.segment.value());
    RXweight = tmpRX.matrix(Mask);
  } 
 
  void MelodicData::create_mask(volume4D<float> &theData, 
				volume<float> &theMask)
  {
    if(opts.use_mask.value() && opts.maskfname.value().size()>0){   // mask provided 
      read_volume(theMask,opts.maskfname.value());
      message("Mask provided : " << opts.maskfname.value()<<endl);
    }
    else{
      if(opts.perf_bet.value() && opts.use_mask.value()){ //use BET
	message("Create mask ... ");
	// set up all strings
	string BET_outputfname = string(Mean_fname)+"_brain";

	string BET_path = opts.binpath + "bet";
	string BET_optarg = "-m -f 0.4"; // see man bet
	string Mask_fname = BET_outputfname+"_mask";

	// Setup external call to BET:

	ostringstream osc;
	osc  << BET_path << " " << Mean_fname << " " 
	     << BET_outputfname << " " << BET_optarg << " > /dev/null ";
	
        message("  Calling BET: " << osc.str() << endl);
	system(osc.str().c_str());
	
	// read back the Mask file   
	read_volume(theMask,Mask_fname);

	message("done" << endl);
      }  
      else{
	if(opts.use_mask.value()){   //just threshold the Mean
	  message("Create mask ... ");
	  float Mmin, Mmax, Mtmp;
	  Mmin = Mean.min(); Mmax = Mean.max();
	  theMask = binarise(Mean,Mmin + opts.threshold.value()* (Mmax-Mmin),Mmax);
          Mtmp = Mmin + opts.threshold.value()* (Mmax-Mmin);
	  message("done" << endl);
	}
	else{ //well, don't threshold then
	  theMask = Mean;
	  theMask = 1.0;
	}
      }
    }
    if(opts.remove_endslices.value()){ 
      // just in case mc introduced something nasty
      message("  Deleting end slices" << endl);
      for(int ctr1=theMask.miny(); ctr1<=theMask.maxy(); ctr1++){
	for(int ctr2=theMask.minx(); ctr2<=theMask.maxx(); ctr2++){   
	  theMask(ctr2,ctr1,Mask.minz()) = 0.0;
	  theMask(ctr2,ctr1,Mask.maxz()) = 0.0;
	}
      }
    }
  } //void create_mask()

  void MelodicData::sort()
  {
    int numComp = mixMatrix.Ncols(), numVox = IC.Ncols(), 
        numTime = mixMatrix.Nrows(), i,j;

    for(int ctr_i = 1; ctr_i <= numComp; ctr_i++){
      if(IC.Row(ctr_i).Sum()>0){
	flipres(ctr_i); };}
    //    cerr << "HERE2" << endl << endl;


    // re-order wrt standard deviation of IC maps
    message("  sorting IC maps" << endl << endl);  
    Matrix tmpscales, tmpICrow, tmpMIXcol;
    tmpscales = stdev(IC,2);
    ICstats = tmpscales;

    //cerr << "SCLAES2 " << tmpscales << endl;

    double max_val, min_val = tmpscales.Minimum()-1;

    for(int ctr_i = 1; ctr_i <= numComp; ctr_i++){

      max_val = tmpscales.Maximum2(i,j);
      ICstats(ctr_i,1)=max_val;

      //cerr << endl << "scales: " << tmpscales << " ctr_i " << ctr_i << 
      //	" i " << i << " j " << j << " max " << max_val << endl << endl;

  
      tmpICrow = IC.Row(ctr_i);
      tmpMIXcol = mixMatrix.Column(ctr_i);
      
      IC.SubMatrix(ctr_i,ctr_i,1,numVox) = IC.SubMatrix(i,i,1,numVox);
      mixMatrix.SubMatrix(1,numTime,ctr_i,ctr_i) = 
	mixMatrix.SubMatrix(1,numTime,i,i);
  
      IC.SubMatrix(i,i,1,numVox) = tmpICrow.SubMatrix(1,1,1,numVox);
      mixMatrix.SubMatrix(1,numTime,i,i) = tmpMIXcol.SubMatrix(1,numTime,1,1);
  
      tmpscales(i,1)=tmpscales(ctr_i,1);
      tmpscales(ctr_i,1)=min_val;
      }

    //cerr << " ICstats " << ICstats << endl << endl;

    ICstats /= ICstats.Column(1).Sum();
    ICstats *= 100;


    //cerr << " ICstats " << ICstats << endl << endl;
    
    if(EVP.Storage()>0){
      tmpscales = ICstats.Column(1).AsMatrix(ICstats.Nrows(),1) * EVP(1,numComp);
      ICstats |= tmpscales;
    }

    if(DataVN.Storage()>0&&stdDev.Storage()>0){
      //cerr << " ICstats " << ICstats << endl << endl;

      Matrix copeP(tmpscales), copeN(tmpscales);
      Matrix max_ICs(tmpscales), min_ICs(tmpscales);

      for(int ctr_i = 1; ctr_i <= numComp; ctr_i++){
	int i,j;
	max_ICs(ctr_i,1) = IC.Row(ctr_i).Maximum2(i,j);
	//cerr << " ICstats " << ICstats << endl << endl;

	//cerr << endl <<(pinv(mixMatrix)*DataVN.Column(j)) << endl;
	copeP(ctr_i,1) = std::abs((pinv(mixMatrix)*DataVN.Column(j)).Row(ctr_i).AsScalar()*stdDev(1,j)*100*(mixMatrix.Column(ctr_i).Maximum()-mixMatrix.Column(ctr_i).Minimum())/meanR(1,j));

	min_ICs(ctr_i,1) = IC.Row(ctr_i).Minimum2(i,j);
	copeN(ctr_i,1) = -1.0*std::abs((pinv(mixMatrix)*DataVN.Column(j)).Row(ctr_i).AsScalar()*stdDev(1,j)*100*(mixMatrix.Column(ctr_i).Maximum()-mixMatrix.Column(ctr_i).Minimum())/meanR(1,j));

      }
      ICstats |= copeP;
      ICstats |= copeN;
    }
    
    mixFFT=calc_FFT(mixMatrix);
    unmixMatrix = pinv(mixMatrix);

    //if(ICstats.Storage()>0){cout << "ICstats: " << ICstats.Nrows() <<"x" << ICstats.Ncols() << endl;}else{cout << "ICstats empty " <<endl;}
  }


  void MelodicData::status(const string &txt)
  {
    cout << "MelodicData Object " << txt << endl;
    if(Data.Storage()>0){cout << "Data: " << Data.Nrows() <<"x" << Data.Ncols() << endl;}else{cout << "Data empty " <<endl;}
    if(pcaE.Storage()>0){cout << "pcaE: " << pcaE.Nrows() <<"x" << pcaE.Ncols() << endl;}else{cout << "pcaE empty " <<endl;}
    if(pcaD.Storage()>0){cout << "pcaD: " << pcaD.Nrows() <<"x" << pcaD.Ncols() << endl;}else{cout << "pcaD empty " <<endl;}
    if(whiteMatrix.Storage()>0){cout << "white: " << whiteMatrix.Nrows() <<"x" << whiteMatrix.Ncols() << endl;}else{cout << "white empty " <<endl;}
    if(dewhiteMatrix.Storage()>0){cout << "dewhite: " << dewhiteMatrix.Nrows() <<"x" << dewhiteMatrix.Ncols() << endl;}else{cout << "dewhite empty " <<endl;}
    if(mixMatrix.Storage()>0){cout << "mix: " << mixMatrix.Nrows() <<"x" << mixMatrix.Ncols() << endl;}else{cout << "mix empty " <<endl;}
    if(unmixMatrix.Storage()>0){cout << "unmix: " << unmixMatrix.Nrows() <<"x" << unmixMatrix.Ncols() << endl;}else{cout << "unmix empty " <<endl;}
    if(IC.Storage()>0){cout << "IC: " << IC.Nrows() <<"x" << IC.Ncols() << endl;}else{cout << "IC empty " <<endl;}
   
  } //void status()

  Matrix MelodicData::calc_FFT(const Matrix& Mat)
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
	if(opts.logPower.value()) tmpPow = log(tmpPow);
	if(res.Storage()==0){res= tmpPow;}else{res|=tmpPow;}
      }
    return res;
  } //Matrix calc_FFT()

  Matrix MelodicData::smoothColumns(const Matrix &inp)
  {
    Matrix temp(inp);
    int ctr1 = temp.Nrows();
    Matrix temp2(temp);
    temp2=0;
    
    temp = temp.Row(4) & temp.Row(3) & temp.Row(2) & temp & temp.Row(ctr1-1) 
      & temp.Row(ctr1-2) &temp.Row(ctr1-3);
    
    double kern[] ={0.0045 , 0.055, 0.25, 0.4, 0.25, 0.055, 0.0045};
    double fac = 0.9090909;

    //  Matrix FFTinp;
    //FFTinp = calc_FFT(inp);
    //write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_FT"),
    //			 FFTinp);    FFTinp = calc_FFT(inp);
  //write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_Sinp"),
    //		 inp);   
    
    for(int cc=1;cc<=temp2.Ncols();cc++){
      //  double all = FFTinp.Column(cc).Rows(1,30).SumAbsoluteValue() / FFTinp.Column(cc).SumAbsoluteValue();
      // if( all > 0.5){
      for(int cr=1;cr<=temp2.Nrows();cr++){
	temp2(cr,cc) = fac*( kern[0] * temp(cr,cc) + kern[1] * temp(cr+1,cc) + 
			     kern[2] * temp(cr+2,cc) + kern[3] * temp(cr+3,cc) + 
			     kern[4] * temp(cr+4,cc) + kern[5] * temp(cr+5,cc) + 
			     kern[6] * temp(cr+6,cc));
      }//}
      //else{
      //	for(int cr=1;cr<=temp2.Nrows();cr++){
      //	  temp2(cr,cc) = temp(cr+3,cc);
      //	}
      //}
    }
//write_ascii_matrix(logger.appendDir(opts.outputfname.value() + "_Sout"),
//		 temp2);   
    return temp2;
  }
}

