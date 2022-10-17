/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melpca.cc - PCA and whitening 

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
#include "utils/log.h"
#include "meloptions.h"
#include "meldata.h"
#include "melodic.h"
#include "newmatap.h"
#include "newmatio.h"
#include "melpca.h"
#include "libprob.h"
#include <sstream>

using namespace Utilities;
using namespace NEWIMAGE;

namespace Melodic{

   void MelodicPCA::perf_pca(const Matrix &Data)
   {    
     message("Starting PCA  ... ");

     SymmetricMatrix Corr;

     if(opts.segment.value().length()>0){
       Matrix Res;
       Res = ones(Data.Nrows(),1)*melodat.get_RXweight();
       Res = SP(melodat.get_Data(),Res);
       Corr = cov(Res.t());
     }
     else{
       Corr = cov(melodat.get_Data().t());
     }

     Matrix tmpE;
     DiagonalMatrix tmpD;
 
     EigenValues(Corr,tmpD,tmpE);

     if(opts.tsmooth.value()){
       message(" temporal smoothing of Eigenvectors " << endl);
       tmpE=melodat.smoothColumns(tmpE);
     }
  
     melodat.set_pcaE(tmpE);
     melodat.set_pcaD(tmpD);
     
     RowVector AdjEV;
     AdjEV = tmpD.AsRow().Reverse();
     SortDescending(AdjEV);
  
     RowVector PercEV(AdjEV);     
     PercEV = cumsum(AdjEV / sum(AdjEV,2).AsScalar());
     write_ascii_matrix(logger.appendDir("eigenvalues_percent"),PercEV);
     melodat.set_EVP(PercEV);
     AdjEV = (AdjEV - min(AdjEV).AsScalar())/(max(AdjEV).AsScalar() - min(AdjEV).AsScalar());
     melodat.set_EV((AdjEV));
     message("done" << endl);
   }
  
  void MelodicPCA::perf_white(const Matrix &Data)
  {
    Matrix RE;
    DiagonalMatrix RD;
    Matrix tmpWhite;
    Matrix tmpDeWhite;

    int N = melodat.get_pcaE().Ncols();    
    if(opts.pca_dim.value() > N){
      message("dimensionality too large - using -dim " << N << 
              " instead " << endl);
      opts.pca_dim.set_T(N);
    }
    if(opts.pca_dim.value() < 0){
      if(opts.remove_meanvol.value()){
	opts.pca_dim.set_T(N-2);
      }else{
	opts.pca_dim.set_T(N-1);
      }
    }
    if(opts.pca_dim.value() ==0){
      opts.pca_dim.set_T(pcadim());
      if(melodat.get_Data().Nrows()<20){
	opts.pca_dim.set_T(N-2);
	message("too few data points for full estimation, using "
		<< opts.pca_dim.value() << " instead"<< endl);
      }
    }
    if(opts.approach.value()==string("jade") && opts.pca_dim.value() > 30){
      message("dimensionality too large for jade estimation - using --dim 30 instead" << endl);
      opts.pca_dim.set_T(30);
    }
    message("Start whitening using  "<< opts.pca_dim.value()<<" dimensions ... " << endl);
    RowVector tmpEVP;
    tmpEVP << melodat.get_EVP();
    float varp = 1.0;
    if(opts.pca_dim.value() <= tmpEVP.Ncols()){
      varp = tmpEVP(opts.pca_dim.value());
    }
    message("  retaining "<< varp*100 <<" percent of the variability " << endl);

    RE = melodat.get_pcaE().Columns(N-opts.pca_dim.value()+1,N);
    RD << abs(melodat.get_pcaD().SymSubMatrix(N-opts.pca_dim.value()+1,N));    

    tmpWhite = sqrt(abs(RD.i()))*RE.t();
    tmpDeWhite = RE*sqrt(RD);
    melodat.set_white(tmpWhite);
    melodat.set_dewhite(tmpDeWhite);
    message(" ... done"<< endl << endl);
   }

  int MelodicPCA::pcadim()
  { 
    message("Estimating the number of sources from the data (PPCA) ..." << endl);

    //First, estimate the smoothness of the data;
    //   set up all strings

    string SM_path = opts.binpath + "smoothest";
    string Mask_fname = logger.appendDir("mask");

    if(opts.segment.value().length()>0){
      Mask_fname =  opts.segment.value();
    } 

    // Setup external call to smoothest:
    ostringstream osc;
    osc  << SM_path << " -d " << melodat.data_dim()
    	 << " -r " << opts.inputfname.value() << " -m " 
    	 << Mask_fname << " > " << logger.appendDir("smoothest");
    
    message("  Calling Smoothest: " << osc.str() << endl);
    system(osc.str().c_str());

    //read back the results
    ifstream in;
    string str;
    float Resel = 1.0;

    in.open(logger.appendDir("smoothest").c_str(), ios::in);
    if(in>0){
      for(int ctr=1; ctr<7; ctr++){ in >> str;}
      in.close();
      if(str!="nan"){
	Resel = atof(str.c_str());
      }
    }

    //cerr << "  Resels : " << Resel << endl << endl;
    
    melodat.set_resels(Resel);
    
    Matrix PPCAest;
  
    //   if(!opts.varnorm.value()){
    SymmetricMatrix Corr;
    if(opts.segment.value().length()>0){
      Matrix Res;
      Res = ones(melodat.get_DataVN().Nrows(),1)*melodat.get_RXweight();
      Res = SP(melodat.get_DataVN(),Res);
      Corr = cov(Res.t());
     }
    else{
      Corr = cov(melodat.get_DataVN().t());
    }
    DiagonalMatrix tmpD;
    Matrix tmpE;
    EigenValues(Corr,tmpD,tmpE);
    // }

    RowVector AdjEV;
    AdjEV << tmpD.AsRow();
    AdjEV = AdjEV.Columns(3,AdjEV.Ncols());
    AdjEV = AdjEV.Reverse();

    RowVector CircleLaw;
    int NumVox = (int) floor(melodat.data_samples()/(2.5*Resel));


    CircleLaw = Feta(int(AdjEV.Ncols()), NumVox);

    for(int ctr=1;ctr<=CircleLaw.Ncols(); ctr++){
      if(CircleLaw(ctr)<5*10e-10){CircleLaw(ctr) = 5*10e-10;}
    } 

    //write_ascii_matrix(logger.appendDir("tmpA1"),AdjEV);
    //AdjEV = AdjEV.Columns(2,AdjEV.Ncols());
    //write_ascii_matrix(logger.appendDir("tmpA2"),AdjEV);

    //cerr<< AdjEV.Nrows() << " x " << AdjEV.Ncols() << endl;
    
    //cerr<< CircleLaw.Nrows() << " x " << CircleLaw.Ncols() << endl;

    float slope;
    slope = CircleLaw.Columns(int(AdjEV.Ncols()/4),AdjEV.Ncols() - 
			      int(AdjEV.Ncols()/4)).Sum() /  
      AdjEV.Columns(int(AdjEV.Ncols()/4),AdjEV.Ncols() - 
		    int(AdjEV.Ncols()/4)).Sum();

    //CircleLaw = slope * (CircleLaw-1) + 1;

    //    write_ascii_matrix(logger.appendDir("claw"),CircleLaw.Columns(1,AdjEV.Ncols()));

    RowVector PercEV(AdjEV);
    PercEV = cumsum(AdjEV / sum(AdjEV,2).AsScalar());

    //    write_ascii_matrix(logger.appendDir("ev"),AdjEV);

    //cerr << int(AdjEV.Ncols()) << "   " << NumVox << "   " << slope << endl;
    AdjEV << SP(AdjEV,pow(CircleLaw.Columns(1,AdjEV.Ncols()),-1));
    //  cerr << "recalculated" << endl;

    SortDescending(AdjEV);
    int maxEV = 1;
    float threshold = 0.98;
    for(int ctr_i = 1; ctr_i < PercEV.Ncols(); ctr_i++){ 
      if((PercEV(ctr_i)<threshold)&&(PercEV(ctr_i+1)>=threshold)){maxEV=ctr_i;}
    }

    if(maxEV<3){maxEV=PercEV.Ncols()/2;}
    RowVector NewEV;
    Matrix temp1;
    temp1 = abs(AdjEV);
    NewEV << temp1.SubMatrix(1,1,1,maxEV);
    PPCAest = ppca_est(NewEV, NumVox);
    RowVector estimators(5);
    estimators = 1.0;
    
    Matrix PPCA2(PPCAest);
    
    for(int ctr=1; ctr<=PPCA2.Ncols(); ctr++){
      PPCA2.Column(ctr) = (PPCA2.Column(ctr) - 
			   min(PPCA2.Column(ctr)).AsScalar()) / 
	( max(PPCA2.Column(ctr)).AsScalar() - 
	  min(PPCA2.Column(ctr)).AsScalar());
    }
    
    int ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,2) < PPCAest(ctr_i+1,2))&&(ctr_i<maxEV))
      {estimators(1)=ctr_i+1;ctr_i++;}
    ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,3) < PPCAest(ctr_i+1,3))&&(ctr_i<maxEV))
      {estimators(2)=ctr_i+1;ctr_i++;}
    ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,4) < PPCAest(ctr_i+1,4))&&(ctr_i<maxEV))
      {estimators(3)=ctr_i+1;ctr_i++;}
    ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,5) < PPCAest(ctr_i+1,5))&&(ctr_i<maxEV))
      {estimators(4)=ctr_i+1;ctr_i++;}
    ctr_i = 1;
    while((ctr_i< PPCAest.Nrows()-1)&&
	  (PPCAest(ctr_i,6) < PPCAest(ctr_i+1,6))&&(ctr_i<maxEV))
      {estimators(5)=ctr_i+1;ctr_i++;}

    int res = 0;
    ColumnVector PPCA;

    if(opts.pca_est.value() == string("lap")){
      res = int(estimators(1));
      PPCA << PPCA2.Column(2);
    }

    if(opts.pca_est.value() == string("bic")){
      res = int(estimators(2));
      PPCA << PPCA2.Column(2);
    }
    if(opts.pca_est.value() == string("mdl")){
      res = int(estimators(3));
      PPCA << PPCA2.Column(4);
    }
    if(opts.pca_est.value() == string("aic")){
      res = int(estimators(5));
      PPCA << PPCA2.Column(6);
    }
    if(res==0){//median estimator
      PPCA = PPCA2.Column(2);

      for(int ctr=1; ctr<=PPCA2.Nrows(); ctr++){ 
	RowVector tmp = PPCA2.SubMatrix(ctr,ctr,2,6);
// 	SortAscending(tmp);
// 	float themean = float(tmp.Sum()/5);
// 	if(std::abs(int(tmp(2)-themean)) < std::abs(int(tmp(3)-themean)))
// 	  PPCA(ctr) = tmp(2);
// 	else
// 	  PPCA(ctr) = tmp(3);

	PPCA(ctr) = float(tmp.Sum()/5);
      }

      ctr_i = 1; 
      while((PPCA(ctr_i) < PPCA(ctr_i+1))&&(ctr_i<maxEV)){
	res=ctr_i+1;ctr_i++;
      }
    }

    //    cerr << estimators << "   "  << res << endl;
			       
    //write_ascii_matrix(logger.appendDir("PPCA2"),PPCA2);
    AdjEV = (AdjEV - min(AdjEV).AsScalar())/(max(AdjEV).AsScalar() - min(AdjEV).AsScalar());


    write_ascii_matrix(logger.appendDir("PPCA"),PPCA);			      
    write_ascii_matrix(logger.appendDir("eigenvalues_adjusted"),AdjEV.t());
    write_ascii_matrix(logger.appendDir("eigenvalues_percent"),PercEV.t());

    melodat.set_EVP(PercEV);
    melodat.set_EV(AdjEV);
    melodat.set_PPCA(PPCA);

    //PPCA << sum(PPCAest.Columns(2,6),2);
    
    //ctr_i = 1;
    
    //while((PPCA(ctr_i) < PPCA(ctr_i+1))&&(ctr_i<maxEV)){
    //  res=ctr_i+1;ctr_i++;
    //} 
    
    //res = int(sum(estimators,2).AsScalar()/5);
    
    // res = int(estimators(1)); // Laplace approximation
    //     SortAscending(estimators);
    //     if(std::abs(int(estimators(2))-res) < std::abs(int(estimators(3))-res))
    //       res = int(estimators(2));
    //     else
    //       res = int(estimators(3));
    
    
    //write_ascii_matrix(logger.appendDir("PPCA"),PPCAest);
    //write_ascii_matrix(logger.appendDir("dimest"),estimators);
    
    return res;
  }

RowVector MelodicPCA::Feta(int n1, int n2)
  {
    float nu = (float) n1/n2; 
    float bm = pow((1-sqrt(nu)),2.0);
    float bp = pow((1+sqrt(nu)),2.0);

    //cerr << "nu, bm, bp " << nu << " " <<bm << " " << bp << endl;

    float lrange = 0.9*bm;
    float urange = 1.1*bp;

    // int dummy;
   
    RowVector eta(30*n1);
    float rangestepsize = (urange - lrange) / eta.Ncols(); 
    for(int ctr_i = 0; ctr_i < eta.Ncols(); ctr_i++){ 
      eta(ctr_i+1) = lrange + rangestepsize * (ctr_i);
    }

    RowVector teta(10*n1);
    teta = 0;
    float stepsize = (bp - bm) / teta.Ncols();
    for(int ctr_i = 0; ctr_i < teta.Ncols(); ctr_i++){ 
      teta(ctr_i+1) = stepsize*(ctr_i);
    }  
    
    //cerr << teta(1)+bm << " " << teta(1000)+bm << endl;
    //cerr << eta(1)<< " " << eta(eta.Ncols())<< endl;
    
    //cerr << "BP1" << endl;
 
    //write_ascii_matrix(logger.appendDir("teta"),teta.t());


    //  RowVector tmp1(teta);
    //     tmp1 = teta + bm;
    
    //     cerr << "tmp1" << endl;
    //     tmp1 =  pow(2*M_PI*nu*(tmp1),-1);
    //     cerr << "tmp1" << endl;
    
    //     RowVector tmp2(teta); 
    //     cerr << "tmp2" << endl;
    //     tmp2 = SP(teta, bp-bm-teta);
    //     cerr << "tmp2" << endl;
    //     tmp2=abs(tmp2);
    //     cerr << "tmp2" << endl;


    RowVector feta(teta);
    feta = SP(pow(2*M_PI*nu*(teta + bm),-1), pow(SP(teta, bp-bm-teta),0.5));
   
    //Matrix location;
    teta = teta + bm;

    //cerr << "teta : " << teta.Nrows() << " x " << teta.Ncols() << endl;
    //cerr << "eta : " << eta.Nrows() << " x " << eta.Ncols() << endl;
    //cerr << "feta : " << feta.Nrows() << " x " << feta.Ncols() << endl;
    //c/err << "vor location (input)" << endl;
    //cin >> dummy;
    //location = SP(teta.t()*ones(1,eta.Ncols()),pow(ones(teta.Ncols(),1)*eta,-1));
    //cerr << "nach location (input)" << endl;
    //cin >> dummy;
    //cerr << " weiter " << endl;

    //for(int ctr_i = 1; ctr_i <= location.Ncols(); ctr_i++){
    //  for(int ctr_j = 1; ctr_j <= location.Nrows(); ctr_j++){
    //	if(location(ctr_j,ctr_i)<1){location(ctr_j,ctr_i)=1;}
    //	else {location(ctr_j,ctr_i)=0;}
    //  }
    // } 
    //write_ascii_matrix(logger.appendDir("location"),location);
    // write_ascii_matrix(logger.appendDir("teta"),teta);
    //write_ascii_matrix(logger.appendDir("eta"),eta);
    //write_ascii_matrix(logger.appendDir("feta"),feta);
 
  
    //RowVector claw; 
    //  claw = n1*(1-sum(SP(stepsize*feta.t()*ones(1,eta.Ncols()),location),1).AsRow());

    RowVector claw(eta);
    claw = 0;
    for(int ctr_i = 1; ctr_i <= eta.Ncols(); ctr_i++){
      double tmpval = 0.0;
      for(int ctr_j = 1; ctr_j <= teta.Ncols(); ctr_j++){
	if(( double(teta(ctr_j))/double(eta(ctr_i)) )<1)
	  tmpval += feta(ctr_j);
      }
      claw(ctr_i) = n1*(1-stepsize*tmpval);
    }
    
    //write_ascii_matrix(logger.appendDir("claw"),claw);
    //cerr << "BP1" << endl;
    RowVector Res(n1); //invert the CDF
    //cerr << "n1=" << n1 << endl;
    for(int ctr_i = 1; ctr_i < eta.Ncols(); ctr_i++){
      if(floor(claw(ctr_i))>floor(claw(ctr_i+1))){
	//	cerr << int(floor(claw(ctr_i))) << " ";
	Res(int(floor(claw(ctr_i)))) = eta(ctr_i);
      }
    }
 
    //cerr << endl;
    //    cerr << " Done with loop " << int(floor(tmp4b(ctr_i))) << endl; 
    //write_ascii_matrix(logger.appendDir("claw-dstn"),Res);
    return Res;
  }


  RowVector MelodicPCA::cumsum(const RowVector& Inp)
  {
    UpperTriangularMatrix UT(Inp.Ncols());
    UT=1.0;
    RowVector Res;
    Res = Inp * UT;
    return Res;
  }


  Matrix MelodicPCA::ppca_est(const RowVector& eigenvalues, const int N)
  {
    RowVector logLambda(eigenvalues);
    logLambda = log(eigenvalues);

    int d = logLambda.Ncols();

    RowVector k(d);
    for(int ctr_i = 1; ctr_i <=d; ctr_i++){
      k(ctr_i)=ctr_i;
    }
   
    RowVector m(d);
    m=d*k-0.5*SP(k,k+1); 

    RowVector loggam(d);
    loggam = 0.5*k.Reverse();
    for(int ctr_i = 1; ctr_i <=d; ctr_i++){
      loggam(ctr_i)=lgam(loggam(ctr_i));
    }
    loggam = cumsum(loggam); 

    RowVector l_probU(d);
    l_probU = -log(2)*k + loggam - cumsum(0.5*log(M_PI)*k.Reverse());

    RowVector tmp1;
    tmp1 = -cumsum(eigenvalues).Reverse()+sum(eigenvalues,2).AsScalar();
    tmp1(1) = 0.95*tmp1(2);
    tmp1=tmp1.Reverse();

    RowVector tmp2;
    tmp2 = -cumsum(logLambda).Reverse()+sum(logLambda,2).AsScalar();
    tmp2(1)=tmp2(2);
    tmp2=tmp2.Reverse();

    RowVector tmp3;
    tmp3 = d-k;
    tmp3(d) = 1.0;

    RowVector tmp4;
    tmp4 = SP(tmp1,pow(tmp3,-1));    
    for(int ctr_i = 1; ctr_i <=d; ctr_i++){
      if(tmp4(ctr_i)<0.01){tmp4(ctr_i)=0.01;}
      if(tmp3(ctr_i)<0.01){tmp3(ctr_i)=0.01;}
      if(tmp1(ctr_i)<0.01){tmp1(ctr_i)=0.01;}
    }

    RowVector l_nu;
    l_nu = SP(-N/2*(d-k),log(tmp4));
    l_nu(d) = 0;

    RowVector l_lam;
    l_lam = -(N/2)*cumsum(logLambda);

    RowVector l_lhood;
    l_lhood = SP(pow(tmp3,-1),tmp2)-log(SP(pow(tmp3,-1),tmp1));

    Matrix t1,t2, t3;
    UpperTriangularMatrix triu(d);
    triu = 1.0;
    for(int ctr_i = 1; ctr_i <= triu.Ncols(); ctr_i++){
      triu(ctr_i,ctr_i)=0.0;
    }
    t1 = (ones(d,1) * eigenvalues);
    t1 = SP(triu,t1.t() - t1);
    t2 = pow(tmp4,-1).t()*ones(1,d);
    t3 = ones(d,1)*pow(eigenvalues,-1);
    t2 = SP(triu, t2.t()-t3.t());
    for(int ctr_i = 1; ctr_i <= t1.Ncols(); ctr_i++){
      for(int ctr_j = 1; ctr_j <= t1.Nrows(); ctr_j++){
	if(t1(ctr_j,ctr_i)<=0){t1(ctr_j,ctr_i)=1;}
      } 
    }
    for(int ctr_i = 1; ctr_i <= t2.Ncols(); ctr_i++){
      for(int ctr_j = 1; ctr_j <= t2.Nrows(); ctr_j++){
	if(t2(ctr_j,ctr_i)<=0){t2(ctr_j,ctr_i)=1;}
      }
    } 
    t1 = cumsum(sum(log(t1),2).AsRow());
    t2 = cumsum(sum(log(t2),2).AsRow());

    RowVector l_Az(d);
    l_Az << (t1+t2);

    RowVector l_lap;
    l_lap = l_probU + l_nu +l_Az + l_lam + 0.5*log(2*M_PI)*(m+k)-0.5*log(N)*k;
 
    RowVector l_BIC;
    l_BIC = l_lam + l_nu - 0.5*log(N)*(m+k);

    RowVector l_RRN;
    l_RRN = -0.5*N*SP(k,log(SP(cumsum(eigenvalues),pow(k,-1))))+l_nu;

    RowVector l_AIC;
    l_AIC = -2*N*SP(tmp3,l_lhood)+ 2*(1+d*k+0.5*(k-1));
    l_AIC = -l_AIC;

    RowVector l_MDL;
    l_MDL = -N*SP(tmp3,l_lhood)+ 0.5*(1+d*k+0.5*(k-1))*log(N);
    l_MDL = -l_MDL;

    Matrix Res;

    Res = eigenvalues.t();
    Res |= l_lap.t();
    Res |= l_BIC.t();
    Res |= l_MDL.t();
    Res |= l_RRN.t();
    Res |= l_AIC.t();
    
   
    return Res;

  }




}


