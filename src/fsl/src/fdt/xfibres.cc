/* Xfibres Diffusion Partial Volume Model  

    Tim Behrens - FMRIB Image Analysis Group
 
    Copyright (C) 2005 University of Oxford  */

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
#include <strstream>
#define WANT_STREAM
#define WANT_MATH
//  #include "newmatap.h"
//  #include "newmatio.h"
#include <string>
#include <math.h>
#include "utils/log.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "stdlib.h"
#include "fibre.h"
#include "xfibresoptions.h"

//#include "bint/model.h"
//#include "bint/lsmcmcmanager.h"
//#include "bint/lslaplacemanager.h"

using namespace FIBRE;
using namespace  Xfibres;
using namespace Utilities;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;




inline float min(float a,float b){
  return a<b ? a:b;}
inline float max(float a,float b){
  return a>b ? a:b;}
inline Matrix Anis()
{ 
  Matrix A(3,3);
  A << 1 << 0 << 0
    << 0 << 0 << 0
    << 0 << 0 << 0;
  return A;
}

inline Matrix Is()
{ 
  Matrix I(3,3);
  I << 1 << 0 << 0
    << 0 << 1 << 0
    << 0 << 0 << 1;
  return I;
}

inline ColumnVector Cross(const ColumnVector& A,const ColumnVector& B)
{
  ColumnVector res(3);
  res << A(2)*B(3)-A(3)*B(2)
      << A(3)*B(1)-A(1)*B(3)
      << A(1)*B(2)-B(1)*A(2);
  return res;
}

inline Matrix Cross(const Matrix& A,const Matrix& B)
{
  Matrix res(3,1);
  res << A(2,1)*B(3,1)-A(3,1)*B(2,1)
      << A(3,1)*B(1,1)-A(1,1)*B(3,1)
      << A(1,1)*B(2,1)-B(1,1)*A(2,1);
  return res;
}

float mod(float a, float b){
  while(a>b){a=a-b;}
  while(a<0){a=a+b;} 
  return a;
}


Matrix form_Amat(const Matrix& r,const Matrix& b)
{
  Matrix A(r.Ncols(),7);
  Matrix tmpvec(3,1), tmpmat;
  
  for( int i = 1; i <= r.Ncols(); i++){
    tmpvec << r(1,i) << r(2,i) << r(3,i);
    tmpmat = tmpvec*tmpvec.t()*b(1,i);
    A(i,1) = tmpmat(1,1);
    A(i,2) = 2*tmpmat(1,2);
    A(i,3) = 2*tmpmat(1,3);
    A(i,4) = tmpmat(2,2);
    A(i,5) = 2*tmpmat(2,3);
    A(i,6) = tmpmat(3,3);
    A(i,7) = 1;
  }
  return A;
}

inline SymmetricMatrix vec2tens(ColumnVector& Vec){
  SymmetricMatrix tens(3);
  tens(1,1)=Vec(1);
  tens(2,1)=Vec(2);
  tens(3,1)=Vec(3);
  tens(2,2)=Vec(4);
  tens(3,2)=Vec(5);
  tens(3,3)=Vec(6);
  return tens;
}



class Samples{
  xfibresOptions& opts;
  Matrix m_dsamples;
  Matrix m_S0samples;
  Matrix m_lik_energy;
  
  vector<Matrix> m_thsamples;
  vector<Matrix> m_phsamples;
  vector<Matrix> m_fsamples;
  vector<Matrix> m_lamsamples;

  //for storing means
  RowVector m_mean_dsamples;
  RowVector m_mean_S0samples;
  vector<Matrix> m_dyadic_vectors;
  vector<RowVector> m_mean_fsamples;
  vector<RowVector> m_mean_lamsamples;
  
  float m_sum_d;
  float m_sum_S0;
  vector<SymmetricMatrix> m_dyad;
  vector<float> m_sum_f;
  vector<float> m_sum_lam;
  int m_nsamps;
  ColumnVector m_vec;
  
  volume<int> m_vol2matrixkey;
  Matrix m_matrix2volkey;
  volume<int> m_beenhere;

  
public:

  Samples( volume<int> vol2matrixkey,Matrix matrix2volkey,int nvoxels):
    opts(xfibresOptions::getInstance()),m_vol2matrixkey(vol2matrixkey),m_matrix2volkey(matrix2volkey){
    
    m_beenhere=m_vol2matrixkey*0;
    int count=0;
    int nsamples=0;
    
    for(int i=0;i<opts.njumps.value();i++){
      count++;
      if(count==opts.sampleevery.value()){
	count=0;nsamples++;
      }
    }
 

    m_dsamples.ReSize(nsamples,nvoxels);
    m_dsamples=0;
    m_S0samples.ReSize(nsamples,nvoxels);
    m_S0samples=0;
    m_lik_energy.ReSize(nsamples,nvoxels);
    
    m_mean_dsamples.ReSize(nvoxels);
    m_mean_dsamples=0;
    m_mean_S0samples.ReSize(nvoxels);
    m_mean_S0samples=0;
    Matrix tmpvecs(3,nvoxels);
    tmpvecs=0;
    m_sum_d=0;
    m_sum_S0=0;
    SymmetricMatrix tmpdyad(3);
    tmpdyad=0;
    m_nsamps=nsamples;
    m_vec.ReSize(3);
    for(int f=0;f<opts.nfibres.value();f++){
      m_thsamples.push_back(m_S0samples);
      m_phsamples.push_back(m_S0samples);
      m_fsamples.push_back(m_S0samples);
      m_lamsamples.push_back(m_S0samples);
      

      m_dyadic_vectors.push_back(tmpvecs);
      m_mean_fsamples.push_back(m_mean_S0samples);
      m_mean_lamsamples.push_back(m_mean_S0samples);

      m_sum_lam.push_back(0);
      m_sum_f.push_back(0);
      m_dyad.push_back(tmpdyad);
    }
 
  }
  
  
  void record(Multifibre& mfib, int vox, int samp){
    m_dsamples(samp,vox)=mfib.get_d();
    m_sum_d+=mfib.get_d();
    m_S0samples(samp,vox)=mfib.get_S0();
    m_sum_S0+=mfib.get_S0();
    m_lik_energy(samp,vox)=mfib.get_likelihood_energy();
    for(int f=0;f<opts.nfibres.value();f++){
      float th=mfib.fibres()[f].get_th();
      float ph=mfib.fibres()[f].get_ph();
      m_thsamples[f](samp,vox)=th;
      m_phsamples[f](samp,vox)=ph;
      m_fsamples[f](samp,vox)=mfib.fibres()[f].get_f();
      m_lamsamples[f](samp,vox)=mfib.fibres()[f].get_lam();
      //for means
      m_vec << sin(th)*cos(ph) << sin(th)*sin(ph)<<cos(th) ;
      m_dyad[f] << m_dyad[f]+m_vec*m_vec.t();
      m_sum_f[f]+=mfib.fibres()[f].get_f();
      m_sum_lam[f]+=mfib.fibres()[f].get_lam();
    }
  }
  
  void finish_voxel(int vox){
    m_mean_dsamples(vox)=m_sum_d/m_nsamps;
    m_mean_S0samples(vox)=m_sum_S0/m_nsamps;

    m_sum_d=0;
    m_sum_S0=0;

    DiagonalMatrix dyad_D; //eigenvalues
    Matrix dyad_V; //eigenvectors
    int nfibs=0;
    for(int f=0;f<opts.nfibres.value();f++){
      
      EigenValues(m_dyad[f],dyad_D,dyad_V);
      int maxeig;
      if(dyad_D(1)>dyad_D(2)){
	if(dyad_D(1)>dyad_D(3)) maxeig=1;
	else maxeig=3;
      }
      else{
	if(dyad_D(2)>dyad_D(3)) maxeig=2;
	else maxeig=3;
      }
      m_dyadic_vectors[f](1,vox)=dyad_V(1,maxeig);
      m_dyadic_vectors[f](2,vox)=dyad_V(2,maxeig);
      m_dyadic_vectors[f](3,vox)=dyad_V(3,maxeig);
      
      if((m_sum_f[f]/m_nsamps)>0.05){
	nfibs++;
      }
      m_mean_fsamples[f](vox)=m_sum_f[f]/m_nsamps;
      m_mean_lamsamples[f](vox)=m_sum_lam[f]/m_nsamps;
      
      m_dyad[f]=0;
      m_sum_f[f]=0;
      m_sum_lam[f]=0;
    }
    m_beenhere(int(m_matrix2volkey(vox,1)),int(m_matrix2volkey(vox,2)),int(m_matrix2volkey(vox,3)))=nfibs;
    cout<<nfibs<<endl;
    cout <<"boobooboo"<<endl;
    cout<<int(m_matrix2volkey(vox,1))<<" "<<int(m_matrix2volkey(vox,2))<<" "<<int(m_matrix2volkey(vox,3))<<endl;
}
  
  
  bool neighbour_initialise(int vox, Multifibre& mfibre){
    int xx = int(m_matrix2volkey(vox,1));
    int yy = int(m_matrix2volkey(vox,2));
    int zz = int(m_matrix2volkey(vox,3));
    int voxn=1,voxbest=1;
    bool ret=false;
    int maxfib=1;
    float maxsf=0;
    for(int x=xx-1;x<=xx+1;x++){
      for(int y=yy-1;y<=yy+1;y++){
	for(int z=zz-1;z<=zz+1;z++){
	  if(m_beenhere(x,y,z)>=maxfib){
	    float sumf=0;
	    voxn=m_vol2matrixkey(x,y,z);
	    for(unsigned int fib=0;fib<m_mean_fsamples.size();fib++){
	      if(voxn!=0)
		sumf+=m_mean_fsamples[fib](voxn);
	      else sumf=0;
	     	    }
	    if(sumf>maxsf){
	      maxsf=sumf;
	      maxfib=m_beenhere(x,y,z);
	      voxbest=voxn;
	      
	      ret=true;
	    }
	  }
	
	  
	}
      } 
    }
    ret=(maxfib>1); //
    cout<<ret<<" ret"<<endl;
    if(ret){
      mfibre.set_d(m_mean_dsamples(voxbest));
      mfibre.set_S0(m_mean_S0samples(voxbest));
      for(int f=0; f<opts.nfibres.value(); f++){
	
	float th;
	float ph;
	ColumnVector vec(3);
	vec(1)= m_dyadic_vectors[f](1,voxbest);
	vec(2)= m_dyadic_vectors[f](2,voxbest);
	vec(3)= m_dyadic_vectors[f](3,voxbest);
	cart2sph(vec,th,ph);
	if(f==0)
	  mfibre.addfibre(th,ph,m_mean_fsamples[f](voxbest),1,false);//no a.r.d. on first fibre
	else
	  mfibre.addfibre(th,ph,m_mean_fsamples[f](voxbest),1,true);

      }
	
      
    }
    return ret;
    
  }
  
  
  void save(const volume<float>& mask){
    volume4D<float> tmp;
    //So that I can sort the output fibres into
    // files ordered by fibre fractional volume..
    vector<Matrix> thsamples_out=m_thsamples;
    vector<Matrix> phsamples_out=m_phsamples;
    vector<Matrix> fsamples_out=m_fsamples;
    vector<Matrix> lamsamples_out=m_lamsamples;
    
    vector<Matrix> dyadic_vectors_out=m_dyadic_vectors;
    vector<Matrix> mean_fsamples_out;
    for(unsigned int f=0;f<m_mean_fsamples.size();f++)
      mean_fsamples_out.push_back(m_mean_fsamples[f]);

    Log& logger = LogSingleton::getInstance();
    tmp.setmatrix(m_dsamples,mask);
    save_volume4D(tmp,logger.appendDir("dsamples"));
    tmp.setmatrix(m_S0samples,mask);
    save_volume4D(tmp,logger.appendDir("S0samples"));
    tmp.setmatrix(m_lik_energy,mask);
    save_volume4D(tmp,logger.appendDir("lik_energy"));
    
    //Sort the output based on mean_fsamples
    // 
    vector<Matrix> sumf;
    for(int f=0;f<opts.nfibres.value();f++){
      Matrix tmp=sum(m_fsamples[f],1);
      sumf.push_back(tmp);
    }  
    for(int vox=1;vox<=m_dsamples.Ncols();vox++){
      vector<pair<float,int> > sfs;
      pair<float,int> ftmp;
      
      for(int f=0;f<opts.nfibres.value();f++){
	ftmp.first=sumf[f](1,vox);
	ftmp.second=f;
	sfs.push_back(ftmp);
      }
      sort(sfs.begin(),sfs.end());
      
      for(int samp=1;samp<=m_dsamples.Nrows();samp++){
	for(int f=0;f<opts.nfibres.value();f++){;
	  thsamples_out[f](samp,vox)=m_thsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	  phsamples_out[f](samp,vox)=m_phsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	  fsamples_out[f](samp,vox)=m_fsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	  lamsamples_out[f](samp,vox)=m_lamsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	}
      }
      
      for(int f=0;f<opts.nfibres.value();f++){
	mean_fsamples_out[f](1,vox)=m_mean_fsamples[sfs[(sfs.size()-1)-f].second](vox);
	dyadic_vectors_out[f](1,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](1,vox);
	dyadic_vectors_out[f](2,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](2,vox);
	dyadic_vectors_out[f](3,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](3,vox);
      }
      
    }
    // save the sorted fibres
    for(int f=0;f<opts.nfibres.value();f++){
      //      element_mod_n(thsamples_out[f],M_PI);
      //      element_mod_n(phsamples_out[f],2*M_PI);
      tmp.setmatrix(thsamples_out[f],mask);
      string oname="th"+num2str(f+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
      tmp.setmatrix(phsamples_out[f],mask);
      oname="ph"+num2str(f+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
      tmp.setmatrix(fsamples_out[f],mask);
      oname="f"+num2str(f+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
      //      tmp.setmatrix(lamsamples_out[f],mask);
      //      oname="lam"+num2str(f+1)+"samples";
      //      save_volume4D(tmp,logger.appendDir(oname));
      tmp.setmatrix(mean_fsamples_out[f],mask);
      oname="mean_f"+num2str(f+1)+"samples";
      save_volume(tmp[0],logger.appendDir(oname));
      tmp.setmatrix(dyadic_vectors_out[f],mask);
      oname="dyads"+num2str(f+1);
      save_volume4D(tmp,logger.appendDir(oname));
    }
  }
  
};









class xfibresVoxelManager{
 
  xfibresOptions& opts;
  
  Samples& m_samples;
  int m_voxelnumber;
  const ColumnVector m_data;
  const ColumnVector& m_alpha;
  const ColumnVector& m_beta;
  const Matrix& m_bvals; 
  Multifibre m_multifibre;
 public:
  xfibresVoxelManager(const ColumnVector& data,const ColumnVector& alpha, 
		      const ColumnVector& beta, const Matrix& b,
		      Samples& samples,int voxelnumber):
    opts(xfibresOptions::getInstance()), 
    m_samples(samples),m_voxelnumber(voxelnumber),m_data(data), 
    m_alpha(alpha), m_beta(beta), m_bvals(b), 
    m_multifibre(m_data,m_alpha,m_beta,m_bvals,opts.nfibres.value(),opts.fudge.value()){ }
  
   
  void initialise(const Matrix& Amat){
    if(!m_samples.neighbour_initialise(m_voxelnumber,m_multifibre)){
      initialise_tensor(Amat);
    }
    m_multifibre.initialise_energies();
    m_multifibre.initialise_props();
  }
  
  void initialise_tensor(const Matrix& Amat){
    //initialising 
    ColumnVector logS(m_data.Nrows()),tmp(m_data.Nrows()),Dvec(7),dir(3);
    SymmetricMatrix tens;   
    DiagonalMatrix Dd;  
    Matrix Vd;  
    float mDd,fsquared;
    float th,ph,f,D,S0;
    for ( int i = 1; i <= logS.Nrows(); i++)
      {
	if(m_data(i)>0){
	  logS(i)=log(m_data(i));
	}
	else{
	  logS(i)=0;
	}
      }
    Dvec = -pinv(Amat)*logS;
    if(  Dvec(7) >  -maxlogfloat  ){ 
      S0=exp(-Dvec(7));
    }
    else{
      S0=m_data.MaximumAbsoluteValue();
    }

    for ( int i = 1; i <= logS.Nrows(); i++)
      {
	if(S0<m_data.Sum()/m_data.Nrows()){ S0=m_data.MaximumAbsoluteValue();  }
	logS(i)=(m_data(i)/S0)>0.01 ? log(m_data(i)):log(0.01*S0);
      }

    Dvec = -pinv(Amat)*logS;
    S0=exp(-Dvec(7));
    if(S0<m_data.Sum()/m_data.Nrows()){ S0=m_data.Sum()/m_data.Nrows();  }
    tens = vec2tens(Dvec);
    EigenValues(tens,Dd,Vd);
    mDd = Dd.Sum()/Dd.Nrows();
    int maxind = Dd(1) > Dd(2) ? 1:2;   //finding maximum eigenvalue
    maxind = Dd(maxind) > Dd(3) ? maxind:3;
    dir << Vd(1,maxind) << Vd(2,maxind) << Vd(3,maxind);
    cart2sph(dir,th,ph);
    th= mod(th,M_PI);
    ph= mod(ph,2*M_PI);
    D = Dd(maxind);

    float numer=1.5*((Dd(1)-mDd)*(Dd(1)-mDd)+(Dd(2)-mDd)*(Dd(2)-mDd)+(Dd(3)-mDd)*(Dd(3)-mDd));
    float denom=(Dd(1)*Dd(1)+Dd(2)*Dd(2)+Dd(3)*Dd(3));
    if(denom>0) fsquared=numer/denom;
    else fsquared=0;
    if(fsquared>0){f=sqrt(fsquared);}
    else{f=0;}
    if(f>=0.95) f=0.95;
    if(f<=0.001) f=0.001;
    if(D<=0) D=2e-3;
    m_multifibre.set_d(D);
    m_multifibre.set_S0(S0);
    if(opts.nfibres.value()>0){
      m_multifibre.addfibre(th,ph,f,1,false);//no a.r.d. on first fibre
      for(int i=2; i<=opts.nfibres.value(); i++){
	 m_multifibre.addfibre();
      }
    
    }
    
    
  }
 

  void runmcmc(){
    int count=0, recordcount=0,sample=1;//sample will index a newmat matrix 
    for( int i =0;i<opts.nburn.value();i++){
      m_multifibre.jump( !opts.no_ard.value() );
      count++;
      if(count==opts.updateproposalevery.value()){
	m_multifibre.update_proposals();
	count=0;
      }
    }
    
    for( int i =0;i<opts.njumps.value();i++){
      m_multifibre.jump(!opts.no_ard.value());
      count++;
     
      if(opts.verbose.value()) 
	{
	  cout<<endl<<i<<" "<<endl<<endl;
	  m_multifibre.report();
	  
	}
      recordcount++;
      if(recordcount==opts.sampleevery.value()){
	m_samples.record(m_multifibre,m_voxelnumber,sample);
	sample++;
	recordcount=0;
      }
      if(count==opts.updateproposalevery.value()){
	m_multifibre.update_proposals();
	count=0;
	
      }
    }
    
    m_samples.finish_voxel(m_voxelnumber);
  }
    
};


  
int main(int argc, char *argv[])
{
  try{  
    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    xfibresOptions& opts = xfibresOptions::getInstance();
    opts.parse_command_line(argc,argv,logger);
    srand(xfibresOptions::getInstance().seed.value());
    Matrix datam, bvals,bvecs,matrix2volkey;
    volume<float> mask;
    volume<int> vol2matrixkey;
    bvals=read_ascii_matrix(opts.bvalsfile.value());
    bvecs=read_ascii_matrix(opts.bvecsfile.value());

    {//scope in which the data exists in 4D format;
      volume4D<float> data;
      read_volume4D(data,opts.datafile.value());
      read_volume(mask,opts.maskfile.value());
      datam=data.matrix(mask); 
      matrix2volkey=data.matrix2volkey(mask);
      vol2matrixkey=data.vol2matrixkey(mask);
    }
    Matrix Amat;
    ColumnVector alpha, beta;
    Amat=form_Amat(bvecs,bvals);
    cart2sph(bvecs,alpha,beta);
    Samples samples(vol2matrixkey,matrix2volkey,datam.Ncols());
    
    for(int vox=1;vox<=datam.Ncols();vox++){
      cout <<vox<<"/"<<datam.Ncols()<<endl;
      xfibresVoxelManager  vm(datam.Column(vox),alpha,beta,bvals,samples,vox);
      vm.initialise(Amat);
      vm.runmcmc();
    }
    
    samples.save(mask);

  }
  catch(Exception& e) 
    {
      cerr << endl << e.what() << endl;
    }
  catch(X_OptionError& e) 
    {
      cerr << endl << e.what() << endl;
    }

  return 0;
}
