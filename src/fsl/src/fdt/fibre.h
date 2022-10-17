/*  Copyright (C) 2005 University of Oxford  */

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

#ifndef __FIBRE_H_
#define __FIBRE_H_


#include <iostream>
#include "stdlib.h"
#include "libprob.h"
#include <cmath>
#include "miscmaths/miscprob.h"

using namespace std; 
using namespace NEWMAT;
using namespace MISCMATHS;

const float maxfloat=1e10;
const float minfloat=1e-10;
const float maxlogfloat=23;
const float minlogfloat=-23;
namespace FIBRE{
  

  class Fibre{
    float m_th;
    float m_ph;
    float m_f;
    float m_lam;
    float m_th_prop;
    float m_ph_prop;
    float m_f_prop;
    float m_lam_prop;
    float m_th_old;
    float m_ph_old;
    float m_f_old;
    float m_lam_old;
    float m_th_prior;
    float m_ph_prior;
    float m_f_prior;
    float m_lam_prior;
    float m_th_old_prior;
    float m_ph_old_prior;
    float m_f_old_prior;
    float m_lam_old_prior;
    float m_prior_en;
    float m_old_prior_en;
    float m_ardfudge;
    int m_th_acc; 
    int m_th_rej;
    int m_ph_acc;
    int m_ph_rej; 
    int m_f_acc;
    int m_f_rej;
    int m_lam_acc; 
    int m_lam_rej;
    bool m_lam_jump;
    ColumnVector m_Signal; 
    ColumnVector m_Signal_old; 
    const float& m_d;
    const ColumnVector& m_alpha;
    const ColumnVector& m_beta;
    const Matrix& m_bvals;
 public:
    //constructors::
    Fibre( const ColumnVector& alpha, const ColumnVector& beta, 
	   const Matrix& bvals,const float& d,const float& ardfudge=1):
     m_ardfudge(ardfudge), m_d(d), m_alpha(alpha), m_beta(beta), m_bvals(bvals){

	m_th=M_PI/2;
	m_th_old=m_th;
	m_ph=0;
	m_ph_old=m_ph;
	m_f=0.01;
	m_f_old=m_f;
	m_lam=10;
	m_lam_old=m_lam;
	m_lam_jump=true;
	m_th_prop=0.2;
	m_ph_prop=0.2;
	m_f_prop=0.2;
	m_lam_prop=1;

	m_th_prior=0;
	compute_th_prior();
	
	m_ph_prior=0;
	compute_ph_prior();
	
	m_f_prior=0;
      compute_f_prior();
      
      m_lam_prior=0;
      compute_lam_prior();

      m_Signal.ReSize(alpha.Nrows());
      m_Signal=0;
      m_Signal_old=m_Signal;

      compute_prior();
      compute_signal();

      m_th_acc=0; m_th_rej=0;
      m_ph_acc=0; m_ph_rej=0;
      m_f_acc=0; m_f_rej=0;
      m_lam_acc=0; m_lam_rej=0;

    }
    Fibre(const ColumnVector& alpha, 
	  const ColumnVector& beta, const Matrix& bvals, const float& d, const float& ardfudge,
	  const float& th, const float& ph, const float& f, 
	  const float& lam, const bool lam_jump=true) : 
      m_th(th), m_ph(ph), m_f(f), m_lam(lam), m_ardfudge(ardfudge),m_lam_jump(lam_jump),  m_d(d),
      m_alpha(alpha), m_beta(beta), m_bvals(bvals)
     {
      m_th_old=m_th;
      m_ph_old=m_ph;
      m_f_old=m_f;
      m_lam_old=m_lam;
      m_th_prop=0.2;
      m_ph_prop=0.2;
      m_f_prop=0.2;
      m_lam_prop=1;
      
      m_th_prior=0;
      compute_th_prior();

      m_ph_prior=0;
      compute_ph_prior();
      
      m_f_prior=0;
      compute_f_prior();
      
      m_lam_prior=0;
      compute_lam_prior();


      m_Signal.ReSize(alpha.Nrows());
      m_Signal=0;
      m_Signal_old=m_Signal;

      compute_prior();
      compute_signal();
	    
      m_th_acc=0; m_th_rej=0;
      m_ph_acc=0; m_ph_rej=0;
      m_f_acc=0; m_f_rej=0;
      m_lam_acc=0; m_lam_rej=0;
     
    }
    Fibre(const Fibre& rhs): 
      m_d(rhs.m_d), m_alpha(rhs.m_alpha), m_beta(rhs.m_beta), m_bvals(rhs.m_bvals){
      *this=rhs;
    }

      
    ~Fibre(){}
    
    inline float get_th() const{ return m_th;}
    inline void set_th(const float th){ m_th=th; }
    
    inline float get_ph() const{ return m_ph;}
    inline void set_ph(const float ph){ m_ph=ph; }
    
    inline float get_f() const{ return m_f;}
    inline void set_f(const float f){ m_f=f; }
    
    inline float get_lam() const{ return m_lam;}
    inline void set_lam(const float lam){ m_lam=lam; }
    
    inline void report() const {
    OUT(m_th);
    OUT(m_ph);
    OUT(m_f);
    OUT(m_lam);
    OUT(m_th_prop);
    OUT(m_ph_prop);
    OUT(m_f_prop);
    OUT(m_lam_prop);
    OUT(m_th_old);
    OUT(m_ph_old);
    OUT(m_f_old);
    OUT(m_lam_old);
    OUT(m_th_prior);
    OUT(m_ph_prior);
    OUT(m_f_prior);
    OUT(m_lam_prior);
    OUT(m_th_old_prior);
    OUT(m_ph_old_prior);
    OUT(m_f_old_prior);
    OUT(m_lam_old_prior);
    OUT(m_prior_en);
    OUT(m_old_prior_en);
    OUT(m_th_acc); 
    OUT(m_th_rej);
    OUT(m_ph_acc);
    OUT(m_ph_rej); 
    OUT(m_f_acc);
    OUT(m_f_rej);
    OUT(m_lam_acc); 
    OUT(m_lam_rej);
    OUT(m_lam_jump);
      

    }
    
    inline const ColumnVector& getSignal() const{  
      return m_Signal;                      
    }
    
    inline void restoreSignal() {
      m_Signal=m_Signal_old;
    }
    inline void setSignal(const ColumnVector& Signal){
      m_Signal=Signal;
    }
    
    inline void setSignal(const int i, const float val){
      m_Signal(i)=val;
    }

    inline float get_prior() const{ return m_prior_en;}
    

    inline void update_proposals(){
      m_th_prop*=sqrt(float(m_th_acc+1)/float(m_th_rej+1));
      m_th_prop=min(m_th_prop,maxfloat);
      m_ph_prop*=sqrt(float(m_ph_acc+1)/float(m_ph_rej+1));
      m_ph_prop=min(m_ph_prop,maxfloat);
      m_f_prop*=sqrt(float(m_f_acc+1)/float(m_f_rej+1));
      m_f_prop=min(m_f_prop,maxfloat);
      m_lam_prop*=sqrt(float(m_lam_acc+1)/float(m_lam_rej+1));
      m_lam_prop=min(m_lam_prop,maxfloat);
      m_th_acc=0; 
      m_th_rej=0;
      m_ph_acc=0; 
      m_ph_rej=0;
      m_f_acc=0; 
      m_f_rej=0;
      m_lam_acc=0; 
      m_lam_rej=0;
    }
    
    
    inline bool compute_th_prior(){
      m_th_old_prior=m_th_prior;
      if(m_th==0){m_th_prior=0;}
      else{
	m_th_prior=-log(fabs(sin(m_th)/2));
      }
      return false; //instant rejection flag
    }
    inline bool compute_ph_prior(){
      m_ph_old_prior=m_ph_prior;
      m_ph_prior=0;
      return false;
    }
    inline bool compute_f_prior(bool can_use_ard=true){
      //note(gamma(lam+1)/(gamma(1)*gamma(lam)) = lam
      // the following is a beta distribution with alpha=0
      m_f_old_prior=m_f_prior;
      if (m_f<=0 | m_f>=1 )
	return true;
      else{
	//m_f_prior=-(log(m_lam) + (m_lam-1)*log(1-m_f));
	if(!can_use_ard ){
	  m_f_prior=0;
	}
	else{
	  if(m_lam_jump)
	    // m_f_prior=log(1-m_f)+2*log(fabs(log(1-m_f))); //marginalised with uniform prior on lambda
	    m_f_prior=std::log(m_f);
	    else
	    m_f_prior=0;
	}
	m_f_prior=m_ardfudge*m_f_prior;
	return false;
      }
    }
    
    inline bool compute_lam_prior(){
      m_lam_old_prior=m_lam_prior;
      if(m_lam <0 | m_lam > 1e16 )
	return true;
      else{
	m_lam_prior=0;
	return false;
      }
    }
    
    inline void compute_prior(){
      m_old_prior_en=m_prior_en;
      m_prior_en=m_th_prior+m_ph_prior+m_f_prior+m_lam_prior;
    }

     void compute_signal(){
       m_Signal_old=m_Signal;
       for (int i = 1; i <= m_alpha.Nrows(); i++){
 	float angtmp=cos(m_ph-m_beta(i))*sin(m_alpha(i))*sin(m_th) + cos(m_alpha(i))*cos(m_th);
 	angtmp=angtmp*angtmp;
 	m_Signal(i)=exp(-m_d*m_bvals(1,i)*angtmp);
       }
     }


    inline bool propose_th(){
      m_th_old=m_th;
      m_th+=normrnd().AsScalar()*m_th_prop;
      bool rejflag=compute_th_prior();//inside this it stores the old prior
      compute_prior();
      compute_signal();
      return rejflag;
    };
    
    inline void accept_th(){
      m_th_acc++;      
    }
    
    inline void reject_th(){
      m_th=m_th_old;
      m_th_prior=m_th_old_prior;
      m_prior_en=m_old_prior_en;
      m_Signal=m_Signal_old;//Is there a better way of doing this??
      m_th_rej++;
    }
    
    inline bool propose_ph(){
      m_ph_old=m_ph;
      m_ph+=normrnd().AsScalar()*m_ph_prop;
       bool rejflag=compute_ph_prior();//inside this it stores the old prior
      compute_prior();
      compute_signal();
      return rejflag;
    };
    
    inline void accept_ph(){
      m_ph_acc++;
    }
    
    inline void reject_ph(){
      m_ph=m_ph_old;
      m_ph_prior=m_ph_old_prior;
      m_prior_en=m_old_prior_en;
      m_Signal=m_Signal_old;//Is there a better way of doing this??
      m_ph_rej++;
    }
    
    inline bool propose_f( bool can_use_ard=true){
      m_f_old=m_f;
      m_f+=normrnd().AsScalar()*m_f_prop;
      bool rejflag=compute_f_prior(can_use_ard);
      compute_prior();
      return rejflag;
    };
    
    inline void accept_f(){
      m_f_acc++;
    }
    
    inline void reject_f(){
      m_f=m_f_old;
      m_f_prior=m_f_old_prior;
      m_prior_en=m_old_prior_en;
      m_f_rej++;
    }
    
    inline bool propose_lam(){
      if(m_lam_jump){
	m_lam_old=m_lam;
	m_lam+=normrnd().AsScalar()*m_lam_prop;
	bool rejflag=compute_lam_prior();
	compute_f_prior();
	compute_prior();
	return rejflag;
      }
      else {return true;}
    };
    
    inline void accept_lam(){
      m_lam_acc++;
    }
    
    inline void reject_lam(){
      m_lam=m_lam_old;
      m_lam_prior=m_lam_old_prior;
      m_prior_en=m_old_prior_en;
      m_lam_rej++;
    }
    
    
    Fibre& operator=(const Fibre& rhs){
      m_th=rhs.m_th;
      m_ph=rhs.m_ph;
      m_f=rhs.m_f;
      m_lam=rhs.m_lam;
      m_th_prop=rhs.m_th_prop;
      m_ph_prop=rhs.m_ph_prop;
      m_f_prop=rhs.m_f_prop;
      m_lam_prop=rhs.m_lam_prop;
      m_th_old=rhs.m_th_old;
      m_ph_old=rhs.m_ph_old;
      m_f_old=rhs.m_f_old;
      m_lam_old=rhs.m_lam_old;
      m_th_prior=rhs.m_th_prior;
      m_ph_prior=rhs.m_ph_prior;
      m_f_prior=rhs.m_f_prior;
      m_lam_prior=rhs.m_lam_prior;
      m_th_old_prior=rhs.m_th_old_prior;
      m_ph_old_prior=rhs.m_ph_old_prior;
      m_f_old_prior=rhs.m_f_old_prior;
      m_lam_old_prior=rhs.m_lam_old_prior;
      m_prior_en=rhs.m_prior_en;
      m_old_prior_en=rhs.m_old_prior_en;
      m_th_acc=rhs.m_th_acc; 
      m_th_rej=rhs.m_th_rej;
      m_ph_acc=rhs.m_ph_acc;
      m_ph_rej=rhs.m_ph_rej; 
      m_f_acc=rhs.m_f_acc;
      m_f_rej=rhs.m_f_rej;
      m_lam_acc=rhs.m_lam_acc; 
      m_lam_rej=rhs.m_lam_rej;
      m_lam_jump=rhs.m_lam_jump;
      m_Signal=rhs.m_Signal; 
      m_Signal_old=rhs.m_Signal_old; 
      return *this;
    }

    friend  ostream& operator<<(ostream& ostr,const Fibre& p);

  };
//overload <<
  inline ostream& operator<<(ostream& ostr,const Fibre& p){
    ostr<<p.m_th<<" "<<p.m_ph<<" "<<p.m_f<<endl;
    return ostr;
  }


  class Multifibre{
    vector<Fibre> m_fibres;
    float m_d;
    float m_d_old;
    float m_d_prop;
    float m_d_prior; 
    float m_d_old_prior;
    float m_d_acc;
    float m_d_rej;
    float m_S0;
    float m_S0_old;
    float m_S0_prop;
    float m_S0_prior;
    float m_S0_old_prior;
    float m_S0_acc;
    float m_S0_rej;
    float m_prior_en;
    float m_old_prior_en;
    float m_likelihood_en;
    float m_old_likelihood_en;
    float m_energy;
    float m_old_energy;
    float m_ardfudge;

    const ColumnVector& m_data;
    const ColumnVector& m_alpha;
    const ColumnVector& m_beta;
    const Matrix& m_bvals;
    
  public:
    Multifibre(const ColumnVector& data,const ColumnVector& alpha, 
	       const ColumnVector& beta, const Matrix& b, int N ,float ardfudge=1):
      m_ardfudge(ardfudge),m_data(data), m_alpha(alpha), m_beta(beta), m_bvals(b){
      
      //      initialise(Amat,N);
    }
    
    ~Multifibre(){}
    
    const vector<Fibre>& fibres() const{
      return m_fibres;
    }
    
    void addfibre(const float th, const float ph, const float f,const float lam, const bool lam_jump=true){
      Fibre fib(m_alpha,m_beta,m_bvals,m_d,m_ardfudge,th,ph,f,lam,lam_jump);
      m_fibres.push_back(fib);
    }

    void addfibre(){
      Fibre fib(m_alpha,m_beta,m_bvals,m_d,m_ardfudge);
      m_fibres.push_back(fib);
    }
    
    void initialise_energies(){
      compute_d_prior();
      compute_S0_prior();
      m_prior_en=0;
      compute_prior();
      m_likelihood_en=0;
      compute_likelihood();
      compute_energy();
    }
    
    void initialise_props(){
      m_S0_prop=m_S0/10; //must have set inital values before this is called;
      m_d_prop=m_d/10;
    }
    
    inline float get_d() const{ return m_d;}
    inline void set_d(const float d){ m_d=d; }
    
    inline float get_energy() const { return m_likelihood_en+m_prior_en;}
    inline float get_likelihood_energy() const { return m_likelihood_en;}
    
    inline float get_S0() const{ return m_S0;}
    inline void set_S0(const float S0){ m_S0=S0; }
    
    inline void report() const{
      OUT(m_d);
      OUT(m_d_old);
      OUT(m_d_prop);
      OUT(m_d_prior); 
      OUT(m_d_old_prior);
      OUT(m_d_acc);
      OUT(m_d_rej);
      OUT(m_S0);
      OUT(m_S0_old);
      OUT(m_S0_prop);
      OUT(m_S0_prior);
      OUT(m_S0_old_prior);
      OUT(m_S0_acc);
      OUT(m_S0_rej);
      OUT(m_prior_en);
      OUT(m_old_prior_en);
      OUT(m_likelihood_en);
      OUT(m_old_likelihood_en);
      OUT(m_energy);
      OUT(m_old_energy);
      for (unsigned int i=0;i<m_fibres.size();i++){cout <<"fibre "<<i<<endl;m_fibres[i].report();}
    }


    inline bool compute_d_prior(){
      m_d_old_prior=m_d_prior;
      if(m_d<0)
	return true;
      else{
	m_d_prior=0;
	return false;
      }
    }
    
    inline bool compute_S0_prior(){
      m_S0_old_prior=m_S0_prior;
      if(m_S0<0) return true;
      else
	{    
	  m_S0_prior=0;
	  return false;
	}
    }

    bool reject_f_sum(){
      float fsum=0;
      for(unsigned int f=0;f<m_fibres.size();f++){
	fsum+=m_fibres[f].get_f();	
      }
      return fsum>1; //true if sum(f) > 1 and therefore, we should reject f
    }
    
    inline void compute_prior(){
      m_old_prior_en=m_prior_en;
      m_prior_en=m_d_prior+m_S0_prior;
      for(unsigned int f=0;f<m_fibres.size(); f++){
	m_prior_en=m_prior_en+m_fibres[f].get_prior();
      } 
    }
    
    void compute_likelihood(){
      m_old_likelihood_en=m_likelihood_en;
      ColumnVector pred(m_alpha.Nrows());
      pred=0;
      float fsum=0;
      for(unsigned int f=0;f<m_fibres.size();f++){
	pred=pred+m_fibres[f].get_f()*m_fibres[f].getSignal();
	fsum+=m_fibres[f].get_f();
	
      }
      for(int i=1;i<=pred.Nrows();i++){
	pred(i)=pred(i)+(1-fsum)*exp(-m_d*m_bvals(1,i));
      }
      pred=pred*m_S0;
      
      float sumsquares=(m_data-pred).SumSquare();
      m_likelihood_en=(m_data.Nrows()/2)*log(sumsquares/2);
    }
    
  
    inline  void  compute_energy(){
      m_old_energy=m_energy;
      m_energy=m_prior_en+m_likelihood_en;
    }
    
    bool test_energy(){
      float tmp=exp(m_old_energy-m_energy);
      return (tmp>unifrnd().AsScalar());
    }
    
    inline void restore_energy(){
      m_prior_en=m_old_prior_en;
      m_likelihood_en=m_old_likelihood_en;
      m_energy=m_old_energy;
    }
    
    inline void restore_energy_no_lik(){
      m_prior_en=m_old_prior_en;
      m_energy=m_old_energy;
    }

    inline bool propose_d(){
      m_d_old=m_d;
      m_d+=normrnd().AsScalar()*m_d_prop;
      bool rejflag=compute_d_prior();//inside this it stores the old prior
      for(unsigned int f=0;f<m_fibres.size();f++)
	m_fibres[f].compute_signal();
      return rejflag;
    };
    
    inline void accept_d(){
      m_d_acc++;      
    }
    
    inline void reject_d(){
      m_d=m_d_old;
      m_d_prior=m_d_old_prior;
      m_prior_en=m_old_prior_en;
      for(unsigned int f=0;f<m_fibres.size();f++)
	m_fibres[f].restoreSignal();
      m_d_rej++;
    }
    
    inline bool propose_S0(){
      m_S0_old=m_S0;
      m_S0+=normrnd().AsScalar()*m_S0_prop;
      bool rejflag=compute_S0_prior();//inside this it stores the old prior
      return rejflag;
    };
    
    inline void accept_S0(){
      m_S0_acc++;      
    }
    
    inline void reject_S0(){
      m_S0=m_S0_old;
      m_S0_prior=m_S0_old_prior;
      m_prior_en=m_old_prior_en;
      m_S0_rej++;
    }
    
    inline void update_proposals(){
      m_d_prop*=sqrt(float(m_d_acc+1)/float(m_d_rej+1));
      m_d_prop=min(m_d_prop,maxfloat);
      m_S0_prop*=sqrt(float(m_S0_acc+1)/float(m_S0_rej+1));
      m_S0_prop=min(m_S0_prop,maxfloat);
      m_d_acc=0; 
      m_d_rej=0;
      m_S0_acc=0; 
      m_S0_rej=0;
      for(unsigned int f=0; f<m_fibres.size();f++){
	m_fibres[f].update_proposals();
      }
      
    }

    
    
    

    void jump( bool can_use_ard=true ){
      
      if(!propose_d()){
	compute_prior();
	compute_likelihood();
	compute_energy();
	if(test_energy()){
	  accept_d();
	}
	else{
	  reject_d();
	  restore_energy();
	}
      }
      else{ 
	reject_d();
      }
      
      if(!propose_S0()){
	compute_prior();
	compute_likelihood();
	compute_energy();
	if(test_energy())
	  accept_S0();
	else{
	  reject_S0();
	  restore_energy();
	}
      }
      else{
	reject_S0();
      } 
      
      for(unsigned int f=0;f<m_fibres.size();f++){
	//cout<<"pre th"<<m_fibres[f].get_th()<<endl;	
	//	cout << m_fibres[f].getSignal()<<endl;
	//cout<<"prior"<<m_prior_en<<endl;
	//cout<<"lik"<<m_likelihood_en<<endl;
	//cout<<"energy"<<m_energy<<endl;
	//cout <<endl;
	if(!m_fibres[f].propose_th()){
	  
	  compute_prior();
	  compute_likelihood();
	  compute_energy();
	  //cout<<"pre th"<<m_fibres[f].get_th()<<endl;	
	  //	cout << m_fibres[f].getSignal()<<endl;
	  //cout<<"prior"<<m_prior_en<<endl;
	  //cout<<"lik"<<m_likelihood_en<<endl;
	  //cout<<"energy"<<m_energy<<endl;
	  	  
	  if(test_energy()){
	    m_fibres[f].accept_th();
	  //cout<<"accepted"<<endl;}
	  }
	  else{
	    m_fibres[f].reject_th();
	    restore_energy();
	    //cout<<"rejected"<<endl;
	  }
	    //cout <<endl<<endl;
	}
	else {
	  m_fibres[f].reject_th();
	}
	
	  if(!m_fibres[f].propose_ph()){
	    compute_prior();
	    compute_likelihood();
	    compute_energy();
	    if(test_energy())
	      m_fibres[f].accept_ph();
	    else{
	      m_fibres[f].reject_ph();
	      restore_energy();
	    }
	  }
	  else{
	    m_fibres[f].reject_ph();
	  }
	  
	  if(!m_fibres[f].propose_f( can_use_ard )){
	    if(!reject_f_sum()){
	      compute_prior();
	      compute_likelihood();
	      compute_energy();
	      if(test_energy()){
		m_fibres[f].accept_f();
	      }
	      else{
		m_fibres[f].reject_f();
		restore_energy();
	      }
	    }
	    else{//else for rejectin fsum>1
	      m_fibres[f].reject_f();
	    }
	  }
	  else{//else for rejecting rejflag returned from propose_f()
	    m_fibres[f].reject_f();
	  }
	  
	  
//	   if(!m_fibres[f].propose_lam()){
// 	    compute_prior();
// 	    compute_energy();
// 	    if(test_energy()){
// 	      m_fibres[f].accept_lam();
// 	    }
// 	    else{
// 	      m_fibres[f].reject_lam();
// 	      restore_energy_no_lik();
// 	    }
// 	  }
// 	  else{
// 	    m_fibres[f].reject_lam();
// 	  }
      
 
      }
      
      
      
    }
    
    
  };
  
  
}

#endif
