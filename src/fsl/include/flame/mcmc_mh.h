/*  mcmc_mh.h

    Mark Woolrich, Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(mcmc_mh_h)
#define mcmc_mh_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "gsoptions.h"
#include "newimage/newimageall.h"
#include "newmat.h"
#include "design.h"

using namespace NEWIMAGE;
using namespace NEWMAT;

namespace Gs {
    
  class Mcmc_Mh
    {
    public:

      // constructor
      Mcmc_Mh(const ColumnVector pcopedata, 
	      const ColumnVector pvarcopedata, 
	      const ColumnVector pdofvarcopedata, 
	      const Design& pdesign, 
	      ColumnVector pgamma_mean, 
	      ColumnVector pgamma_S, 
	      ColumnVector pbeta_b, 
	      ColumnVector pbeta_c, 
	      Matrix& pgamma_samples,
	      Matrix& pbeta_samples,
	      Matrix& pphi_samples,
	      ColumnVector& plikelihood_samples,
	      vector<ColumnVector>& pss_samples,
	      int pnsamples) :
	copedata(pcopedata),
	varcopedata(pvarcopedata),
	dofvarcopedata(pdofvarcopedata),
	design(pdesign),
	gamma_mean(pgamma_mean),
	gamma_S(),
	gamma_latest(design.getnevs()),
	gamma_samples(pgamma_samples),
	gamma_nrejected(design.getnevs()),
	gamma_naccepted(design.getnevs()),
	gamma_proposal_std(design.getnevs()),
	beta_b(pbeta_b),
	beta_c(pbeta_c),
	beta_latest(design.getngs()),
	beta_samples(pbeta_samples),
	beta_nrejected(design.getngs()),
	beta_naccepted(design.getngs()),
	beta_proposal_std(design.getngs()),
	beta_prior_energy_old(design.getngs()),
	phi_latest(design.getntpts()),
	phi_samples(pphi_samples),
	phi_nrejected(design.getntpts()),
	phi_naccepted(design.getntpts()),
	phi_proposal_std(design.getntpts()),
	phi_prior_energy_old(design.getntpts()),
	likelihood_energy_old(0.0),
	likelihood_samples(plikelihood_samples),
	//	ss_samples(pss_samples),
	ngs(design.getngs()),
	nevs(design.getnevs()),
	ntpts(design.getntpts()),
	nsamples(pnsamples),
	opts(GsOptions::getInstance()),
	delta(design.getntpts()),
	sampcount(0),
	subsampcount(0),
	sumovere(design.getntpts()),
	prec_ontwo(design.getntpts()),
	logprec_ontwo(design.getntpts()),
	uncertainty_in_varcopes(opts.dofvarcopefile.value() != string(""))
	{ 
	  reshape(gamma_S,pgamma_S,nevs,nevs);
	}

      // load data from file in from file and set up starting values
      void setup();

      // runs the chain
      void run();

      // jumps
      void jump();

      // sample chain
      void sample(int samp);
      
      // DIC
/*       void dic(float& DIC, float& pd); */

      // getters
      const int getnsamples() const { return nsamples; }   

      const ColumnVector& getgamma_naccepted() const { return gamma_naccepted; }
      const ColumnVector& getgamma_nrejected() const { return gamma_nrejected; } 
      const ColumnVector& getbeta_naccepted() const { return beta_naccepted; }
      const ColumnVector& getbeta_nrejected() const { return beta_nrejected; } 
      const ColumnVector& getphi_naccepted() const { return phi_naccepted; }
      const ColumnVector& getphi_nrejected() const { return phi_nrejected; }  
      // Destructor
      virtual ~Mcmc_Mh() {}

      ColumnVector c_samples; 
    private:
    
      void beta_jump();
      void phi_jump();
      void gamma_jump();
      void all_jump();

      float likelihood_energy(const int echanged, const float gamma_old, const bool betachanged);

      float likelihood_energy_phichanged(const int t);

      float beta_prior_energy(int g);
      float phi_prior_energy(int g);

/*       void sample_sumsquares(int samp); */
/*       float sumsquare_residuals(const Matrix& pdm, const ColumnVector& pdata, const ColumnVector& ppes); */

      Mcmc_Mh();
      const Mcmc_Mh& operator=(Mcmc_Mh& mcmc_mh);     
      Mcmc_Mh(Mcmc_Mh& mcmc_mh);

      const ColumnVector copedata;
      const ColumnVector varcopedata;
      const ColumnVector dofvarcopedata;
      const Design& design;

      // mean
      ColumnVector gamma_mean;

      // Covariance:
      Matrix gamma_S;

      ColumnVector gamma_latest;
      Matrix& gamma_samples;
      ColumnVector gamma_nrejected;
      ColumnVector gamma_naccepted;
      ColumnVector gamma_proposal_std;

      ColumnVector beta_b;
      ColumnVector beta_c;
      ColumnVector beta_latest;
      Matrix& beta_samples;
      ColumnVector beta_nrejected;
      ColumnVector beta_naccepted;
      ColumnVector beta_proposal_std;
      ColumnVector beta_prior_energy_old;

      ColumnVector phi_latest;
      Matrix& phi_samples;
      ColumnVector phi_nrejected;
      ColumnVector phi_naccepted;
      ColumnVector phi_proposal_std;
      ColumnVector phi_prior_energy_old;

      float likelihood_energy_old;

      ColumnVector& likelihood_samples;

/*       vector<ColumnVector>& ss_samples; */

      int ngs;
      int nevs;
      int ntpts;

      int nsamples;

      GsOptions& opts;

      float c_latest;

      ColumnVector delta;
      
      int sampcount;
      int subsampcount;

      ColumnVector sumovere;
      ColumnVector prec_ontwo;
      ColumnVector logprec_ontwo;

      bool uncertainty_in_varcopes;
    };
}   
#endif

