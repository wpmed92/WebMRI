/*  mcmc.h

    Mark Woolrich, Tim Behrens - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(mcmc_h)
#define mcmc_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "gsoptions.h"
#include "newimage/newimageall.h"
#include "newmat.h"
#include "design.h"
#include "gam.h"

using namespace NEWIMAGE;
using namespace NEWMAT;

namespace Gs {
    
  class Mcmc
    {
    public:

      // constructor
      Mcmc(const ColumnVector pcopedata, 
	   const ColumnVector pvarcopedata, 
	   const Design& pdesign, 
	   ColumnVector pgamma_mean, 
	   ColumnVector pgamma_S, 
	   ColumnVector pbeta_b, 
	   ColumnVector pbeta_c, 
	   ColumnVector pm_mean, 
	   ColumnVector pm_var, 
	   Matrix& pgamma_samples, 
	   int pnsamples, bool pswitched = false) :
	copedata(pcopedata),
	varcopedata(pvarcopedata),
	design(pdesign),
	gamma_mean(pgamma_mean),
	gamma_S(),
	gamma_samples(pgamma_samples),
	beta_b(pbeta_b),
	beta_c(pbeta_c),
	m_mean(pm_mean),
	m_var(pm_var),
	ngs(design.getngs()),
	nevs(design.getnevs()),
	ntpts(design.getntpts()),
	nsamples(pnsamples),
	opts(GsOptions::getInstance()),
	gam(Gam::getInstance()),	
	switched(pswitched),
	prior_dominating(false)
	{ 
	  reshape(gamma_S,pgamma_S,nevs,nevs);
	}

      // load data from file in from file and set up starting values
      void setup();

      // runs the chain
      void run();

      // jumps
      void jump(bool relax);

      // sample chain
      void sample(int samp);

      // getters
      const int getnsamples() const { return nsamples; }      
      const bool is_prior_dominating() const { return prior_dominating; }

      // Destructor
      virtual ~Mcmc() {}

      ColumnVector c_samples; 
    private:
    
      void beta_jump(bool relax);
      void m_jump(bool relax);
      void gamma_jump(bool relax);

      void beta_jump_switched(bool relax);
      void m_jump_switched(bool relax);
      void gamma_jump_switched(bool relax);
      
      Mcmc();
      const Mcmc& operator=(Mcmc& mcmc);     
      Mcmc(Mcmc& mcmc);

      const ColumnVector copedata;
      const ColumnVector varcopedata;
      const Design& design;

      ColumnVector gamma_mean;
      Matrix gamma_S;
      ColumnVector gamma_latest;
      Matrix& gamma_samples;

      ColumnVector beta_b;
      ColumnVector beta_c;
      ColumnVector beta_latest;

      ColumnVector m_mean;
      ColumnVector m_var;
      ColumnVector m_latest;

      int ngs;
      int nevs;
      int ntpts;

      int nsamples;

      GsOptions& opts;

      Gam& gam;

      bool switched;

      float c_latest;

      bool prior_dominating;
      
    };
}   
#endif

