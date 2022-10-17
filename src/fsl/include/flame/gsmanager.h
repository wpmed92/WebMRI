/*  gsmanager.h

    Mark Woolrich, Tim Behrens - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(gsmanager_h)
#define gsmanager_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "gsoptions.h"
#include "newimage/newimageall.h"
#include "design.h"
#include "miscmaths/volume.h"
#include "miscmaths/volumeseries.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

namespace Gs {
    
  // Give this class a file containing
  class Gsmanager
    {
    public:

      // constructor
      Gsmanager() : 
	opts(GsOptions::getInstance()),
/* 	ntapsi(16), */
/* 	ntpsi(23), */
/* 	psics(ntpsi), */
/* 	apsics(ntapsi), */
	nmaskvoxels(0)
/* 	doflut(50) */
	{
/* 	  doflut = 0; */
/* 	  for(int d=1;d <=doflut.Nrows();d++) */
/* 	    { */
/* 	      doflut(d)=0.0005*pow(d,4.0)-0.0448*pow(d,3.0)+1.3134*pow(d,2.0)-9.1104*d+26.1860; */
/* 	    } */

/* 	  psics << -.038057080835217922E0<< */
/* 	    .49141539302938713E0<< */
/* 	    -.056815747821244730E0<< */
/* 	    .008357821225914313E0<< */
/* 	    -.001333232857994342E0<< */
/* 	    .000220313287069308E0<< */
/* 	    -.000037040238178456E0<< */
/* 	    .000006283793654854E0<< */
/* 	    -.000001071263908506E0<< */
/* 	    .000000183128394654E0<< */
/* 	    -.000000031353509361E0<< */
/* 	    .000000005372808776E0<< */
/* 	    -.000000000921168141E0<< */
/* 	    .000000000157981265E0<< */
/* 	    -.000000000027098646E0<< */
/* 	    .000000000004648722E0<< */
/* 	    -.000000000000797527E0<< */
/* 	    .000000000000136827E0<< */
/* 	    -.000000000000023475E0<< */
/* 	    .000000000000004027E0<< */
/* 	    -.000000000000000691E0<< */
/* 	    .000000000000000118E0<< */
/* 	    -.000000000000000020E0; */
	  
/* 	  apsics <<-.0204749044678185E0<< */
/* 	    -.0101801271534859E0<< */
/* 	    .0000559718725387E0<< */
/* 	    -.0000012917176570E0<< */
/* 	    .0000000572858606E0<< */
/* 	    -.0000000038213539E0<< */
/* 	    .0000000003397434E0<< */
/* 	    -.0000000000374838E0<< */
/* 	    .0000000000048990E0<< */
/* 	    -.0000000000007344E0<< */
/* 	    .0000000000001233E0<< */
/* 	    -.0000000000000228E0<< */
/* 	    .0000000000000045E0<< */
/* 	    -.0000000000000009E0<< */
/* 	    .0000000000000002E0<< */
/* 	    -.0000000000000000E0; */
	}

      // load data from file in from file and set up starting values
      void setup();

      // initialise
      void initialise();

      void run();

      // saves results in logging directory
      void save();

      // does ols
      void ols(); 

      // Destructor
      virtual ~Gsmanager() {}
 
    private:

      const Gsmanager& operator=(Gsmanager& par);     
      Gsmanager(Gsmanager& des);
      
/*       void variational_bayes(int vox, ColumnVector& gammeanout, Matrix& gamcovout); */

      void multitfit(const Matrix& x, ColumnVector& m, SymmetricMatrix& covar, float& v, bool fixmean=false) const;

/*       float logtpdf(const float& v, const ColumnVector& x2, const float& phi) const; */
/*       float digamma(float x) const; */
/*       float csevl(const float x, const ColumnVector& cs, const int n) const; */
/*       float mgradpt(const float v, const ColumnVector& xsq, const int P) const; */
/*       float mdofls(const ColumnVector& xsq, const float phi, const int P) const; */

      float marg_posterior_energy(float x, const ColumnVector& y, const Matrix& z, const ColumnVector& S);

      float solveforbeta(const ColumnVector& y, const Matrix& z, const ColumnVector& S);

      bool pass_through_to_mcmc(float zlowerthresh, float zupperthresh, int vox);

      void ols_contrasts(const ColumnVector& gammean, const Matrix& gamS, int vox);

      void em_contrasts(const ColumnVector& gammean, const Matrix& gamS, int vox);

      void t_ols_contrast(const ColumnVector& gammean, const Matrix& gamS, const RowVector& tcontrast, double& cope, double& varcope, double& t, double& dof, double& z, bool lookupdof, int vox);
	
      void f_ols_contrast(const ColumnVector& gammean, const Matrix& gamS, const Matrix& fcontrast, double& f, double& dof1, double& dof2, double& z, bool lookupdof, int vox);
	
      void mcmc_contrasts(const Matrix& gamsamples, int vox);
	
      void t_mcmc_contrast(const Matrix& gamsamples, const RowVector& tcontrast, double& cope, double& varcope, double& t, double& dof, double& z, int vox);
	
      void f_mcmc_contrast(const Matrix& gamsamples, const Matrix& fcontrast, double& f, double& dof1, double& dof2, double& z, int vox);

      void marginal_em();
      void strap_on();

/*       int doflookuptable(int dofin) const {  */
/* 	if(dofin > doflut.Nrows()) */
/* 	  return dofin; */
/* 	else */
/* 	  return int(doflut(dofin));  */
/*       } */

      // inputs
      VolumeSeries copedata;
      VolumeSeries varcopedata;
      VolumeSeries dofvarcopedata;
      volume<float> mask;
      volume<float> zmask;

      // intermediates
      Design design;

      vector<Volume> beta_b;
      vector<Volume> beta_c;
      vector<Volume> beta_mean;

      VolumeSeries cov_pes;

      // outputs
/*       vector<volume<float> > pes; */
/*       vector<volume<float> > ts; */
/*       vector<volume<float> > tdofs; */
/*       vector<volume<float> > zts; */
/*       vector<volume<float> > tcopes; */
/*       vector<volume<float> > tvarcopes;       */
/*       vector<volume<float> > fs; */
/*       vector<volume<float> > fdof1s; */
/*       vector<volume<float> > fdof2s; */
/*       vector<volume<float> > zfs; */

      vector<Volume> pes;
      vector<Volume> ts;
      vector<Volume> tdofs;
      vector<Volume> zts;
      vector<Volume> zolsts;
      vector<Volume> zemupperts;
      vector<Volume> zemlowerts;
      vector<Volume> tcopes;
      vector<Volume> tvarcopes;      
      vector<Volume> fs;
      vector<Volume> fdof1s;
      vector<Volume> fdof2s;
      vector<Volume> zfs;
      vector<Volume> zolsfs;
      vector<Volume> zemupperfs;
      vector<Volume> zemlowerfs;

      // intermediates
      int ngs;
      int nevs;
      int ntpts;
      int xsize;
      int ysize;
      int zsize;

      GsOptions& opts;

//        const int ntapsi;
//        const int ntpsi;
//        ColumnVector psics;
//        ColumnVector apsics;

      int nmaskvoxels;

      volumeinfo volinfo;
      Volume sparsemask;

      bool dofpassedin;
/*       ColumnVector doflut; */
    };
}   
#endif







