/*  gsmanager.cc

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

#include "gsmanager.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
//#include "mcmc.h"
#include "mcmc_mh.h"
#include "miscmaths/t2z.h"
#include "miscmaths/f2z.h"

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace Gs {

float Gsmanager::marg_posterior_energy(float x, const ColumnVector& y, const Matrix& z, const ColumnVector& S)
{
  Tracer tr("Gsmanager::marg_posterior_energy");  

  // x is log(variance)

  float ex=exp(x);

  // ex is variance

  if(ex <= 0)
    {
      float ret =1e32;
      return ret;
    }

  DiagonalMatrix iU(ntpts);
  iU = 0;
  for(int t=1;t<=ntpts;t++)
    { 
      if(ex+S(t) == 0) return 1e32;
      iU(t) = 1.0/(S(t)+ex);
    }

  SymmetricMatrix ziUz;
  ziUz << z.t()*iU*z;
  ColumnVector gam=ziUz.i()*z.t()*iU*y;  

//   OUT(iU);
//   OUT(iU.LogDeterminant().LogValue());
//   OUT(log(priorsum));

  float logprior = -x;
  float logdfxbydx = x;

  float ret = -(0.5*iU.LogDeterminant().LogValue()-0.5*ziUz.LogDeterminant().LogValue()-0.5*(y.t()*iU*y-gam.t()*ziUz*gam).AsScalar()+logprior+logdfxbydx);

//   OUT(ret);
//   exit(1);
  return ret;

}

float Gsmanager::solveforbeta(const ColumnVector& y, const Matrix& z, const ColumnVector& S)
{ 
  Tracer tr("Gsmanager::solveforbeta");

  // Value of golden section (1 + sqrt(5))/2.0
  double phi = 1.6180339887499;
  double cphi = 1.0 - 1.0/phi;
  double TOL = 1.0e-6;	// Maximal fractional precision
  double TINY = 1.0e-10;         // Can't use fractional precision when minimum is at 0
  
  // Bracket the minimum
  double br_min=log(1e-10);
  double br_mid=0;
  double br_max=log(1e10);

  int dir=1;
  double pt=1e-10;

  // Use Brent's algorithm to find minimum
  // Initialise the points and function values
  double w = br_mid;   	// Where second from minimum is
  double v = br_mid;   	// Previous value of w
  double x = v;   	// Where current minimum is
  double e = 0.0; 	// Distance moved on step before last

  // x is log(variance)
  double fx = marg_posterior_energy(pt+x*dir, y, z, S);

  double fv = fx;
  double fw = fx;
  double d = 0.0;
  double tol1 = TOL*abs(x) + TINY;
   
  float prec = 0.0001;
  int niters=100;
  double u = 0.0;
  for(int n = 1;n<=niters;n++)
    {
      double xm = 0.5*(br_min+br_max);
      // Make sure that tolerance is big enough
      tol1 = TOL*abs(x) + TINY;

      // Decide termination on absolute precision required by options(2)
      if (abs(x - xm) <= prec && br_max-br_min < 4*prec)
	break;
      
      // Check if step before last was big enough to try a parabolic step.
      // Note that this will fail on first iteration, which must be a golden
      // section step.
      if (abs(e) > tol1)
	{
	  // Construct a trial parabolic fit through x, v and w
	  double r = (fx - fv) * (x - w);
	  double q = (fx - fw) * (x - v);

	  double p = (x - v)*q - (x - w)*r;

	  q = 2.0 * (q - r);

	  if (q > 0.0 && p != 0)
	    p = -p;
	  q = abs(q);
	  // Test if the parabolic fit is OK
	  
	  if (abs(p) >= abs(0.5*q*e) || p <= q*(br_min-x) || p >= q*(br_max-x))
	    {

	      // No it isn't, so take a golden section step
	      if (x >= xm)
		e = br_min-x;
	      else
		e = br_max-x;

	      d = cphi*e;
	    }
	  else
	    {
	      
	      // Yes it is, so take the parabolic step
	      e = d;
	      d = p/q;
	      u = x+d;
	      if (u-br_min < 2*tol1 || br_max-u < 2*tol1)
		d = sign(xm-x)*tol1;
	      
	    }
	}
      else
	{
	  
	  // Step before last not big enough, so take a golden section step
	  if (x >= xm)
	    {
	      e = br_min - x;
	    }
	  else
	    {
	      e = br_max - x;
	    }

	  d = cphi*e;
	}
      
      // Make sure that step is big enough
      if (abs(d) >= tol1)
	{
	  u = x+d;
	}
      else
	{
	  u = x + sign(d)*tol1;
	}


      // Evaluate function at u
      double fu = marg_posterior_energy(pt+u*dir, y, z, S);

      // Reorganise bracket
      if (fu <= fx)
	{
	  if (u >= x)
	    br_min = x;
	  else
	    br_max = x;
	  
	  v = w; w = x; x = u;
	  fv = fw; fw = fv; fx = fu;
	}
      else
	{
	  if (u < x)
	    br_min = u;   
	  else
	    br_max = u;
	  
	  if (fu <= fw || w == x)
	    {
	      v = w; w = u;
	      fv = fw; fw = fu;
	    }
	  else if (fu <= fv || v == x || v == w)
	    {
	      v = u;
	      fv = fu;
	    }
	}      
    }

//    OUT(exp(x));

  return Max(1e-10,exp(x));
}

  void Gsmanager::multitfit(const Matrix& x, ColumnVector& m, SymmetricMatrix& covar, float& v, bool fixmean/*=false*/) const
  {
    Tracer_Plus trace("Gsmanager::multitfit");

    int n = x.Ncols();
    int P = x.Nrows();

    v = 10;

    Matrix mx;
    if(!fixmean)
      {	
	Matrix mn;
	remmean(x,mx,mn,2);
	m = mn.t().AsColumn();
      }
    else
      {
	mx=x;
	for(int ctr = 1; ctr <= x.Nrows(); ctr++) {
	  for (int ctr2 =1; ctr2 <= x.Ncols(); ctr2++) {
	    mx(ctr,ctr2)-=m(ctr);
	  }
	}
      }


    covar = cov(mx.t());
    float tmp = pow(covar.Determinant(),(1.0/P));    

    covar=covar/tmp;

    // xsq(i) = x(i)'*inv(c)*x(i)
    ColumnVector xsq(n);
    SymmetricMatrix invcovar = covar.i();

    float cosi = 0.0;    
    for(int i = 1; i <= n; i++)
      {
	xsq(i) = (mx.Column(i).t()*invcovar*mx.Column(i)).AsScalar();
	cosi += xsq(i);
      }
    cosi /= n-1;

    float phi = tmp;

    for(int i = 1; i <= 50; i++)
      {
	float newphi = 0.0;
	for(int j = 1; j <= n; j++)
	  {
	    float tau = phi*(v+P)/(v*phi+xsq(j));
	    newphi += tau*xsq(j);
	  }
	phi = newphi/float(n*P);

// 	write_ascii_matrix(xsq,'xsq');
// 	OUT(P);
// 	OUT(phi);
	//	v = mdofls(xsq,phi,P);
	v =  2.0/(1.0-phi/cosi);
//	OUT(v);
      } 

    covar = covar*phi;
  }

  void Gsmanager::setup()
  {
    Tracer_Plus trace("Gsmanager::setup");
    
    if(GsOptions::getInstance().debuglevel.value()==2)
      {
	cout << "******************************************" << endl
	     << "SETUP" << endl << "******************************************"
	     << endl;
      }

    dofpassedin = opts.dofvarcopefile.value() != string("");

    // read data
//     read_volume4D(copedata, opts.copefile.value(),volinfo);
//     read_volume4D(varcopedata, opts.varcopefile.value());
    copedata.read(opts.copefile.value());
    varcopedata.read(opts.varcopefile.value());

    // mask:
    read_volume(mask, opts.maskfile.value(),volinfo);

    sparsemask.read(opts.maskfile.value());

    zmask = mask;
    zmask.binarise(1e-16);  // newimage call
    nmaskvoxels = int(zmask.sum());

    sparsemask.threshold(1e-16);

    if(sparsemask.getVolumeSize() != nmaskvoxels) throw Exception("masks not consistent");
    
    // threshold using mask:
    copedata.setPreThresholdPositions(sparsemask.getPreThresholdPositions());
    copedata.thresholdSeries();

    varcopedata.setPreThresholdPositions(sparsemask.getPreThresholdPositions());
    varcopedata.thresholdSeries();

//     if(opts.dofvarcopefile.value() != string(""))
//       read_volume4D(dofvarcopedata, opts.dofvarcopefile.value());
    if(dofpassedin)
      {
	dofvarcopedata.read(opts.dofvarcopefile.value());
	dofvarcopedata.setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	dofvarcopedata.thresholdSeries();
      }
    else
      {
	dofvarcopedata = varcopedata;
	dofvarcopedata = 0;
      }

    if(GsOptions::getInstance().debuglevel.value()==2)
      {
	cout << "******************************************" << endl
	     << "Data Read Complete" << endl << "******************************************"
	     << endl;
      }

    // setup design
    design.setup();    
    
    // check design and data compatability
    if(design.getntpts() != copedata.tsize())
      {
	cout << "design.getntpts()=" << design.getntpts() << endl;
	cout << "copedata.tsize()=" << copedata.tsize() << endl;
	throw Exception("Cope data and design have different numbers of time points");
      }
    if(design.getntpts() != varcopedata.tsize())
      {	
	cout << "dm.ntpts()=" << design.getntpts() << endl;
	cout << "copedata.tsize()=" << varcopedata.tsize() << endl;
	
	throw Exception("Varcope data and design have different numbers of time points");
      }
    
    ngs = design.getngs();
    nevs = design.getnevs();
    ntpts = design.getntpts();
    xsize = mask.xsize();
    ysize = mask.ysize();
    zsize = mask.zsize();
    
    cout << "nevs=" << nevs << endl;
    cout << "ntpts=" << ntpts << endl;
    cout << "ngs=" << ngs << endl;
    cout << "nvoxels=" << nmaskvoxels << endl;

    initialise();

  }

  void Gsmanager::initialise()
  {
    Tracer_Plus trace("Gsmanager::initialise");

    pes.resize(design.getnevs());
    ts.resize(design.getnumtcontrasts());
    tdofs.resize(design.getnumtcontrasts());
    zts.resize(design.getnumtcontrasts());
    tcopes.resize(design.getnumtcontrasts());
    tvarcopes.resize(design.getnumtcontrasts());
    fs.resize(design.getnumfcontrasts());
    fdof1s.resize(design.getnumfcontrasts());
    fdof2s.resize(design.getnumfcontrasts());
    zfs.resize(design.getnumfcontrasts());
       
    if(GsOptions::getInstance().ols.value() || GsOptions::getInstance().justols.value())
      {
	zolsts.resize(design.getnumtcontrasts());
	zolsfs.resize(design.getnumfcontrasts());
      }

    zemlowerts.resize(design.getnumtcontrasts());
    zemupperts.resize(design.getnumtcontrasts());
    zemlowerfs.resize(design.getnumfcontrasts());
    zemupperfs.resize(design.getnumfcontrasts());

    beta_b.resize(design.getngs());
    beta_c.resize(design.getngs());
    beta_mean.resize(design.getngs());

    cov_pes.ReSize(nevs*nevs,nmaskvoxels);   
    cov_pes = 0;

    for(int g = 0; g < ngs; g++)
      {
	beta_mean[g].ReSize(nmaskvoxels);		
	beta_mean[g] = 0;
	beta_b[g].ReSize(nmaskvoxels);		
	beta_b[g] = 0;
	beta_c[g].ReSize(nmaskvoxels);		
	beta_c[g] = 0;
      }

    for(int e = 0; e < nevs; e++)
      {
	pes[e].ReSize(nmaskvoxels);
	pes[e] = 0;
      }

    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	ts[t].ReSize(nmaskvoxels);
	ts[t] = 0;
	tdofs[t].ReSize(nmaskvoxels);

	tdofs[t] = 0;
	zts[t].ReSize(nmaskvoxels);
	zts[t] = 0;
	if(GsOptions::getInstance().ols.value() || GsOptions::getInstance().justols.value())
	  {
	    zolsts[t].ReSize(nmaskvoxels);
	    zolsts[t] = 0;
	  }
	zemupperts[t].ReSize(nmaskvoxels);
	zemupperts[t] = 0;
	zemlowerts[t].ReSize(nmaskvoxels);
	zemlowerts[t] = 0;
	tcopes[t].ReSize(nmaskvoxels);
	tcopes[t] = 0;
	tvarcopes[t].ReSize(nmaskvoxels);
	tvarcopes[t] = 0;
      }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	fs[f].ReSize(nmaskvoxels);
	fs[f] = 0;
	fdof1s[f].ReSize(nmaskvoxels);
	fdof1s[f] = 0;
	fdof2s[f].ReSize(nmaskvoxels);
	fdof2s[f] = 0;
	zfs[f].ReSize(nmaskvoxels);
	zfs[f] = 0;
	if(GsOptions::getInstance().ols.value() || GsOptions::getInstance().justols.value())
	  {
	    zolsfs[f].ReSize(nmaskvoxels);
	    zolsfs[f] = 0;
	  }
	zemupperfs[f].ReSize(nmaskvoxels);
	zemupperfs[f] = 0;
	zemlowerfs[f].ReSize(nmaskvoxels);
	zemlowerfs[f] = 0;
      }

    
  } 

  void Gsmanager::run()
  {
    Tracer_Plus trace("Gsmanager::run");
    
    if(GsOptions::getInstance().justols.value() || GsOptions::getInstance().ols.value())
      {
	ols();
      }
      
    if(!GsOptions::getInstance().justols.value())
      {
	// call strap on
	  
	if(GsOptions::getInstance().em.value())
	  marginal_em();
	else
	  strap_on();

	if(!opts.fixed.value())
	  {
	    // MCMC
	    vector<volume4D<float> > gamma_samples(nevs);
	    vector<volume<float> > gamma_naccepted(nevs);
	    vector<volume<float> > gamma_nrejected(nevs);
	    vector<volume4D<float> > beta_samples(ngs);
	    vector<volume<float> > beta_naccepted(ngs);
	    vector<volume<float> > beta_nrejected(ngs);
	    vector<volume4D<float> > phi_samples(ntpts);
	    vector<volume<float> > phi_naccepted(ntpts);
	    vector<volume<float> > phi_nrejected(ntpts);
	    vector<volume4D<float> > ss_samples(design.getnumfcontrasts()+1);
    
	    volume<float> DIC(xsize,ysize,zsize);
	    volume<float> pd(xsize,ysize,zsize);
	    DIC = 0;
	    pd = 0;

	    int nsamples = (opts.njumps.value()-opts.burnin.value())/opts.sampleevery.value();
	    cout << "njumps = " << opts.njumps.value() << endl;
	    cout << "burnin = " << opts.burnin.value() << endl;
	    cout << "sampleevery = " << opts.sampleevery.value() << endl;
	    cout << "nsamples = " << nsamples << endl << endl;
        

	    if(GsOptions::getInstance().verbose.value())
	      {
		for(int e = 0; e < nevs; e++)
		  {	
		    gamma_samples[e].reinitialize(xsize,ysize,zsize,nsamples);
		    gamma_samples[e] = 0;
		    gamma_naccepted[e].reinitialize(xsize,ysize,zsize);
		    gamma_naccepted[e] = 0;
		    gamma_nrejected[e].reinitialize(xsize,ysize,zsize);
		    gamma_nrejected[e] = 0;
		  }

		if(dofpassedin)
		  for(int t = 0; t < ntpts; t++)
		    {
		      phi_samples[t].reinitialize(xsize,ysize,zsize,nsamples);
		      phi_samples[t] = 0;
		      phi_naccepted[t].reinitialize(xsize,ysize,zsize);
		      phi_naccepted[t] = 0;
		      phi_nrejected[t].reinitialize(xsize,ysize,zsize);
		      phi_nrejected[t] = 0;
		    }
	
		for(int g = 0; g < ngs; g++)
		  {
		    beta_samples[g].reinitialize(xsize,ysize,zsize,nsamples);
		    beta_samples[g] = 0;
		    beta_naccepted[g].reinitialize(xsize,ysize,zsize);
		    beta_naccepted[g] = 0;
		    beta_nrejected[g].reinitialize(xsize,ysize,zsize);
		    beta_nrejected[g] = 0;
		  }

		for(int f = 0; f < design.getnumfcontrasts()+1; f++)
		  {
		    ss_samples[f].reinitialize(xsize,ysize,zsize,nsamples);
		    ss_samples[f] = 0;
		  }

	      }

	    Matrix gamsamples(nevs, nsamples);
	    Matrix betasamples(ngs, nsamples);
	    Matrix phisamples(ntpts, nsamples);
	    ColumnVector likelihood_samples(nsamples);
	    vector<ColumnVector> sssamples(design.getnumfcontrasts()+1);

	    for(int f = 0; f < design.getnumfcontrasts()+1; f++)
	      {
		sssamples[f].ReSize(nsamples);
		sssamples[f] = 0;
	      }

	    int vox = 0;
	    int vox2 = 0;
	    int voxout = 0;

	    cout << "Metropolis Hasting Sampling" << endl;
	    cout << "Number of voxels=" << nmaskvoxels << endl;
	    cout << "Percentage done:" << endl;

	    //     VolumeSeries test = copedata;
	    //     test = 0;

	    for(int x = 0; x < xsize; x++)
	      for(int y = 0; y < ysize; y++)
		for(int z = 0; z < zsize; z++)
		  {
		    if(mask(x,y,z))
		      {
			vox++;
		
			//		test.setSeries(copedata.getSeries(vox),vox);

			if(zmask(x,y,z))
			  {
			    vox2++;
			    gamsamples = 0;
			    betasamples = 0;
			    phisamples = 0;

			    if(vox2 > voxout*nmaskvoxels/100.0)
			      {
				//cout<<(voxout+1)<<'%' <<'\r';		
				cout << " " << (voxout+1);
				cout.flush();
				voxout++;
			      }

			    if(GsOptions::getInstance().debuglevel.value()==2)
			      {
				cout << "--------------" << endl;
				cout << "x=" << x << "y=" << y << "z=" << z << endl;
			      }

			    ColumnVector gammamean(nevs);
			    for(int e = 0; e < nevs; e++)
			      {
				gammamean(e+1) = pes[e](vox);
			      }
			    ColumnVector betab(ngs);
			    ColumnVector betac(ngs);
			    for(int g = 0; g < ngs; g++)
			      {
				betab(g+1) = beta_b[g](vox);
				betac(g+1) = beta_c[g](vox);
			      }

			    srand(GsOptions::getInstance().seed.value());

			    Mcmc_Mh mcmc_mh(copedata.getSeries(vox), varcopedata.getSeries(vox),dofvarcopedata.getSeries(vox), design, gammamean, cov_pes.getSeries(vox), betab, betac, gamsamples, betasamples, phisamples, likelihood_samples, sssamples, nsamples);
			    mcmc_mh.setup();
			    mcmc_mh.run();
			    
			    if(!GsOptions::getInstance().fixmeanfortfit.value())
			      {
				for(int e = 0; e < nevs; e++)
				  {
				    pes[e](vox) = mean(gamsamples.Row(e+1).t()).AsScalar();
				  }
			      }
			    for(int g = 0; g < ngs; g++)
			      {
				beta_mean[g](vox) = mean(betasamples.Row(g+1).t()).AsScalar();
			      }

			    mcmc_contrasts(gamsamples,vox);

			    if((std::abs(zts[0](vox)-zemlowerts[0](vox))>3))
			    {
			      cout << endl << "WARNING: FLAME stage 2 has given an abnormally large difference to stage 1" << endl;
			      cout << "x=" << x << ",y=" << y << ",z=" << z << endl;
			      OUT(vox); 
			      OUT(zts[0](vox));
			      OUT(zemlowerts[0](vox));
			      OUT(beta_mean[0](vox));
			      OUT(beta_b[0](vox)/beta_c[0](vox));

			      for(int e = 0; e < nevs; e++)
				{
				  write_ascii_matrix(LogSingleton::getInstance().appendDir("gamma_"+num2str(vox+1)+"_"+num2str(e+1)),gamsamples.Row(e+1).t());
				}
			      for(int g = 0; g < ngs; g++)
				{
				  write_ascii_matrix(LogSingleton::getInstance().appendDir("beta_"+num2str(vox+1)+"_"+num2str(g+1)),betasamples.Row(g+1).t());
				}
			    }

			    //			mcmc_mh.dic(DIC(x,y,z), pd(x,y,z));

			    if(GsOptions::getInstance().verbose.value())
			      {
				for(int e = 0; e < nevs; e++)
				  {
				    gamma_samples[e].setvoxelts(gamsamples.Row(e+1).t(),x,y,z);
				    gamma_naccepted[e](x,y,z) = mcmc_mh.getgamma_naccepted()(e+1);
				    gamma_nrejected[e](x,y,z) = mcmc_mh.getgamma_nrejected()(e+1);
				  }	

				for(int f = 0; f < design.getnumfcontrasts()+1; f++)
				  {
				    ss_samples[f].setvoxelts(sssamples[f],x,y,z);
				  }
		    
				for(int g = 0; g < ngs; g++)
				  {
				    beta_samples[g].setvoxelts(betasamples.Row(g+1).t(),x,y,z);
				    beta_naccepted[g](x,y,z) = mcmc_mh.getbeta_naccepted()(g+1);
				    beta_nrejected[g](x,y,z) = mcmc_mh.getbeta_nrejected()(g+1);
				  }

				if(dofpassedin)
				  for(int t = 0; t < ntpts; t++)
				    {
				      phi_samples[t].setvoxelts(phisamples.Row(t+1).t(),x,y,z);
				      phi_naccepted[t](x,y,z) = mcmc_mh.getphi_naccepted()(t+1);
				      phi_nrejected[t](x,y,z) = mcmc_mh.getphi_nrejected()(t+1);
				    }

			      }
			  }	      	    
		      }
		  }
	  

	    //     VolumeInfo inf = sparsemask.getInfo();
	    //     inf.v = test.tsize();
	    //     test.setInfo(inf);
	    //     test.setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	    //     test.unthresholdSeries();
	    //     test.writeAsFloat("test");
	
	    if(GsOptions::getInstance().verbose.value())
	      {
		for(int g = 0; g < ngs; g++)
		  {
		    save_volume4D(beta_samples[g], LogSingleton::getInstance().appendDir("beta"+num2str(g+1)+"_samples"));	
		    save_volume(beta_naccepted[g], LogSingleton::getInstance().appendDir("beta"+num2str(g+1)+"_naccepted"));
		    save_volume(beta_nrejected[g], LogSingleton::getInstance().appendDir("beta"+num2str(g+1)+"_nrejected"));	
		  }

		if(dofpassedin)
		  for(int t = 0; t < ntpts; t++)
		    {
		      save_volume4D(phi_samples[t], LogSingleton::getInstance().appendDir("phi"+num2str(t+1)+"_samples"));	
		      save_volume(phi_naccepted[t], LogSingleton::getInstance().appendDir("phi"+num2str(t+1)+"_naccepted"));
		      save_volume(phi_nrejected[t], LogSingleton::getInstance().appendDir("phi"+num2str(t+1)+"_nrejected"));	
		    }
	
		for(int e = 0; e < nevs; e++)
		  {
		    save_volume4D(gamma_samples[e], LogSingleton::getInstance().appendDir("gamma"+num2str(e+1)+"_samples"));
		    save_volume(gamma_naccepted[e], LogSingleton::getInstance().appendDir("gamma"+num2str(e+1)+"_naccepted"));
		    save_volume(gamma_nrejected[e], LogSingleton::getInstance().appendDir("gamma"+num2str(e+1)+"_nrejected"));
	    
		  }
	
		// 	for(int f = 0; f < design.getnumfcontrasts()+1; f++)
		// 	  {
		// 	    save_volume4D(ss_samples[f], LogSingleton::getInstance().appendDir("ss"+num2str(f+1)+"_samples"));
		// 	  }

		save_volume(DIC, LogSingleton::getInstance().appendDir("dic"));
		save_volume(pd, LogSingleton::getInstance().appendDir("pd"));
	      }
	  }
	cout << endl;
      }
  }

void Gsmanager::save()
{
    Tracer_Plus trace("Gsmanager::save");

    varcopedata.CleanUp();
    dofvarcopedata.CleanUp();

    // need to save residuals:
    VolumeSeries res = copedata;
    res = 0;
    
    const Matrix& dm = design.getdm();
    int vox = 0;
    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(mask(x,y,z))
	      {
		vox++;
		ColumnVector petmp(nevs);
		for(int e = 0; e < nevs; e++)
		  {
		    petmp(e+1) = pes[e](vox);
		  }

		res.setSeries(copedata.getSeries(vox)-dm*petmp,vox);
	      }
	  }

    VolumeInfo inf = sparsemask.getInfo();
    inf.v = res.tsize();
    res.setInfo(inf);
    res.setPreThresholdPositions(sparsemask.getPreThresholdPositions());
    res.unthresholdSeries();
    res.writeAsFloat(LogSingleton::getInstance().appendDir("res4d"));
    res.CleanUp();
    
    copedata.CleanUp();
    save_volume(mask, LogSingleton::getInstance().appendDir("mask"),volinfo);

    VolumeInfo info = sparsemask.getInfo();
    
    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	info.intent_code = NIFTI_INTENT_TTEST;
	ts[t].setInfo(info);
	ts[t].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	ts[t].unthreshold();
	ts[t].writeAsFloat(LogSingleton::getInstance().appendDir("tstat"+num2str(t+1)));
	ts[t].CleanUp();

	info.intent_code = NIFTI_INTENT_NONE;
	tdofs[t].setInfo(info);
	tdofs[t].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	tdofs[t].unthreshold();
	tdofs[t].writeAsFloat(LogSingleton::getInstance().appendDir("tdof_t"+num2str(t+1)));
	tdofs[t].CleanUp();

	info.intent_code = NIFTI_INTENT_ZSCORE;
	zts[t].setInfo(info);
	zts[t].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	zts[t].unthreshold();
	zts[t].writeAsFloat(LogSingleton::getInstance().appendDir("zstat"+num2str(t+1)));
	zts[t].CleanUp();

	if(GsOptions::getInstance().ols.value())
	  {
	    info.intent_code = NIFTI_INTENT_ZSCORE;
	    zolsts[t].setInfo(info);
	    zolsts[t].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	    zolsts[t].unthreshold();
	    zolsts[t].writeAsFloat(LogSingleton::getInstance().appendDir("zolststat"+num2str(t+1)));
	    zolsts[t].CleanUp();
	  }

	info.intent_code = NIFTI_INTENT_ZSCORE;
	zemupperts[t].setInfo(info);
	zemupperts[t].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	zemupperts[t].unthreshold();
	zemupperts[t].writeAsFloat(LogSingleton::getInstance().appendDir("zemuppertstat"+num2str(t+1)));
	zemupperts[t].CleanUp();

	info.intent_code = NIFTI_INTENT_ZSCORE;
	zemlowerts[t].setInfo(info);
	zemlowerts[t].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	zemlowerts[t].unthreshold();
	zemlowerts[t].writeAsFloat(LogSingleton::getInstance().appendDir("zemlowertstat"+num2str(t+1)));
	zemlowerts[t].CleanUp();

	info.intent_code = NIFTI_INTENT_ESTIMATE;
	tcopes[t].setInfo(info);
	tcopes[t].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	tcopes[t].unthreshold();
	tcopes[t].writeAsFloat(LogSingleton::getInstance().appendDir("cope"+num2str(t+1)));
	tcopes[t].CleanUp();

	info.intent_code = NIFTI_INTENT_ESTIMATE;
	tvarcopes[t].setInfo(info);
	tvarcopes[t].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	tvarcopes[t].unthreshold();
	tvarcopes[t].writeAsFloat(LogSingleton::getInstance().appendDir("varcope"+num2str(t+1)));
	tvarcopes[t].CleanUp();

      }
// 	save_volume(ts[t], LogSingleton::getInstance().appendDir("tstat"+num2str(t+1)),volinfo);
// 	save_volume(tdofs[t], LogSingleton::getInstance().appendDir("tdof_t"+num2str(t+1)),volinfo);
// 	save_volume(zts[t], LogSingleton::getInstance().appendDir("zstat"+num2str(t+1)),volinfo);
// 	save_volume(tcopes[t], LogSingleton::getInstance().appendDir("cope"+num2str(t+1)),volinfo);
// 	save_volume(tvarcopes[t], LogSingleton::getInstance().appendDir("varcope"+num2str(t+1)),volinfo);
//       }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	info.intent_code = NIFTI_INTENT_FTEST;
	fs[f].setInfo(info);
	fs[f].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	fs[f].unthreshold();
	fs[f].writeAsFloat(LogSingleton::getInstance().appendDir("fstat"+num2str(f+1)));
	fs[f].CleanUp();
	
	if(GsOptions::getInstance().ols.value())
	  {
	    info.intent_code = NIFTI_INTENT_ZSCORE;
	    zolsfs[f].setInfo(info);
	    zolsfs[f].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	    zolsfs[f].unthreshold();
	    zolsfs[f].writeAsFloat(LogSingleton::getInstance().appendDir("zolsfstat"+num2str(f+1)));
	    zolsfs[f].CleanUp();
	  }

	info.intent_code = NIFTI_INTENT_ZSCORE;
	zemupperfs[f].setInfo(info);
	zemupperfs[f].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	zemupperfs[f].unthreshold();
	zemupperfs[f].writeAsFloat(LogSingleton::getInstance().appendDir("zemupperfstat"+num2str(f+1)));
	zemupperfs[f].CleanUp();

	info.intent_code = NIFTI_INTENT_ZSCORE;
	zemlowerfs[f].setInfo(info);
	zemlowerfs[f].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	zemlowerfs[f].unthreshold();
	zemlowerfs[f].writeAsFloat(LogSingleton::getInstance().appendDir("zemlowerfstat"+num2str(f+1)));
	zemlowerfs[f].CleanUp();

	info.intent_code = NIFTI_INTENT_NONE;
	fdof1s[f].setInfo(info);
	fdof1s[f].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	fdof1s[f].unthreshold();
	fdof1s[f].writeAsFloat(LogSingleton::getInstance().appendDir("fdof1_f"+num2str(f+1)));
	fdof1s[f].CleanUp();

	info.intent_code = NIFTI_INTENT_NONE;
	fdof2s[f].setInfo(info);
	fdof2s[f].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	fdof2s[f].unthreshold();
	fdof2s[f].writeAsFloat(LogSingleton::getInstance().appendDir("fdof2_f"+num2str(f+1)));
	fdof2s[f].CleanUp();

	info.intent_code = NIFTI_INTENT_ZSCORE;
	zfs[f].setInfo(info);
	zfs[f].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	zfs[f].unthreshold();
	zfs[f].writeAsFloat(LogSingleton::getInstance().appendDir("zfstat"+num2str(f+1)));
	zfs[f].CleanUp();
      }

//     for(int f = 0; f < design.getnumfcontrasts(); f++)
//       {
// 	save_volume(fs[f], LogSingleton::getInstance().appendDir("fstat"+num2str(f+1)),volinfo);
// 	save_volume(fdof1s[f], LogSingleton::getInstance().appendDir("fdof1_f"+num2str(f+1)),volinfo);
// 	save_volume(fdof2s[f], LogSingleton::getInstance().appendDir("fdof2_f"+num2str(f+1)),volinfo);
// 	save_volume(zfs[f], LogSingleton::getInstance().appendDir("zfstat"+num2str(f+1)),volinfo);
//       }

    for(int e = 0; e < nevs; e++)
      {
	info.intent_code = NIFTI_INTENT_ESTIMATE;
	pes[e].setInfo(info);
	pes[e].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	pes[e].unthreshold();
	pes[e].writeAsFloat(LogSingleton::getInstance().appendDir("pe"+num2str(e+1)));
	pes[e].CleanUp();
      }

//     for(int e = 0; e < nevs; e++)
//       {
// 	save_volume(pes[e], LogSingleton::getInstance().appendDir("pe")+num2str(e+1),volinfo);
//       }


    for(int g = 0; g < ngs; g++)
      {
	info.intent_code = NIFTI_INTENT_ESTIMATE;
	beta_mean[g].setInfo(info);
	beta_mean[g].setPreThresholdPositions(sparsemask.getPreThresholdPositions());
	beta_mean[g].unthreshold();
	beta_mean[g].writeAsFloat(LogSingleton::getInstance().appendDir("mean_random_effects_var"+num2str(g+1)));
	beta_mean[g].CleanUp();
      }

//     save_volume4D(beta_mean, LogSingleton::getInstance().appendDir("mean_random_effects_var"),volinfo);

    
  }

  void Gsmanager::ols()
  {
    Tracer_Plus trace("Gsmanager::ols");

    const Matrix& dm = design.getdm();

    Matrix pinvdm = (dm.t()*dm).i();

    Matrix pinvdmdm = pinvdm*dm.t();

    if(nevs >= ntpts)
      {
	throw Exception("Singular design. Number of EVs > number of time points. ");
      }
 
    int vox = 0;
    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(mask(x,y,z))
	      {
		vox++;
		ColumnVector Y = copedata.getSeries(vox);

		ColumnVector petmp = pinvdmdm*Y;

		for(int e = 0; e < nevs; e++)
		  {
		    pes[e](vox) = petmp(e+1);
		  }

		ColumnVector r = Y-dm*petmp;
		float r2 = (r.t()*r).AsScalar();
		
		Matrix covariance = (pinvdm*r2/(ntpts-nevs));
		ColumnVector covariance2;
		reshape(covariance2, covariance, nevs*nevs, 1);

		cov_pes.setSeries(covariance2,vox);
		
		ols_contrasts(petmp,covariance,vox);

	      }
	  }
    
  }

  void Gsmanager::strap_on()
  {
    Tracer_Plus trace("Gsmanager::strap_on");

    const Matrix& dm = design.getdm();    

    if(!opts.fixed.value() && nevs >= ntpts)
      {
	throw Exception("Singular matrix.");
      }   
    
//     OUT(dm);
    Matrix pinvdm = pinv(dm);

    // setup design matrices for groups
    vector<Matrix> zg(ngs);
    vector<Matrix> zzg(ngs);
    vector<float> traceRg(ngs);
    for(int g = 1; g <= ngs; g++)
      {
	zg[g-1].ReSize(design.getntptsingroup(g),nevs);
	zg[g-1]=0;
	int t2=1;
	for(int t = 1; t <= ntpts; t++)
	  {
	    if(design.getgroup(t)==g)
	      {		
		zg[g-1].Row(t2++) = dm.Row(t);
	      }
	  }	

	traceRg[g-1] = (Identity(design.getntptsingroup(g))-zg[g-1]*pinv(zg[g-1])).Trace();
//  	OUT(zg[g-1]);
//  	OUT(traceRg[g-1]);

	zzg[g-1] = zg[g-1].t()*zg[g-1];

      }
    
    int vox = 0;
    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(mask(x,y,z))
	      {
		vox++;

// 		cout << x << "," << y << "," << z << endl;
		// setup data
		ColumnVector Y = copedata.getSeries(vox);
		ColumnVector S = varcopedata.getSeries(vox);
//  		OUT(Y.t());
//  		OUT(S.t());
				

		// setup data for groups
		vector<ColumnVector> Yg(ngs);
		vector<ColumnVector> Sg(ngs);
		for(int g = 1; g <= ngs; g++)
		  {
		    Yg[g-1].ReSize(design.getntptsingroup(g));
		    Yg[g-1] = 0;
		    Sg[g-1].ReSize(design.getntptsingroup(g));
		    Sg[g-1] = 0; 
		    int t2=1;
		    for(int t = 1; t <= ntpts; t++)
		      {
			if(design.getgroup(t)==g)
			  {
			    Yg[g-1](t2) = Y(t);
			    Sg[g-1](t2) = S(t);
			    t2++;
			  }
		      }
//  		    OUT(Yg[g-1].t());
//  		    OUT(Sg[g-1].t());
		  }
		
		// calc betas
		ColumnVector beta(ngs);
		beta = 0.0;

		if(!opts.fixed.value())
		  {
		    for(int g = 1; g <= ngs; g++)
		      {
			// 		    if(x==16 && y==41 && z==48)
			// 		      {
			// 			write_ascii_matrix(Yg[g-1],"/usr/people/woolrich/tmp/Y");
			// 			write_ascii_matrix(zg[g-1],"/usr/people/woolrich/tmp/z");
			// 			write_ascii_matrix(Sg[g-1],"/usr/people/woolrich/tmp/Sg");
			// 		      }
			beta(g) = solveforbeta(Yg[g-1],zg[g-1],Sg[g-1]);
		      }

		    //		    OUT(log(beta(1)));
		  }

		// calc gam
		DiagonalMatrix iU(ntpts);
		iU = 0;
		
		for(int t=1;t<=ntpts;t++)
		  {		    
		    iU(t) = 1.0/(S(t)+beta(design.getgroup(t)));
		  }

		ColumnVector gam=(dm.t()*iU*dm).i()*dm.t()*iU*Y;  
// for MCMC:	       

 		//OUT(gam.t());
 		//OUT(beta.t());
		for(int g = 0; g < ngs; g++)
		  {
		    float betamean = beta(g+1);
		    float betavar = 1;
		    beta_b[g](vox) = Sqr(betamean)/betavar;
		    beta_c[g](vox) = betamean/betavar;
		  }

		// store results for gam:
		for(int e = 0; e < nevs; e++)
		  {
		    pes[e](vox) = gam(e+1);			
		  }

		// store results for beta:
		for(int g = 0; g < ngs; g++)
		  {
		    beta_mean[g](vox) = beta(g+1);
		  }
		
		
		// store results for cov		

		Matrix gamcovariance = (dm.t()*iU*dm).i();

		ColumnVector gamcovariance2;
		reshape(gamcovariance2, gamcovariance, nevs*nevs, 1);		
		cov_pes.setSeries(gamcovariance2,vox);	        

		em_contrasts(gam,gamcovariance,vox);  

//  		OUT(beta(1));
// 		OUT(S.t());
// 		OUT(gam(1)/std::sqrt(gamcovariance(1,1)));
//  		OUT(gam);
//  		OUT(gamcovariance);
// 		OUT(beta_b[0](1)/beta_c[0](1));
// 		OUT(beta_b[0](1));
// 		OUT(beta_c[0](1))

		if(GsOptions::getInstance().zlowerthreshold.value() > 0)
		  {
		    // check variance ratio
		    ColumnVector vr(ngs);
		    for(int g = 0; g < ngs; g++)
		      {
			vr(g+1) = beta(g+1)/mean(Sg[g]).AsScalar();
		      }		    
		    bool likeols = min(vr).AsScalar()>10; 
		    
		    if(likeols)
		      {
			// don't do MCMC if variance ratio large
			// assume OLS DOF
			zmask(x,y,z) = false;			
		      }
		    else
		      {	

			zmask(x,y,z) = pass_through_to_mcmc(GsOptions::getInstance().zlowerthreshold.value(),GsOptions::getInstance().zupperthreshold.value(),vox);
		      }
		  }

	      }
	  }
    nmaskvoxels = int(zmask.sum()); 
    OUT(nmaskvoxels);
  }

  void Gsmanager::marginal_em()
  {
    Tracer_Plus trace("Gsmanager::marginal_em");

    const Matrix& dm = design.getdm();    

    if(nevs >= ntpts)
      {
	throw Exception("Singular matrix. ");
      }   
    
//     OUT(dm);
    Matrix pinvdm = pinv(dm);

    // setup design matrices for groups
    vector<Matrix> zg(ngs);
    vector<Matrix> zzg(ngs);
    vector<float> traceRg(ngs);
    for(int g = 1; g <= ngs; g++)
      {
	zg[g-1].ReSize(design.getntptsingroup(g),nevs);
	zg[g-1]=0;
	int t2=1;
	for(int t = 1; t <= ntpts; t++)
	  {
	    if(design.getgroup(t)==g)
	      {		
		zg[g-1].Row(t2++) = dm.Row(t);
	      }
	  }
	

	traceRg[g-1] = (Identity(design.getntptsingroup(g))-zg[g-1]*pinv(zg[g-1])).Trace();
//  	OUT(zg[g-1]);
//  	OUT(traceRg[g-1]);

	zzg[g-1] = zg[g-1].t()*zg[g-1];

      }
    
    int vox = 0;
    for(int x = 0; x < xsize; x++)
      for(int y = 0; y < ysize; y++)
	for(int z = 0; z < zsize; z++)
	  {
	    if(mask(x,y,z))
	      {
		vox++;

		// setup data
		ColumnVector Y = copedata.getSeries(vox);
		ColumnVector S = varcopedata.getSeries(vox);
//  		OUT(Y.t());
//  		OUT(S.t());
		
		// setup data for groups
		vector<ColumnVector> Yg(ngs);
		vector<ColumnVector> Sg(ngs);
		for(int g = 1; g <= ngs; g++)
		  {
		    Yg[g-1].ReSize(design.getntptsingroup(g));
		    Yg[g-1] = 0;
		    Sg[g-1].ReSize(design.getntptsingroup(g));
		    Sg[g-1] = 0; 
		    int t2=1;
		    for(int t = 1; t <= ntpts; t++)
		      {
			if(design.getgroup(t)==g)
			  {
			    Yg[g-1](t2) = Y(t);
			    Sg[g-1](t2) = S(t);
			    t2++;
			  }
		      }
//  		    OUT(Yg[g-1].t());
//  		    OUT(Sg[g-1].t());
		  }
		
		// initialize using ols
		// initialize gam
		ColumnVector gam = pinvdm*Y;
//  		OUT(gam.t());
		
		// initialize beta
		ColumnVector beta(ngs);
		beta = 0;
		for(int g = 1; g <= ngs; g++)
		  {
		    ColumnVector r = Yg[g-1]-zg[g-1]*gam;
		    float r2 = (r.t()*r).AsScalar();
		    float tmp = r2/(traceRg[g-1])-mean(Sg[g-1]).AsScalar();
		    
		    beta(g) = Max(0.01,tmp);
		  }
//  		OUT(beta.t());

		ColumnVector m(ntpts);
		m = 0;		
		ColumnVector v(ntpts);
		v = 0;
	
		// iterate em
		for(int i = 0; i < 100; i++)
		  {
		    // e step
		    for(int t = 1; t <= ntpts; t++)
		      {			    
			m(t) = (Y(t)*beta(design.getgroup(t)) + S(t)*(dm.Row(t)*gam).AsScalar())/(beta(design.getgroup(t)) + S(t));
			v(t) = (beta(design.getgroup(t))*S(t))/(beta(design.getgroup(t)) + S(t));
		      }
		    
		    // m step for gam
		    gam=pinvdm*m;
		    
		    // m step for beta
		    for(int g = 1; g <= ngs; g++)
		      {			
			ColumnVector mg(design.getntptsingroup(g));
			mg = 0;
			ColumnVector vg(design.getntptsingroup(g));
			vg = 0;

			int t2=1;
			for(int t = 1; t <= ntpts; t++)
			  {
			    if(design.getgroup(t)==g)
			      {
				mg(t2) = m(t);
				vg(t2) = v(t);
				t2++;
			      }
			  }
			beta(g)=(sum(vg)+mg.t()*mg-gam.t()*zzg[g-1]*gam).AsScalar()/(traceRg[g-1]);
		      }
		  }
		
//  		OUT(gam.t());
//  		OUT(beta.t());

		// for MCMC:	       
		for(int g = 0; g < ngs; g++)
		  {
		    float betamean = beta(g+1);
		    float betavar = abs(betamean)*2;
		    beta_b[g](vox) = Sqr(betamean)/betavar;
		    beta_c[g](vox) = betamean/betavar;
		  }

		// store results for gam:
		for(int e = 0; e < nevs; e++)
		  {
		    pes[e](vox) = gam(e+1);			
		  }

		// store results for beta:
		for(int g = 0; g < ngs; g++)
		  {
		    beta_mean[g](vox) = beta(g+1);
		  }
		
		// store results for cov		
		DiagonalMatrix V(ntpts);
		for(int t = 1; t <= ntpts; t++)
		  {
		    V(t) = S(t)+beta(design.getgroup(t));
		  }

		Matrix gamcovariance = (dm.t()*V.i()*dm).i();
// 		OUT(gamcovariance);
		ColumnVector gamcovariance2;
		reshape(gamcovariance2, gamcovariance, nevs*nevs, 1);		
		cov_pes.setSeries(gamcovariance2,vox);	        

		em_contrasts(gam,gamcovariance,vox);  
		if(GsOptions::getInstance().zlowerthreshold.value() > 0)
		  {
		    // check variance ratio
		    ColumnVector vr(ngs);
		    for(int g = 0; g < ngs; g++)
		      {
			vr(g+1) = beta(g+1)/mean(Sg[g]).AsScalar();
		      }		    
		    bool likeols = min(vr).AsScalar()>10; 
		    
		    if(likeols)
		      {
			// don't do MCMC if variance ratio large
			// assume OLS DOF
			zmask(x,y,z) = false;			
		      }
		    else
		      {			
			zmask(x,y,z) = pass_through_to_mcmc(GsOptions::getInstance().zlowerthreshold.value(),GsOptions::getInstance().zupperthreshold.value(),vox);
		      }
		  }
	      }
	  }
    nmaskvoxels = int(zmask.sum());    
  }

  void Gsmanager::ols_contrasts(const ColumnVector& mn, const Matrix& covariance, int vox)
  {
    Tracer_Plus trace("Gsmanager::ols_contrasts");    
    
    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	RowVector tcon = design.gettcontrast(t+1);
	t_ols_contrast(mn,covariance,tcon,tcopes[t](vox), tvarcopes[t](vox), ts[t](vox), tdofs[t](vox), zts[t](vox), false, vox);
	zolsts[t](vox) = zts[t](vox);
      }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](vox), fdof1s[f](vox), fdof2s[f](vox), zfs[f](vox), false, vox);
	zolsfs[f](vox) = zfs[f](vox);
      }    
  }

  void Gsmanager::em_contrasts(const ColumnVector& mn, const Matrix& covariance, int vox)
  {
    Tracer_Plus trace("Gsmanager::em_contrasts");    
    
    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	RowVector tcon = design.gettcontrast(t+1);
	
	double tdofupper;
	double tdoflower;

	// call first with highest possible DOF
	t_ols_contrast(mn,covariance,tcon,tcopes[t](vox), tvarcopes[t](vox), ts[t](vox), tdofupper, zemupperts[t](vox), true, vox);
	
	// call with OLS DOF
	if(!opts.fixed.value())
	  t_ols_contrast(mn,covariance,tcon,tcopes[t](vox), tvarcopes[t](vox), ts[t](vox), tdoflower, zemlowerts[t](vox), false, vox);

	if(opts.fixed.value())
	  {
	    tdofs[t](vox) = tdofupper;
	    zts[t](vox) = zemupperts[t](vox);
	  }
	else
	  {
	    tdofs[t](vox) = tdoflower;
	    // set z-score to lower z bound	
	    zts[t](vox) = zemlowerts[t](vox);
	  }
	
      }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	double fdof2upper;
	double fdof2lower;

	// call first with highest possible DOF
	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](vox), fdof1s[f](vox), fdof2upper, zemupperfs[f](vox), true, vox);

	// call with OLS DOF
	f_ols_contrast(mn,covariance,design.getfcontrast(f+1),fs[f](vox), fdof1s[f](vox), fdof2lower, zemlowerfs[f](vox), false, vox);

	if(opts.fixed.value())
	  {
	    fdof2s[f](vox) = fdof2upper;
	    zfs[f](vox) = zemupperfs[f](vox);
	  }
	else
	  {
	    fdof2s[f](vox) = fdof2lower;
	    // set z-score to lower z bound
	    zfs[f](vox) = zemlowerfs[f](vox);
	  }
      }
  }

  bool Gsmanager::pass_through_to_mcmc(float zlowerthresh, float zupperthresh, int vox)
  {
    Tracer_Plus trace("Gsmanager::pass_through_to_mcmc");   
    
    bool ret = false;
    
    // both emupper and emlower need to be jointly above
    // or below the threshold region to avoid the need 
    // for MCMC
    for(int t = 0; !ret && t < design.getnumtcontrasts(); t++)
      {		
	ret = !(((zemupperts[t](vox)) > (zupperthresh) &&
		  (zemlowerts[t](vox)) > (zupperthresh)) ||
	  ((zemupperts[t](vox)) < (zlowerthresh) &&
	   (zemlowerts[t](vox)) < (zlowerthresh)));
      }

    for(int f = 0; !ret && f < design.getnumfcontrasts(); f++)
      {
	ret = !(((zemupperfs[f](vox)) > (zupperthresh) &&
		  (zemlowerfs[f](vox)) > (zupperthresh)) ||
		((zemupperfs[f](vox)) < (zlowerthresh) &&
		 (zemlowerfs[f](vox)) < (zlowerthresh)));
      }

    return ret;
  }

  void Gsmanager::t_ols_contrast(const ColumnVector& mn, const Matrix& covariance, const RowVector& tcontrast, double& cope, double& varcope, double& t, double& dof, double& z, bool lookupdof, int vox)
  {
    Tracer_Plus trace("Gsmanager::t_ols_contrast");

    varcope = (tcontrast*covariance*tcontrast.t()).AsScalar();
    cope = (tcontrast*mn).AsScalar();

    if(lookupdof)
      {
	if(opts.fixed.value() && dofpassedin)
	  {
	    dof = dofvarcopedata.getSeries(vox).Sum() - nevs;
//  	    OUT(dofvarcopedata.getSeries(vox).Sum());
//  	    OUT(nevs);
//  	    OUT(dofvarcopedata.getSeries(vox));
//  	    OUT(dof);
	  }
	
	else
	  //dof = doflookuptable(ntpts - nevs);	
	  dof = 1000;
      }
    else
      {
	dof = float(ntpts - nevs);
      }

//      OUT(dof);
//      OUT(varcope);
//      OUT(cope);
    t = cope/sqrt(varcope);
 //      OUT(t);
 
    z = T2z::getInstance().convert(t,int(dof));
//      OUT(z);
  }

  void Gsmanager::f_ols_contrast(const ColumnVector& mn, const Matrix& covariance, const Matrix& fcontrast, double& f, double& dof1, double& dof2, double& z, bool lookupdof, int vox)
  {
    Tracer_Plus trace("Gsmanager::f_ols_contrast");

    dof1 = float(fcontrast.Nrows());

    if(lookupdof)
      {
	if(opts.fixed.value() && dofpassedin)
	  dof2 = dofvarcopedata.getSeries(vox).Sum() - nevs;
	
	else
	  //dof = doflookuptable(ntpts - nevs);	
	  dof2 = 1000;
      }
    else
      {
	dof2 = float(ntpts - nevs);
      }
    
    f = (mn.t()*fcontrast.t()*(fcontrast*covariance*fcontrast.t()).i()*fcontrast*mn/dof1).AsScalar();

    z = F2z::getInstance().convert(f,int(dof1),int(dof2));
  }
 
  void Gsmanager::mcmc_contrasts(const Matrix& gamsamples, int vox)
  {
    Tracer_Plus trace("Gsmanager::mcmc_contrasts");   
    
    for(int t = 0; t < design.getnumtcontrasts(); t++)
      {
	RowVector tcon = design.gettcontrast(t+1);
	t_mcmc_contrast(gamsamples, tcon, tcopes[t](vox), tvarcopes[t](vox), ts[t](vox), tdofs[t](vox), zts[t](vox), vox);		    
      }

    for(int f = 0; f < design.getnumfcontrasts(); f++)
      {
	f_mcmc_contrast(gamsamples, design.getfcontrast(f+1), fs[f](vox), fdof1s[f](vox), fdof2s[f](vox), zfs[f](vox), vox);		    
      }
       
  }

  void Gsmanager::t_mcmc_contrast(const Matrix& gamsamples, const RowVector& tcontrast, double& cope, double& varcope, double& t, double& dof, double& z, int vox)
  {
    Tracer_Plus trace("Gsmanager::t_mcmc_contrast");
 
    //gamsamples(nevs, nsamples);
  
    ColumnVector m;
    SymmetricMatrix covar;

    Matrix tcsamples = tcontrast*gamsamples;

    float tmpdof;

    ColumnVector gammean(nevs);gammean=0;
    for(int e = 0; e < nevs; e++)
      gammean(e+1) = pes[e](vox);

    m = tcontrast*gammean;
    multitfit(tcsamples, m, covar, tmpdof, true); 
       
    dof = double(tmpdof);
 
    varcope = covar(1,1);
    cope = m(1);

    t = cope/sqrt(varcope);
    z = T2z::getInstance().convert(t,int(dof));
  }

  void Gsmanager::f_mcmc_contrast(const Matrix& gamsamples, const Matrix& fcontrast, double& f, double& dof1, double& dof2, double& z, int vox)
  {
    Tracer_Plus trace("Gsmanager::f_mcmc_contrast");
    
    //gamsamples(nevs, nsamples);
    
    ColumnVector m;
    SymmetricMatrix covar;
    
    float tmpdof2;
    ColumnVector gammean(nevs);gammean=0;
    for(int e = 0; e < nevs; e++)
      gammean(e+1) = pes[e](vox);

    m = gammean;
    multitfit(gamsamples, m, covar, tmpdof2, true);

    dof2 = double(tmpdof2);
    dof1 = float(fcontrast.Nrows());

    f = (m.t()*fcontrast.t()*(fcontrast*covar*fcontrast.t()).i()*fcontrast*m/dof1).AsScalar();

    z = F2z::getInstance().convert(f,int(dof1),int(dof2));
  }

//   void Gsmanager::variational_bayes(int vox, ColumnVector& gammeanout, Matrix& gamcovout)
//   {
//     Tracer_Plus trace("Gsmanager::variational_bayes");
    
//     const Matrix& dm = design.getdm();

//     ColumnVector beta_b(ngs);
//     ColumnVector beta_c(ngs);
//     ColumnVector Y(ntpts+nevs);
//     ColumnVector m(ntpts+nevs);
//     SymmetricMatrix Rinv(ntpts+nevs);
//     SymmetricMatrix R(ntpts+nevs);
//     DiagonalMatrix Q(ntpts+nevs);
//     SymmetricMatrix S(ntpts+nevs);
//     beta_b = 1;
//     beta_c = 1;
//     Y = 0;
//     m = 0;
//     Rinv = 0;
//     R = 0;
//     Q = 0;
//     S = 0;

//     Y.Rows(1,ntpts) = copedata.getSeries(vox);
//     //		m = m_mean & gam_mean;

//     // Setup S
//     for(int t = 1; t <= ntpts; t++)
//       {
// 	S(t,t) = 1.0/varcopedata(t,vox);
		    
// 	for(int e = 1; e <= nevs; e++)
// 	  {
// 	    S(t,ntpts+e) = -dm(t,e)/varcopedata(t,vox);
// 	  }
//       }
	
//     SymmetricMatrix sumovert(nevs);
//     sumovert = 0.0;
//     for(int t = 1; t <= ntpts; t++)
//       {
// 	SymmetricMatrix temp;
// 	temp << dm.Row(t).t()*dm.Row(t)/varcopedata(t,vox);
// 	sumovert += temp;
//       }		    
//     S.SymSubMatrix(ntpts+1,ntpts+nevs) << sumovert;

//     for(int i =1; i <= 40; i++)
//       {		    		    
// 	// Update m and gamma
// 	for(int t = 1; t <= ntpts; t++)
// 	  {
	    
// 	    Q(t,t) = beta_b(design.getgroup(t))/beta_c(design.getgroup(t));
// 	  }
// 	R = S+Q;
// 	Rinv = R.i();
// 	m = Rinv*Q*Y;

// 	// Update beta
// 	for(int g = 1; g <= ngs; g++)
// 	  {
// 	    float sumovert = 0; 
// 	    for(int t = 1; t <= ntpts; t++)
// 	      {
// 		if(design.getgroup(t)==g)
// 		  {			    
// 		    sumovert += Sqr(m(t)-Y(t)) + Rinv(t,t);
// 		  }
// 	      } 
// 	    beta_b(g) = 1e-6+0.5*design.getntptsingroup(g);
// 	    beta_c(g) = 1e-6+0.5*sumovert;

// 	  }

//       }
		
//     gammeanout = m.Rows(ntpts+1,ntpts+nevs);

//     gamcovout = Rinv.SymSubMatrix(ntpts+1,ntpts+nevs); //gamcovout is COVARIANCE of gam    
	
//   }


}
