/*  macfnirt.cc

    Mark Jenkinson, FMRIB Image Analysis Group

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

// Motion Artefact Correction FNIRT
// Uses the Andersson/Ashburner-style framework

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "newimage/warpfns.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="macfnirt (Version 1.0)\nCopyright(c) 2005, University of Oxford (Mark Jenkinson)";
string examples="macfnirt [options] -i <input image> -o <output image>";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debug(string("--debug"), false, 
		     string("switch on debugging ouput"), 
		     false, no_argument);
Option<int> niter(string("-n,--niter"), 10,
		  string("number of main iterations (default=10)"),
		  false, requires_argument);
Option<int> nrefvol(string("--refvol"), 0,
		    string("number of reference volume (0 to N-1: default=N/2)"),
		    false, requires_argument);
Option<float> dtheta(string("--dtheta"), 0.5/180*M_PI,
		  string("delta theta in radians (default=0.5 degrees)"),
		  false, requires_argument);
Option<float> dtrans(string("--dtrans"), 0.5,
		  string("delta translation in mm (default=0.5mm)"),
		  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input filename"),
		  true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename"),
		  true, requires_argument);
Option<string> basiswarpname(string("-s,--shiftmap"), string(""),
		  string("input (B0) shiftmap filename"),
		  false, requires_argument);
int nonoptarg;


volume4D<float> basiswarp; // for convenience, make this global

////////////////////////////////////////////////////////////////////////////

// Local functions


int affine_transform(const volume<float>& vin, volume<float>& vout, 
		     const ColumnVector& params)
{
  Matrix aff(4,4);
  compose_aff(params,Min(6,params.Nrows()),vin.cog(),aff,construct_rotmat_euler);
  vout = vin; // needed to define reference volume size
  affine_transform(vin,vout,aff,1.0);
  return 0;
}


int transform(const volume<float>& vin, volume<float>& vout, 
	      const ColumnVector& params)
{
  //////////////////////////////////////////////////////////
  // WARNING:: THIS HAS NOT BEEN TESTED FOR WARPS YET ... //
  //////////////////////////////////////////////////////////
  if (params.Nrows()<=6) { return affine_transform(vin,vout,params); }
  // calculate affine component (also needed for realignment of B0 field)
  Matrix aff(4,4);
  compose_aff(params.Rows(1,6),6,vin.cog(),aff,construct_rotmat_euler);
  // set up deformation warp
  volume4D<float> warpvol;
  // initialise with identity warp and treat bases as relative warps (in mm)
  affine2warp(Identity(4),warpvol,vin);    // vin defines warp FOV and resolution
  volume<float> currwarp=vin*0.0;
  for (int n=7; n<=params.Nrows(); n++) {
    // transform B0 field by inverse of affine trans (undistorted -> distorted frame)
    for (int idx=0; idx<3; idx++) {
      affine_transform(basiswarp[3*(n-7)+idx],currwarp,aff.i(),1.0);
      warpvol[idx]+=((float) params(n))*currwarp;
    }
  }
  // now get affine warp component
  volume4D<float> affwarpvol, finalwarpvol;
  affine2warp(aff,affwarpvol,vin);  // need vin to define warp FOV and resolution

  // concatenate warps so the undistortion is done first
  concat_warps(warpvol,affwarpvol,finalwarpvol);
  
  if (debug.value()) { print_volume_info(warpvol); }
  // apply the required warp
  vout = vin;  // set up required FOV and resolution for output
  apply_warp(vin,vout,warpvol);
  return 0;
}


volume<float> param_derivative(volume4D<float>& dv, const volume<float>& vin, 
			       const ColumnVector& params, const ColumnVector &d_param) 
{
  int nparams=params.Nrows();
  if (dv.tsize()<nparams) {
    cerr << "ERROR:: too few volumes in dv passed into param_derivative" << endl;
    volume<float> dummy;
    return dummy;
  }
  float dparam=dtheta.value();
  ColumnVector newparams;

  // calculate volume at current param setting
  volume<float> v0;
  transform(vin,v0,params);

  // calculate perturbed volume and take differences to get derivative
  for (int n=1; n<=nparams; n++) {
    newparams = params;
    dparam = d_param(n);
    newparams(n) += dparam;
    transform(vin,dv[n-1],newparams);
    dv[n-1] -= v0;
    dv[n-1] *= (1.0/dparam);  // normalise to get correctly scaled derivative
  }
  return v0;
}


int do_work(int argc, char* argv[]) 
{
  volume4D<float> vin;
  read_volume4D(vin,inname.value());
  if (basiswarpname.set()) {
    volume<float> shiftmap;
    read_volume(shiftmap,basiswarpname.value());
    shift2warp(shiftmap,basiswarp,"y");
  } else {
    // construct a null basis (which should be ignored)
    basiswarp.addvolume(vin[0]*0.0f);
    basiswarp.addvolume(vin[0]*0.0f);
    basiswarp.addvolume(vin[0]*0.0f);
  }

  // set up refvol (n1)
  int n1, n2, realn1=0;
  n1 = vin.mint() - 1;
  if (nrefvol.set()) {
    n1 = nrefvol.value();
  }
  if ( (n1<vin.mint()) || (n1>vin.maxt()) ) {
    n1 = MISCMATHS::round(vin.tsize()/2.0);
  }
  if (n1!=0) {
    // big cheat to make bookkeeping easy - swap volumes n1 and 0
    volume<float> dummy;
    dummy = vin[0];
    vin[0] = vin[n1];
    vin[n1] = dummy;
    realn1 = n1;
    n1 = 0;
  }
  if (verbose.value()) { cout << "Refvol number = " << realn1 << endl; }

  // set up the initial parameters
  int nparams=6*(vin.maxt());  // 6 * one less than number of timepoints
  int nwarps=0;
  if (basiswarpname.set()) {
    nwarps = basiswarp.tsize() / 3;
    nparams += nwarps;
  }
  ColumnVector params(nparams);
  if (verbose.value()) { cout << "Number of warps = " << nwarps << " and number of params = " << nparams << endl; } 

  // set up output
  volume4D<float> vout;
  vout = vin;
  
  // set up required variables
  volume4D<float> dv;
  for (int n=1; n<=6+nwarps; n++) { dv.addvolume(0.0f*vin[n1]); }
  volume<float> vn2_0;
  Matrix xtx_all(nparams,nparams), xtx(6,6); 
  ColumnVector xty_all(nparams), dp(nparams), xty(6);
  xtx=0.0;  xty=0.0;
  float dp_thresh, dp_norm, d_theta, d_trans, tol_factor;
  ColumnVector param_tol(nparams), param_subset(6+nwarps), param_tol_subset(6+nwarps);
 
  params=0.0;  // dumb initialisation
  
  // set up parameter tolerances and thresholds (for convergence)
  dp_thresh = 0.01;
  tol_factor = 1.0;
  d_theta = dtheta.value();
  d_trans = dtrans.value();
  param_tol = 1.0;  // default for all warps
  for (n2=1; n2<=vin.maxt(); n2++) {
    // reset rigid body param tolerances (for rotation and translation separately)
    param_tol(n2*6-5)=d_theta; param_tol(n2*6-4)=d_theta; param_tol(n2*6-3)=d_theta; 
    param_tol(n2*6-2)=d_trans; param_tol(n2*6-1)=d_trans; param_tol(n2*6)=d_trans; 
  }
  dp_norm=10*dp_thresh;

  // loop with convergence criterion (first attempt at criterion)
  int iter=1;
  while ((iter<=niter.value()) && (dp_norm>dp_thresh)) {
    
    // generate the required matrices (including timepoints and voxels)
    for (n2=1; n2<=vin.maxt(); n2++) {
      // generate the derivative images
      param_subset = params.Rows(n2*6-5,n2*6);
      param_tol_subset = param_tol.Rows(n2*6-5,n2*6);
      if (nwarps>0) {
	param_subset = param_subset & params.Rows(vin.maxt()*6 + 1,params.Nrows());
	param_tol_subset = param_tol_subset & param_tol.Rows(vin.maxt()*6 + 1,params.Nrows());
      }
      vn2_0=param_derivative(dv,vin[n2],param_subset,param_tol_subset*tol_factor);
      if (debug.value() && (iter==1)) { save_volume4D(dv,outname.value()+"_deriv"); }
      
      // generate the required dot products
      // i.e. dv.t() * (T(vin[n2])-vin[n1])  and (dv.t()*dv).i()

      volume<float> dummy, vn1_0;
      // calculate the current state of the reference volume (only depends on warp)
      ColumnVector v1_params;
      v1_params = param_subset;
      v1_params.Rows(1,6) = v1_params.Rows(1,6) * 0;
      transform(vin[n1],vn1_0,v1_params);
      // now calculate the correlations between derivatives and signal (xtx and xty)
      for (int n=1; n<=6; n++) {
	dummy = dv[n-1] * ( vn2_0 - vn1_0 );
	xty(n) = dummy.sum() / dummy.nvoxels();
	for (int m=n; m<=6; m++) {
	  dummy = dv[n-1] * dv[m-1];
	  xtx(n,m) = dummy.sum() / dummy.nvoxels();
	  if (m>n) { xtx(m,n) = dummy.sum() / dummy.nvoxels(); }
	}
      }
      // put the relevant submatrix elements in place
      xtx_all.SubMatrix(n2*6-5,n2*6,n2*6-5,n2*6) = xtx.SubMatrix(1,6,1,6);
      xty_all.Rows(n2*6-5,n2*6) = xty.Rows(1,6);
      // calculate terms in interactions of warp and rigid body terms
      for (int n=1; n<=nwarps; n++) {
	for (int m=1; m<=6; m++) {
	  dummy = dv[5+n] * dv[m-1];
	  xtx_all(n2*6-m+1,nparams-nwarps+n) += dummy.sum() / dummy.nvoxels();
	  xtx_all(nparams-nwarps+n,n2*6-m+1) += dummy.sum() / dummy.nvoxels();
	}
	dummy = dv[5+n] * ( vn2_0 - vn1_0 );
	xty_all(nparams-nwarps+n) += dummy.sum() / dummy.nvoxels();
	// diagonal terms (warp with warp on volume n2)
	for (int m=n; m<=nwarps; m++) {
	  dummy = dv[5+n] * dv[5+m];
	  xtx_all(nparams-nwarps+n,nparams-nwarps+m) += dummy.sum() / dummy.nvoxels();
	  xtx_all(nparams-nwarps+m,nparams-nwarps+n) += dummy.sum() / dummy.nvoxels();
	}
      }
    }
    
    // do the fit and get (all) the new params
    dp = xtx_all.i() * xty_all;
    params -= dp;
    
    if (debug.value()) {
      cout << "det(x.t()*x) = " << Determinant(xtx_all) << endl;
      cout << "Update on params: dp = " << -dp.t() << endl;
    }
    
    // calculate normalised parameter step
    dp_norm = 0.0;
    for (int nn=1; nn<=dp.Nrows(); nn++) { dp_norm += Sqr(dp(nn)/param_tol(nn)); }
    dp_norm = sqrt(dp_norm/dp.Nrows());
    // rescale tolerances to get a better match to the size of step taken
    if (dp_norm>1) { tol_factor *= Min(dp_norm,2.0); }
    if (dp_norm<1) { tol_factor *= Max(dp_norm,0.5); }

    if (verbose.value()) { 
      cout << "Params (iteration " << iter << ") = " << params.t() << endl;
      cout << "dp_norm = " << dp_norm << endl;
      cout << "tol_factor = " << tol_factor << endl;
    }
    if (debug.value()) {
      Matrix aff(4,4);
      compose_aff(params,Min(6,params.Nrows()),vin[n2].cog(),aff,construct_rotmat_euler);
      cout << "Affmat = " << endl << aff << endl;
    }
    
    iter++;  // end of while loop
  }
  
  
  // generate output volume
  for (n2=1; n2<=vin.maxt(); n2++) {
    // generate the derivative matrices
    param_subset = params.Rows(n2*6-5,n2*6);
    if (nwarps>0) {
      param_subset = param_subset & params.Rows(vin.maxt()*6 + 1,params.Nrows());
    }
    if (verbose.value()) { 
      cout << "Final parameters = " << param_subset.t() << endl;
      Matrix aff(4,4);
      compose_aff(param_subset,6,vin[n2].cog(),aff,construct_rotmat_euler);
      cout << "Final affine matrix = " << endl << aff << endl; 
    }
    transform(vin[n2],vout[n2],param_subset);
  }
  
  if (realn1!=0) {
    // big cheat - swap vols realn1 and 0 for an easy life
    volume<float> dummy;
    dummy = vout[realn1];
    vout[realn1] = vout[0];
    vout[0] = dummy;
  }

  // save output volume
  save_volume4D(vout,outname.value());
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(basiswarpname);
    options.add(niter);
    options.add(nrefvol);
    options.add(dtheta);
    options.add(dtrans);
    options.add(verbose);
    options.add(debug);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions

  return do_work(argc,argv);
}

