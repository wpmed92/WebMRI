/*  groupttest.cc

    Mark Jenkinson, FMRIB Image Analysis Group
    Ana Juric, Mental Health Research Institute,
       Centre for Neuroscience, University of Melbourne

    Copyright (C) 2004 University of Oxford  */

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

// Calculates the surface normals for a mask, using a smoothed
//  gradient calculation (all non-surface points get zero ouput)

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include <vector>
#include <algorithm>
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"
#include "miscmaths/t2z.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="groupttest (Version 1.1)\nCopyright(c) 2004, University of Oxford (Mark Jenkinson)";
string examples="groupttest --na=<number in group A> --nb=<number in group B> -m <maskvol> -o <groupres> [options] <list of images for group A> <list of images for group B>\ne.g.   groupttest --na=15 --nb=15 -m maskvol -o groupres groupA/*.hdr* groupB/*.hdr*";

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
Option<bool> conservativetest(string("--conservative"), false,
			      string("use conservative FDR correction factor"),
			      false, no_argument);
Option<int>  numa(string("--na"),0,
		  string("number of members of group A (normals)"),
		  true, requires_argument);
Option<int>  numb(string("--nb"),0,
		  string("number of members of group B (patients)"),
		  true, requires_argument);
Option<string> ordername(string("--order"), string(""),
		      string("~\toutput image of order values"),
		      false, requires_argument);
Option<string> maskname(string("-m"), string(""),
		      string("input mask filename"),
		      true, requires_argument);
Option<string> outname(string("-o"), string(""),
		       string("output base filename"),
		       true, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Support functions

int save_as_image(const string& filename, const volume<float>& mask, 
		  const Matrix& valmat)
{
    // put values back into volume format
    if (verbose.value()) { cerr << "Saving results to " << filename << endl; }
    volume4D<float> outvals;
    outvals.addvolume(mask);
    outvals.setmatrix(valmat.t(),mask);
    return save_volume4D(outvals,filename);
}

Matrix get_coord_matrix(const volume<float>& mask)
{
  // construct a matrix of index values 1 -> Ntot
  volume4D<float> outvals;
  outvals.addvolume(mask);
  Matrix index = outvals.matrix(mask);
  int Ntot = index.Ncols();
  for (int j=1; j<=Ntot; j++) {
    index(1,j) = j;
  }
  outvals.setmatrix(index,mask);
  // go through volume and set up new matrix *with coordinates*
  Matrix coords(Ntot,3);
  for (int z=mask.minz(); z<=mask.maxz(); z++) {
    for (int y=mask.miny(); y<=mask.maxy(); y++) {
      for (int x=mask.minx(); x<=mask.maxx(); x++) {
	if (mask(x,y,z)>0.5) {
	  int idx = MISCMATHS::round(outvals(x,y,z,0));
	  coords(idx,1) = x;
	  coords(idx,2) = y;
	  coords(idx,3) = z;
	}
      }
    }
  }
  return coords;
}


volume<float> calc_edge_mask(const volume<float>& vmask)
{
  volume<float> vtmp = vmask;
  bool atedge;
  if (verbose.value()) { cerr << "Extracting Edge Voxels" << endl; }
  for (int z=vmask.minz(); z<=vmask.maxz(); z++) {
    for (int y=vmask.miny(); y<=vmask.maxy(); y++) {
      for (int x=vmask.minx(); x<=vmask.maxx(); x++) {
	atedge = false;
	if ( (vmask(x,y,z)>0.5) ) {
	  if (vmask(x,y,z-1)<0.5) atedge=true;
	  else { 
	    if (vmask(x,y-1,z)<0.5) atedge=true;
	    else { 
	      if (vmask(x-1,y,z)<0.5) atedge=true;
	      else {
		if (vmask(x+1,y,z)<0.5) atedge=true;
		else {
		  if (vmask(x,y+1,z)<0.5) atedge=true;
		  else {
		    if (vmask(x,y,z+1)<0.5) atedge=true;
		  }
		}
	      }
	    }
	  }
	}
	if (atedge) {
	  vtmp(x,y,z)=1;
	} else {
	  vtmp(x,y,z)=0;
	}
      }
    }
  }
  return vtmp;
}

 
//Function written by Ana Juric 
// Gentleman and Jenkins approximation for the t-distribution p-values (Biometrika, 55(3), p 571, 1968)
//  NB: gives coefficents (c1,..,c5) for:
//    p(|t|<X) = 1 - (c5*X^5 + c4*X^4 + c3*X^3 + c2*X^2 + c1*X + 1)^(-8)

double tTesting(double degreesOfFreedom, int coefficientNum)
{
        double coefficientMatrix[5][7]={
                {0.09979441, -0.5818210, 1.390993, -1.222452, 2.151185, -5.537409, 11.42343},
                {0.04431742,-0.2206018, -0.03317253, 5.679969, -12.96519, -5.166733, 13.49862},
                {0.009694901, -0.1408854, 1.889930, -12.75532, 25.77532, -4.233736, 14.39630},
                {-0.00009187228, 0.03789901, -1.280346, 9.249528, -19.08115, -2.777816, 16.46132},
                {0.0005796020, -0.02763334, 0.4517029, -2.657697, 5.127212, -0.5657187, 21.83269} };
 
        double coefficient;
 
	double v=degreesOfFreedom;
	double c6, c5, c4, c3, c2, c1, c0;
	c6 = coefficientMatrix[coefficientNum][6];
	c5 = coefficientMatrix[coefficientNum][5];
	c4 = coefficientMatrix[coefficientNum][4];
	c3 = coefficientMatrix[coefficientNum][3];
	c2 = coefficientMatrix[coefficientNum][2];
	c1 = coefficientMatrix[coefficientNum][1];
	c0 = coefficientMatrix[coefficientNum][0];

	// old version
	/*
        coefficient=
                (((coefficientMatrix[coefficientNum][4]*(pow(degreesOfFreedom,(-4))))+
                  (coefficientMatrix[coefficientNum][3]*(pow(degreesOfFreedom,(-3))))+
                  (coefficientMatrix[coefficientNum][2]*(pow(degreesOfFreedom,(-2))))+
                  (coefficientMatrix[coefficientNum][1]*(pow(degreesOfFreedom,(-1))))+
                  (coefficientMatrix[coefficientNum][0]))
                /((coefficientMatrix[coefficientNum][6]*(pow(degreesOfFreedom,(-2))))+
                  (coefficientMatrix[coefficientNum][5]*(pow(degreesOfFreedom,(-1))))+1));
	*/ 

	// new version - note that both denom & numerator are multiplied by v^4 in order to
	//  have positive powers of v only (not v^(-4), etc.)
	coefficient = (c4 + v*(c3 + v*(c2 + v*(c1 + v*c0)))) 
	  / (v*v*(c6 + v*(c5 + v)));
        return coefficient;
}


double pvalue(double tX, double dof) {
  // return the ONE SIDED t-test p-values: p(t>X)
  //    based on the two-sided formula:
  //    p(|t|>X) = (c5*X^5 + c4*X^4 + c3*X^3 + c2*X^2 + c1*X + 1)^(-8)
  double p1, p, x;
  // Code fragment by Ana Juric
  // "initialises the coeffMatrix with the relevent values"
  double c[5];
  for(int anaj=0; anaj<5; anaj++) { c[anaj]=tTesting(dof,anaj); } 

  x = fabs(tX);
  p =  pow((1 + x*(c[0] + x*(c[1] + x*(c[2] + x*(c[3] + x*c[4]))))),-8.0);
  if (tX>0) {
    p1 = p/2.0;
  } else {
    p1 = 1 - p/2.0;
  }
  return p1;
} 



vector<int> get_sortindex(const Matrix& vals)
{
  // return the mapping of old indices to new indices in the
  //   new *ascending* sort of vals
  int length=vals.Nrows();
  vector<pair<double, int> > sortlist(length);
  for (int n=0; n<length; n++) {
    sortlist[n] = pair<double, int>((double) vals(n+1,1),n+1);
  }
  sort(sortlist.begin(),sortlist.end());  // O(N.log(N))
  vector<int> idx(length);
  for (int n=0; n<length; n++) {
    idx[sortlist[n].second-1] = n+1;
  }
  return idx;
}

////////////////////////////////////////////////////////////////////////////

// Main function - this does all the work

int do_work(int argc, char* argv[], int nonoptarg) 
{
  string basename = fslbasename(outname.value());

  volume<float> vmask, vtmp;
  read_volume(vmask,maskname.value());
  if (verbose.value()) print_info(vmask,"vmask");
  vmask = calc_edge_mask(vmask);

  // get ready to read in flow images
  int Ntot = MISCMATHS::round(vmask.sum());
  int N1=numa.value();
  int N2=numb.value();
  if (verbose.value()) { cerr << "Ntot = " << Ntot << " ; Na,Nb = " << N1 << " , " << N2 << endl; }
  Matrix newcol(Ntot,1), bigmatrix(Ntot,N1+N2), tvalmat(Ntot,1);
  Matrix pmat(Ntot,1), logqmat(Ntot,1);
  Matrix meana(Ntot,1), meanb(Ntot,1);

  // read in images and accumulate values into bigmatrix
  for (int n=1; n<=(N1+N2); n++) {
    volume4D<float> vstat;
    string filename = argv[nonoptarg + n - 1];
    if (verbose.value()) { cerr << "Reading file " << filename << endl; }
    read_volume4D(vstat,filename);
    if (verbose.value()) { print_info(vstat,"vstat"); }
    newcol = vstat.matrix(vmask);
    for (int j=1; j<=Ntot; j++) {
      bigmatrix(j,n) = newcol(1,j);
    }
  }

  // Calculate t-values and relevant means and standard deviations
  for (int j=1; j<=Ntot; j++) {
    double sumx1=0, sumx2=0, sumx1sq=0, sumx2sq=0, meanx1=0, meanx2=0, sdx1=0, sdx2=0;
    for (int m=1; m<=N1; m++) {
      double x1 = bigmatrix(j,m);
      sumx1 += x1;
      sumx1sq += x1*x1;
    }
    for (int m=N1+1; m<=(N1+N2); m++) {
      double x2 = bigmatrix(j,m);
      sumx2 += x2;
      sumx2sq += x2*x2;
    }
    meanx1 = sumx1 / N1;
    meanx2 = sumx2 / N2;
    meana(j,1) = meanx1;
    meanb(j,1) = meanx2;
    sdx1 = sqrt(sumx1sq/N1 - meanx1*meanx1);
    sdx2 = sqrt(sumx2sq/N2 - meanx2*meanx2);
    double sediff = sqrt(sdx1*sdx1 / N1 + sdx2*sdx2 / N2);
    double tval = (meanx1 - meanx2)/sediff;
    // Welch's degree of freedom (for unequal variances)
    double dof = Sqr(sdx1*sdx1/N1 + sdx2*sdx2/N2) / 
      ( Sqr(sdx1*sdx1/N1)/(N1-1) + Sqr(sdx2*sdx2/N2)/(N2-1) );
    // store t value
    tvalmat(j,1) = tval;
    // convert t values to p values
    pmat(j,1) = pvalue(tval,dof);
  }

  // save the image results
  save_as_image(basename+"_meanA",vmask,meana);
  save_as_image(basename+"_meanB",vmask,meanb);
  save_as_image(basename+"_tvals",vmask,tvalmat);
  save_as_image(basename+"_pvals",vmask,pmat);
  // save the matrix result (coords + group means)
  {
    Matrix coords = get_coord_matrix(vmask);
    Matrix save_result = ( coords | meana ) | meanb ;
    write_ascii_matrix(save_result,basename+"_matrix");
  }

  // calculate FDR threshold required to make each voxel significant
  // FDR formula is: p = n*q / (N * C)  
  //   where n=order index, N=total number of p values, 
  //         C=1 for the simple case, and 
  //         C=1/1 + 1/2 + 1/3 + ... + 1/N for the most general correlation
  // We use the inverse formula: q_{min} = N*C*p / n

  if (verbose.value()) { cerr << "Calculating FDR values" << endl; } 
  float C=1.0;
  if (conservativetest.value()) {
    for (int n=2; n<=Ntot; n++) { C+=1.0/((double) n); }
  }

  vector<int> norder = get_sortindex(pmat);
  
  for (int j=1; j<=Ntot; j++) {
    double qval = pmat(j,1) * Ntot * C / norder[j-1];
    // qval isn't a probability - it can be greater than 1, but don't show these
    if (qval>1.0) qval=1.0;
    logqmat(j,1) = -log10(qval);
  }
  
  // save the FDR log(q_min) results
  save_as_image(basename+"_qvals",vmask,logqmat);

  // save the order values, if requested
  if (ordername.set()) {
    Matrix ordermat(Ntot,1);
    for (int j=1; j<=Ntot; j++) { ordermat(j,1) = norder[j-1]; }
    save_as_image(ordername.value(),vmask,ordermat);
  }

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
    options.add(maskname);
    options.add(outname);
    options.add(numa);
    options.add(numb);
    options.add(ordername);
    options.add(conservativetest);
    options.add(verbose);
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

  return do_work(argc,argv,nonoptarg);
}

