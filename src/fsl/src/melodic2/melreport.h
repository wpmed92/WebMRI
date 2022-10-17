/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melreport.h - report generation

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


#ifndef __MELODICREPORT_h
#define __MELODICREPORT_h
#include "newimage/newimageall.h"
#include "utils/log.h"
#include "melpca.h"
#include "meloptions.h"
#include "meldata.h"
#include "melgmix.h"
#include "melodic.h"
#include "newmatap.h"
#include "newmatio.h"
#include <time.h>
#include <sstream>
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"
#include "utils/options.h"

using namespace Utilities;
using namespace NEWIMAGE;
using namespace MISCPLOT;
using namespace MISCPIC;


namespace Melodic{
  
  class MelodicReport
    {
    public:
      MelodicReport(MelodicData &pmelodat, MelodicOptions &popts, Log &plogger):  
	melodat(pmelodat),
	opts(popts),
	logger(plogger)
	{
	  if( bool(opts.genreport.value()) ){
	  const time_t tmptime = time(NULL);
 
	  report.makeDir(logger.appendDir("report"),"00index.html");
	  report << "<HTML>" << endl << endl 
		 << "<TITLE>MELODIC Report</TITLE>"<< endl <<endl
		 << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
		 << "/doc/images/fsl-bg.jpg\">" << endl 
		 << endl << "<hr><CENTER> " << endl 
		 << "<H1>MELODIC Report</H1>" 
		 << endl
		 << report.getDir() << "/" << report.getLogFileName() 
		 << "<br>" << endl
		 << ctime(&tmptime) << "<br>" << endl;
	  /* ifstream reptest(string("../report.log").c_str()); */
/* 	  if(!reptest) */
/* 	    report << "<A HREF=\"/" << logger.appendDir("report.log")  */
/* 		   << "\">report.log</A><br>" << endl; */
/* 	  else */
/* 	    report << "<A HREF=\"" << logger.appendDir("melodic.log")  */
/* 		   << "\">melodic.log</A><br>" << endl; */
/* 	  reptest.close(); */
	  report << "</CENTER>" << endl << endl
		 << "<hr><H2>Components:</H2> <p>" << endl << endl;
	  }
	}

      ~MelodicReport(){
	if( bool(opts.genreport.value()) ){
	  report << "<HR><FONT SIZE=1>This page produced automatically by "
		 << "<A HREF=\"http://www.fmrib.ox.ac.uk/fsl/melodic/index.html\"> MELODIC </A>" 
		 << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
		 << "FMRIB Software Library</A>.</FONT>" << endl
		 << "</BODY></HTML>" << endl;
	}
      }

      inline void analysistxt(){
	if( bool(opts.genreport.value()) ){

	  report << "<hr><h2>Analysis methods</h2> <p>"<<endl;
	  report << "Analysis was carried out using MELODIC (Multivariate Exploratory Linear Decomposition into Independent Components) Version ln(11), part of FSL (FMRIB's Software Library, <A HREF=\"http://www.fmrib.ox.ac.uk/fsl/\">www.fmrib.ox.ac.uk/fsl</A>), an implementation for the estimation of a Probabilistic Independent Component Analysis model [Beckmann 2004]."<<endl;
	  
	  report << "The following melodic pre-processing was applied to the input data file: "<< endl;

	  if(opts.use_mask.value())
	    report << " masking of non-brain voxels;";
	  
	  report << " voxel-wise de-meaning of the data;" << endl;
	  
	  if(opts.varnorm.value())
	    report << " normalisation of the voxel-wise variance; ";
	  
	  report << "<br>"<<endl;
	  
	  report << " Pre-processed data was whitened and projected into a " 
		 << melodat.get_mix().Ncols()<< "-dimensional subspace using ";
	  if(melodat.get_PPCA().Storage()>0){	    
	    report << "probabilistic Principal Component Analysis where the number of dimensions was estimated using ";
	    if(opts.pca_est.value() == string("lap"))
	      report << "the Laplace approximation to the Bayesian evidence of the model order [Minka 2000, Beckmann 2004]. " << endl;
	    else
	      if(opts.pca_est.value() == string("bic"))
		report << "the <em> Bayesian Information Criterion</em> (BIC) [Kass 1993]. " << endl;
	      else
		if(opts.pca_est.value() == string("mdl"))
		  report << "<em> Minimum Description Length</em> (MDL) [Rissanen 1978]. " << endl;
		else
		  if(opts.pca_est.value() == string("aic"))
		    report << "the <em> Akaike Information Criterion</em> (AIC) [Akaike 1969]. " << endl;
		  else
		    report << "a committee of approximations to Bayesian the model order [Beckmann 2004]. " << endl;

	  }	  
	  else
	    report << "Principal Component Analysis. ";
	  
	  
	  report << " The whitened observations were decomposed into a set of time-courses and spatial maps by optimising for non-Gaussian spatial source distributions using a fixed-point iteration technique [Hyv&auml;rinen 1999]. " << endl;
	  
	  report << "Estimated Component maps were divided by the standard deviation of the residual noise";
	  
	  if(opts.perf_mm.value())
	    report << " and thresholded by fitting a mixture model to the histogram of intensity values [Beckmann 2004]. <p>" << endl;
	  else
	    report <<".<p>" << endl;
	 
	  refstxt(); 
	}
      }

      inline void refstxt(){

	if( bool(opts.genreport.value()) ){
	  report << "<h3>References</h3> <p>"<<endl;
	  
	  report << "[Hyv&auml;rinen 1999] A. Hyv&auml;rinen. Fast and Robust Fixed-Point Algorithms for Independent Component Analysis. IEEE Transactions on Neural Networks 10(3):626-634, 1999.<br> " << endl;
 report << "[Beckmann 2004] C.F. Beckmann and S.M. Smith. Probabilistic Independent Component Analysis for Functional Magnetic Resonance Imaging. IEEE Transactions on Medical Imaging 23(2):137-152 2004. <br>" << endl;

	  /* if(opts.perf_mm.value()){ */
/* 	    report << "[Bullmore 1996] E. Bullmore <em>et. al.</em> Statistical methods of estimation and inference for functional MR image analysis. Magnetic Resonance in Medicine, 35(2):261-177, 1996. <br>" << endl; */
	   
/* 	  } */

	  if(melodat.get_PPCA().Storage()>0){	    
	    report << "[Everson 2000] R. Everson and S. Roberts. Inferring the eigenvalues of covariance matrices from limited, noisy data. IEEE Trans Signal Processing, 48(7):2083-2091, 2000<br>"<<endl;
	    report << "[Tipping 1999] M.E. Tipping and C.M.Bishop. Probabilistic Principal component analysis. J Royal Statistical Society B, 61(3), 1999. <br>" << endl;
	    /*  report << "[Beckmann 2001] C.F. Beckmann, J.A. Noble and S.M. Smith. Investigating the intrinsic dimensionality of FMRI data for ICA. In Seventh Int. Conf. on Functional Mapping of the Human Brain, 2001. <br>" << endl;*/
	    if(opts.pca_est.value() == string("lap"))
	      report << "[Minka 2000] T. Minka. Automatic choice of dimensionality for PCA. Technical Report 514, MIT Media Lab Vision and Modeling Group, 2000. <BR>"<< endl;
	    else
	      if(opts.pca_est.value() == string("bic"))
		report << "[Kass 1995] R.E. Kass and A. E. Raftery. Bayes factors. Journal of the American Statistical Association, 90:733-795, 1995 <br>" << endl;
	      else
		if(opts.pca_est.value() == string("mdl"))
		  report << "[Rissanen 1978]. J. Rissanen. Modelling by shortest data description. Automatica, 14:465-471, 1978. <br>" << endl;
		else
		  if(opts.pca_est.value() == string("aic"))
		    report << "[Akaike 1974]. H. Akaike. A new look at statistical model identification. IEEE Transactions on Automatic Control, 19:716-723, 1974. <br>" << endl;
		  else
		    report << "[Minka 2000]. T. Minka. Automatic choice of dimensionality for PCA. Technical Report 514, MIT Media Lab Vision and Modeling Group, 2000. <BR>" << endl;

	}
	}
      }

      inline void addtxt(string what){
	if( bool(opts.genreport.value()) ){
	  report << what << endl;
	}
      }
      
      inline void addpar(string what){
	if( bool(opts.genreport.value()) ){
	  report << "<p>" << what << endl;
	}
      }
      
      inline void addlink(string where, string what){
	if( bool(opts.genreport.value()) ){
	  report << "<A HREF=\"" << where << "\"> " << what << "</A> ";
	}
      }

      inline void addpic(string what, string link = ""){
	if( bool(opts.genreport.value()) ){
	  if( link.length() > 0)
	    report << "<A HREF=\"" << link << "\"> ";

	  report << "<img BORDER=0 SRC=\"" << what<< ".png\"><p>";
	  if( link.length() > 0)
	    report << "</A> ";
	}
      }

      inline string getDir(){
	return report.getDir();
      }

      void IC_rep(MelGMix &mmodel, int cnum, int dim, Matrix ICstats);
      void IC_simplerep(string prefix, int cnum, int dim);

      void PPCA_rep();

    private:
      MelodicData &melodat;
      MelodicOptions &opts;
      Log &logger;
      Log report;

      Log IChtml;
      Log IChtml2;

      void IC_rep_det(MelGMix &mmodel, int cnum, int dim);
      
      string int2str(int n)
	{
	  ostringstream os;
	  //    os.fill(' ');
	  //    os.width(width);
	  os.setf(ios::internal, ios::adjustfield);
	  os << n;
	  return os.str();
	}
      

      string float2str(float f, int width, int prec, int scientif)
	{
	    ostringstream os;
	    int redw = int(std::abs(std::log10(std::abs(f))))+1;
	    if(width>0)
	      os.width(width);
	    if(scientif>0)
	      os.setf(ios::scientific);
	    os.precision(redw+std::abs(prec));
	    os.setf(ios::internal, ios::adjustfield);
	    os << f;
	    return os.str();
	}
    };   

}
#endif
