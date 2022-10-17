/*  film_gls_res.cc

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

#include <iostream>
#include <fstream>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "glim.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/volumeseries.h"
#include "miscmaths/volume.h"
#include "utils/log.h"
#include "AutoCorrEstimator.h"
#include "paradigm.h"
#include "FilmGlsOptionsRes.h"
#include "glimGls.h"
#include <vector>
#include <string>

using namespace NEWMAT;
using namespace FILM;
using namespace Utilities;

int main(int argc, char *argv[])
{
  try{
    rand();
    // parse command line to find out directory name for logging:
    ofstream out2;
    FilmGlsOptionsRes& globalopts = FilmGlsOptionsRes::getInstance();
    globalopts.parse_command_line(argc, argv, out2);

    // Setup logging:
    Log& logger = Log::getInstance();
    logger.makeDir(globalopts.datadir);

    // parse command line again to output arguments to logfile
    globalopts.parse_command_line(argc, argv, logger.str());

    // load non-temporally filtered data
    VolumeSeries x;
    x.read(globalopts.inputfname);

    // if needed output the 12th volume for use later
    Volume epivol;
    if(globalopts.smoothACEst)
      {
	epivol = x.getVolume(12).AsColumn();
	epivol.setInfo(x.getInfo());
	
	epivol.writeAsInt(logger.getDir() + "/" + globalopts.epifname);
      }

    // This also removes the mean from each of the time series:
    x.thresholdSeries(globalopts.thresh, true);

    // if needed later also threshold the epi volume
    if(globalopts.smoothACEst)
      {
	epivol.setPreThresholdPositions(x.getPreThresholdPositions());
	epivol.threshold();
      }

    int sizeTS = x.getNumVolumes();
    int numTS = x.getNumSeries();
    
    // Load paradigm: 
    Paradigm parad;
    parad.load(globalopts.paradigmfname, "", "", false, sizeTS);

    // Sort out detrending:
    if(globalopts.detrend)
      {
	MISCMATHS::detrend(x, 2);
      }

    if(globalopts.verbose)
      {
	logger.out("Gc", parad.getDesignMatrix());
      }

    vector<Matrix> dms;
    for(int i = 1; i <= numTS; i++)
      {
	dms.push_back(parad.getDesignMatrix());
      }
    
    // Setup OLS GLM for temporally filtered data:
    int numParams = parad.getDesignMatrix().Ncols();
    GlimGls glimGls(numTS, sizeTS, numParams);
    
    VolumeSeries residuals(sizeTS, numTS, x.getInfo(), x.getPreThresholdPositions());
    AutoCorrEstimator acEst(residuals);

    if(!globalopts.noest)
      {
	int numiters = globalopts.numiters+1;
	int iters = 1;

	// iters==1 is for high freq removal
	if(!globalopts.highfreqremoval) 
	  iters++;

	for(; iters <= numiters; iters++)
	  { 
	    cerr << "iters = " << iters << endl;
	    // Loop through voxels:
	    cerr << "Calculating residuals..." << endl; 
	    for(int i = 1; i <= numTS; i++)
	      {						    
		glimGls.setData(x.getSeries(i), dms[i-1], i);
		residuals.setSeries(glimGls.getResiduals());
	      }
	    cerr << "Completed" << endl; 

	    cerr << "Estimating residual autocorrelation..." << endl; 
	
	    // Estimate Autocorrelation:	    
	    acEst.calcRaw();

	    if(iters==1)
	      {
		// iters==1 is for high freq removal
		acEst.tukey(10);
	      }
	    else
	      {
		if(globalopts.fitAutoRegressiveModel)
		  {
		    acEst.fitAutoRegressiveModel();
		  }
		else if(globalopts.tukey)
		  {    
		    // Smooth raw estimates:
		    if(globalopts.smoothACEst)
		      {
			acEst.spatiallySmooth(logger.getDir() + "/" + globalopts.epifname, epivol, globalopts.ms, globalopts.epifname, globalopts.susanpath, globalopts.epith);		    
		      }
		    if(globalopts.tukeysize == 0)
		      globalopts.tukeysize = (int)(2*sqrt(sizeTS))/2;
		
		    acEst.tukey(globalopts.tukeysize);
		  }
		else if(globalopts.multitaper)
		  {
		    acEst.multitaper(globalopts.multitapersize);
		  }
		else if(globalopts.pava)
		  {
		    // Smooth raw estimates:
		    if(globalopts.smoothACEst)
		      { 
			acEst.spatiallySmooth(logger.getDir() + "/" + globalopts.epifname, epivol, globalopts.ms, globalopts.epifname, globalopts.susanpath, globalopts.epith);		    
		      }
		
		    acEst.pava();
		  }
	    
	      }
	    cerr << "Completed" << endl; 

	    // Loop through voxels:
	    cerr << "Prewhitening..." << endl;   
	    int co = 1;
	    
	    for(int i = 1; i <= numTS; i++)
	      {					
		ColumnVector I(sizeTS);
		I = 0;
		I(1) = 1;
	    
		ColumnVector xw(sizeTS);
		ColumnVector xprew(sizeTS);
	    
		acEst.setDesignMatrix(dms[i-1]);
	    
		// Use autocorr estimate to prewhiten data:
		xprew = x.getSeries(i);
	    
		Matrix designmattw;
	    
		// iters==1 is for high freq removel
		acEst.preWhiten(xprew, xw, i, designmattw, bool(iters==1));
	    
		if(co > 1000)
		  {
		    co = 1;
		    cerr << (float)i/(float)numTS << ",";
		  }
		else
		  co++;
	   
		x.setSeries(xw,i);
		dms[i-1] = designmattw;
	      }  
      	
	    cerr << "Completed" << endl;
	
	    // Add param number to "pe" to create filename:
	    ostringstream osc;
	    osc << iters - 1;
	
	    VolumeSeries& threshac = acEst.getEstimates();
	    int cutoff = sizeTS/2;
	    if(globalopts.tukey)
	      cutoff = globalopts.tukeysize;
	    threshac = threshac.Rows(1,cutoff);
	    VolumeInfo vinfo = x.getInfo();
	    vinfo.v = cutoff;
	    threshac.unthresholdSeries(vinfo,x.getPreThresholdPositions());
	    threshac.writeAsFloat(logger.getDir() + "/threshac" + osc.str().c_str());
	    threshac.thresholdSeries();

	    if(globalopts.verbose)
	      {	
		cerr << "Saving results... " << endl;
	    
		ColumnVector& countLargeE = acEst.getCountLargeE();
	
		logger.out(string("countLargeE") + strc, countLargeE);

		residuals.unthresholdSeries(x.getInfo(),x.getPreThresholdPositions());
		residuals.writeAsFloat(logger.getDir() + "/res4d" + osc.str().c_str());
		residuals.thresholdSeries();
		cerr << "Completed" << endl;
	      }
	  }
      }
    else // no estimation of autocorrelations
      {
	// do nothing
      }

    // Do once more to compute real param ests:
    cerr << "Computing parameter estimates..." << endl;
    for(int i = 1; i <= numTS; i++)
      {						    
	glimGls.setData(x.getSeries(i), dms[i-1], i);
	
	if(globalopts.verbose)
	  {
	    residuals.setSeries(glimGls.getResiduals(),i);
	  }
      }
    cerr << "Completed" << endl;

    // Write out necessary data:
    cerr << "Saving results... " << endl;

    residuals.unthresholdSeries(x.getInfo(),x.getPreThresholdPositions());
    residuals.writeAsFloat(logger.getDir() + "/res4d");

    // write out design matrices - a volume Series for each param
    if(globalopts.verbose_dms)
      {
	VolumeInfo vinfo = x.getInfo();
	for(int j = 1; j <= numParams; j++)
	  {	  
	    ostringstream osc;
	    osc << j;
	
	    VolumeSeries dmsmat(sizeTS, numTS);
	
	    for(int i = 1; i <= numTS; i++)
	      {
		dmsmat.setSeries(dms[i-1].Column(j),i);
	      }
	    dmsmat.setInfo(vinfo);
	    dmsmat.setPreThresholdPositions(x.getPreThresholdPositions());
	    dmsmat.unthresholdSeries();
	    dmsmat.writeAsFloat(logger.getDir() + "/dms" + osc.str().c_str());
	  }
	// output x
	x.unthresholdSeries();
	x.writeAsFloat(logger.getDir() + "/wx");
      }

    x.CleanUp();
    dms.clear();

    if(globalopts.verbose)
      {	
	VolumeInfo vinfo = x.getInfo();

       	// Save E
	VolumeSeries& E = acEst.getE();
	vinfo.v = acEst.getZeroPad();
	E.setInfo(vinfo);
	E.setPreThresholdPositions(x.getPreThresholdPositions());
	E.unthresholdSeries();	
	E.writeAsFloat(logger.getDir() + "/E");

      }

    glimGls.Save(x.getInfo(), x.getPreThresholdPositions());
    cerr << "Completed" << endl;
  }  
  catch(Exception p_excp) 
    {
      cerr << p_excp.what() << endl;
    }
  catch(...) 
    {
      cerr << "Image error" << endl;
    } 
  return 0;
}

