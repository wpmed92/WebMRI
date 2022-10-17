/*  probtrackxOptions.h

    Tim Behrens, FMRIB Image Analysis Group

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

#if !defined(probtrackxOptions_h)
#define probtrackxOptions_h

#include <string> 
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

//#include "newmatall.h"
using namespace Utilities;

namespace TRACT {

class probtrackxOptions {
 public:
  static probtrackxOptions& getInstance();
  ~probtrackxOptions() { delete gopt; }
  
  Option<int> verbose;
  Option<bool> help;
  Option<string> basename;
  Option<string> maskfile;
  Option<string> seedfile; 
  Option<string> mode;
  Option<string> targetfile;
  Option<bool> simpleout;
  Option<bool> s2tout;
  FmribOption<bool> matrix1out;
  FmribOption<bool> matrix2out;
  FmribOption<bool> maskmatrixout;
  Option<string> outfile;
  Option<string> rubbishfile;
  Option<string> seeds_to_dti;
  FmribOption<string> skipmask;
  Option<string> seedref;
  Option<string> mask2;
  Option<string> waypoints;
  Option<string> meshfile;
  FmribOption<string> lrmask;
  Option<string> logdir; 
  Option<bool> forcedir;
  Option<int> nparticles;
  Option<int> nsteps;
  Option<float> c_thr;
  FmribOption<float> fibthresh;
  Option<float> steplength;
  Option<bool> loopcheck;
  Option<bool> usef;
  Option<bool> randfib;
  Option<bool> modeuler;
  Option<int> rseed;
  void parse_command_line(int argc, char** argv,Log& logger);
  void modecheck();
  void modehelp();
  void matrixmodehelp();
  void status();
 private:
  probtrackxOptions();  
  const probtrackxOptions& operator=(probtrackxOptions&);
  probtrackxOptions(probtrackxOptions&);

  OptionParser options; 
      
  static probtrackxOptions* gopt;
  
};


 inline probtrackxOptions& probtrackxOptions::getInstance(){
   if(gopt == NULL)
     gopt = new probtrackxOptions();
   
   return *gopt;
 }

 inline probtrackxOptions::probtrackxOptions() :
  verbose(string("-V,--verbose"), 0, 
	  string("verbose level, [0-2]"), 
	  false, requires_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   basename(string("-s,--samples"), string("DTI"),
	       string("basename for samples files"),
	       true, requires_argument),  
   maskfile(string("-m,--mask"), string("mask"),
	    string("Bet binary mask file in diffusion space"),
	    true, requires_argument),
   seedfile(string("-x,--seed"), string("Seed"),
	    string("Seed volume or voxel depending on mode"),
	    true, requires_argument),
   mode(string("--mode"), string(""),
	string("use --mode=simple for single seed voxel"),
	    false, requires_argument),
  targetfile(string("--targetmasks"), string("cmasks"),
	    string("File containing a list of target masks - required for seeds_to_targets mode"),
	    false, requires_argument),
  simpleout(string("--opd"), false,
	    string("output path distribution"),
	    false, no_argument), 
  s2tout(string("--os2t"), false,
	 string("output seeds to targets"),
	 false, no_argument),
  matrix1out(string("--omatrix1"), false,
	  string("output matrix1"),
	  false, no_argument), 
  matrix2out(string("--omatrix2"), false,
	  string("output matrix2"),
	  false, no_argument), 
  maskmatrixout(string("--omaskmatrix"), false,
	    string("output maskmatrix"),
	  false, no_argument), 

  outfile(string("-o,--out"), string(""),
	    string("Output file (only for single seed voxel mode)"),
	    false, requires_argument),
   rubbishfile(string("--rubbish"), string(""),
	    string("Rubbish file"),
	    false, requires_argument),
   seeds_to_dti(string("--xfm"), string(""),
	      string("Transform Matrix taking seed space to DTI space default is to use the identity"),false, requires_argument),
  skipmask(string("--no_integrity"), string(""),
	   string("no explanation needed"),
	   false, requires_argument),
  seedref(string("--seedref"), string(""),
	 string("Reference vol to define seed space in simple mode - diffusion space assumed if absent"),
	 false, requires_argument),
  mask2(string("--mask2"), string(""),
	 string("second mask in twomask_symm mode."),
       false, requires_argument),
 waypoints(string("--waypoints"), string(""),
	 string("Waypoint mask or ascii list of waypoint masks."),
       false, requires_argument), 
 meshfile(string("--mesh"), string(""),
	 string(""),
       false, requires_argument),
  lrmask(string("--lrmask"), string(""),
	 string("low resolution binary brain mask for stroring connectivity distribution in matrix2 mode."),
       false, requires_argument),
  logdir(string("--dir"), string(""),
	    string("Directory to put the final volumes in - code makes this directory"),
	    false, requires_argument),
  forcedir(string("--forcedir"), false,
	 string("Use the actual directory name given - i.e. don't add + to make a new directory"),
	 false, no_argument),
  nparticles(string("-P,--nsamples"), 5000,
	 string("Number of samples"),
	 false, requires_argument),
   nsteps(string("-S,--nsteps"), 2000,
	    string("Number of steps per sample"),
	    false, requires_argument),
   c_thr(string("-c,--cthr"), 0.2, 
	 string("Curvature threshold"), 
	 false, requires_argument),
  fibthresh(string("--fibthresh"), 0.01, 
	    string("volume fraction before subsidary fibre orientations are considered (def 0.01)"), 
	 false, requires_argument),
   steplength(string("--steplength"), 0.5, 
	 string("steplength"), 
	 false, requires_argument),
   loopcheck(string("-l,--loopcheck"), false, 
	 string("perform loopchecks on paths - slower, but allows lower curvature threshold"), 
	 false, no_argument),
   usef(string("-f,--usef"), false, 
	 string("Use anisotropy to constrain tracking"), 
	 false, no_argument),
  randfib(string("--randfib"), false, 
	 string("Select randomly from the first two fibres when seeding"), 
	 false, no_argument),
  modeuler(string("--modeuler"), false, 
	   string("Use modified euler streamlining"), 
	   false, no_argument),
   rseed(string("--rseed"), 12345,
	 string("Random seed"),
	 false, requires_argument), 
   options("probtrack","probtrack -s <basename> -m <maskname> -x <seedfile> -o <output> --targetmasks=<textfile>\n probtrack --help\n")
   {
     
    
     try {
       options.add(verbose);
       options.add(help);
       options.add(basename);
       options.add(maskfile);
       options.add(seedfile); 
       options.add(mode);
       options.add(targetfile);
       options.add(skipmask);
       options.add(mask2);
       options.add(waypoints);
       options.add(meshfile);
       options.add(lrmask);
       options.add(seedref);
       options.add(logdir); 
       options.add(forcedir); 
       options.add(simpleout);
       options.add(s2tout);
       options.add(matrix1out);
       options.add(matrix2out);
       options.add(maskmatrixout);
       options.add(outfile);
       options.add(rubbishfile);
       options.add(seeds_to_dti);
       options.add(nparticles);
       options.add(nsteps);
       options.add(c_thr);
       options.add(fibthresh);
       options.add(steplength);
       options.add(loopcheck);
       options.add(usef);
       options.add(randfib);
       options.add(modeuler);
       options.add(rseed);
     }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
     
   }
}

#endif







