/*  ranopts.h

    Tim Behrens & Steve Smith (FMRIB) & Tom Nichols (UMich)

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

#if !defined(ranopts_h)
#define ranopts_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"


using namespace Utilities;

namespace RANDOMISE {

class ranopts {
 public:
  static ranopts& getInstance();
  ~ranopts() { delete gopt; }
  
  Option<bool> verbose;
  Option<bool> help;
  Option<bool> demean_data;
  Option<bool> one_samp;
  Option<string> in_fileroot;
  Option<string> maskname;
  Option<string> out_fileroot;
  Option<string> dm_file;
  Option<string> confound_file;
  Option<string> tc_file;
  //  Option<string> fc_file;
  Option<string> logdir;
  Option<bool> how_many_perms;
  Option<int> n_perm;
  //  Option<string> test_stat;
  Option<float> cluster_thresh;
  Option<float> clustermass_thresh;
  //  Option<float> p_thresh;
  Option<float> var_sm_sig;
  
  void parse_command_line(int argc, char** argv,Log& logger);
  
 private:
  ranopts();  
  const ranopts& operator=(ranopts&);
  ranopts(ranopts&);

  OptionParser options; 
      
  static ranopts* gopt;
  
};

 inline ranopts& ranopts::getInstance(){
   if(gopt == NULL)
     gopt = new ranopts();
   
   return *gopt;
 }

 inline ranopts::ranopts() :
   verbose(string("-V"), false, 
	   string("switch on diagnostic messages"), 
	   false, no_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   demean_data(string("-D"), false,
	string("demean data temporally before model fitting"),
	false, no_argument),
   one_samp(string("-1"), false,
	string("perform 1-sample group-mean test instead of generic permutation test"),
	false, no_argument),
   in_fileroot(string("-i"), "",
       string("~<in_root>\tinput file root"),
       true, requires_argument),
   maskname(string("-m"), "",
       string("~<mask>\tmask image"),
       false, requires_argument),
   out_fileroot(string("-o"), string(""),
	    string("~<out_root>\toutput file root"),
	    true, requires_argument),  
   dm_file(string("-d"), string(""),
	    string("~<design_matrix>\tdesign matrix file"),
	    false, requires_argument),  
   confound_file(string("-x"), string(""),
	    string("~<con_file>\tconfound matrix file"),
	    false, requires_argument),  
   tc_file(string("-t"), string(""),
	    string("~<t_contrasts>\tt contrasts file"),
	    false, requires_argument),  
     //   fc_file(string("-f"), string(""),
     //	    string("~<f_contrasts>\tf contrasts file"),
     //	    false, requires_argument),  
   logdir(string("-l"), string("logdir"),
	    string("~<logdir>\tlog directory"),
	    false, requires_argument),  
   how_many_perms(string("-q"), false,
	    string("print out how many unique permutations would be generated and exit"),
	    false, no_argument),  
   n_perm(string("-n"), 5000,
	    string("~<n_perm>\tnumber of permutations (default 5000, set to 0 for exhaustive)"),
	    false, requires_argument),  
   //   test_stat(string("-s"), string("cluster"),
   //    string("~<test_stat>\ttest statistic - options are max,voxelwise,cluster (default cluster)"),
   //     false, requires_argument),
   //   cluster_thresh(string("-c"), 1,
   //  string("~<c_thresh>\tcluster-forming threshold"),
   //  false, requires_argument),
   cluster_thresh(string("-c"), -1,
	  string("~<c_thresh>\tcarry out cluster-based thresholding (as well as max and voxelwise)"),
	  false, requires_argument),
   clustermass_thresh(string("-C"), -1,
	  string("~<cmass_thresh>\tcarry out cluster-mass-based thresholding (as well as max and voxelwise)"),
	  false, requires_argument),
   //p_thresh(string("-p"), 0.05,
   //    string("~<p_thresh>\tp threshold (default 0.05)"),
   //     false, requires_argument),
   var_sm_sig(string("-v"), 0,
	    string("~<std>\tuse variance smoothing (std is in mm)"),
	     false, requires_argument),
  
  options("randomise v0.11 - Tim Behrens & Steve Smith (FMRIB) & Tom Nichols (UMich)", "randomise -i data -o output -d design.mat")
{
     
    
     try {
       options.add(verbose);
       options.add(help);
       options.add(demean_data);
       options.add(one_samp);
       options.add(in_fileroot);
       options.add(maskname);
       options.add(out_fileroot);
       options.add(dm_file);
       options.add(confound_file);
       options.add(tc_file);
       //       options.add(fc_file);
       options.add(logdir);
       options.add(how_many_perms);
       options.add(n_perm);
       //       options.add(test_stat);       
       options.add(cluster_thresh);
       options.add(clustermass_thresh);
       //options.add(p_thresh);
       options.add(var_sm_sig);
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

