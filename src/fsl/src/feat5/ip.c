/* {{{ Copyright etc. */

/*  ip - FMRI image preprocessing

    Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

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

/* }}} */
/* {{{ includes */

#include "libss/libss.h"
#include "libss/libavw.h"

/* }}} */
/* {{{ usage */

void usage(char *progname)
{
  printf("Usage: %s <input> <output> <%%threshold> [options]\n\n",progname);
  printf("-s <sigma> : spatial lowpass Gaussian linear filtering (sigma in mm, not voxels)\n");
  printf("-S <width> : spatial lowpass median filtering (width is full cube width in voxels)\n");
  printf("-m <mask_filename> : create binary mask, apply to fmri data and save mask to file\n");
  printf("-i <target mean> : global intensity normalisation (separate intensity rescaling for each 3D volume)\n");
  printf("-I <target mean> : grand mean intensity normalisation (one intensity scaling for the whole 4D data set)\n");
  printf("-t <hp_sigma> <lp_sigma> : Bandpass temporal filtering; nonlinear highpass and Gaussian linear lowpass (with sigmas in volumes, not seconds); set either sigma<0 to skip that filter\n");
  printf("-T <hp_sigma> <lp_sigma> : Same as -t, but only uses previous timepoints wrt each timepoint for filtering\n\n");
  exit(1);
}

/* }}} */
/* {{{ main(argc, argv) */

int main(argc, argv)
  int   argc;
  char  *argv [];
{
/* {{{ vars */

int          i;
double       percent_thresh;
image_struct im, mask;

/* }}} */

  /* {{{ process initial arguments and find threshold */

if (argc<4)
     usage(argv[0]);

avw_read(argv[1],&im);

percent_thresh = atof(argv[3]);

if ( (percent_thresh<0) || (percent_thresh>100) )
{
  printf("<%%threshold> must lie between 0 and 100\n");
  usage(argv[0]);
}

percent_thresh/=100.0;

find_thresholds (&im, percent_thresh);
/*printf("thresholds: 0=%f 2=%f t=%f 98=%f 100=%f\n",im.min,im.thresh2,im.thresh,im.thresh98,im.max);*/

/* }}} */
  /* {{{ process options */

for (i=4; i<argc; i++) {
  if (!strcmp(argv[i], "-s"))
    /* {{{ spatial Gaussian filtering */

{
  double sigma;

  i++;
  if (i==argc) usage(argv[0]);
  sigma = atof(argv[i]);

  if (sigma<0)
    {
      printf("Spatial sigma must be >=0\n");
      usage(argv[0]);
    }

  spatfilt(im,sigma);
}

/* }}} */
  else if (!strcmp(argv[i], "-S"))
    /* {{{ spatial median filtering */

{
  int width;

  i++;
  if (i==argc) usage(argv[0]);
  width = atoi(argv[i]);

  if (width<3)
    {
      printf("Width must be >=3\n");
      usage(argv[0]);
    }

  spatmedian(im,width);
}

/* }}} */
  else if (!strcmp(argv[i], "-m"))
    /* {{{ threshold and create mask */

{
  i++;
  if (i==argc) usage(argv[0]);

  mask=im;
  mask.t=1;
  mask.dt=DT_UNSIGNED_CHAR;
  mask.i=malloc(sizeof(unsigned char)*mask.x*mask.y*mask.z);
  mask.min=0;
  mask.max=1;

  fmrimask(im,mask);

  avw_write(argv[i],mask);
}

/* }}} */
  else if ( (!strcmp(argv[i], "-i")) || (!strcmp(argv[i], "-I")) )
    /* {{{ intensity normalisation */

{
  double target_mean;
  int grand_mean=0;

  if (!strcmp(argv[i], "-I"))
    grand_mean=1;

  i++;
  if (i==argc) usage(argv[0]);
  target_mean = atof(argv[i]);

  intnorm(im,grand_mean,target_mean);

  im.thresh=target_mean*percent_thresh;
}

/* }}} */
  else if ( (!strcmp(argv[i], "-t")) || (!strcmp(argv[i], "-T")) )
    /* {{{ temporal filtering */

{
  double hp_sigma, lp_sigma;
  int backwards=0;

  if (!strcmp(argv[i], "-T"))
    backwards=1;

  i++;
  if (i==argc) usage(argv[0]);
  hp_sigma = atof(argv[i]);

  i++;
  if (i==argc) usage(argv[0]);
  lp_sigma = atof(argv[i]);

  tempfilt(im,hp_sigma,lp_sigma,backwards);
}

/* }}} */
  else
    usage(argv[0]);
}

/* }}} */
  /* {{{ output and exit */

avw_write(argv[2],im);

printf("THRESHOLD %f\n",im.thresh);

return(0);

/* }}} */
}

/* }}} */
