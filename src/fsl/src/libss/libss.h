
/*  libss.h - basic C include stuff

    Stephen Smith, FMRIB Image Analysis Group

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

#include "fslio/fslio.h"

#ifdef  __cplusplus
extern "C" {
#endif

#include "libssbase.h"

typedef struct {
  FDT *i;                                    /* image data */
  int x, y, z, t,                            /* image dimensions */
    dt, bpv;                                 /* datatype and bytes-per-voxel */
  double xv, yv, zv,                         /* voxel size (always all +ve) */
    xv0, yv0, zv0,                           /* voxel size as set in original header */
    tr,                                      /* seconds between volumes */
    min, max, thresh2, thresh98, thresh,     /* various brightness thresholds */
    lthresh, uthresh,                        /* lower and upper thresholds - some procs will ignore data outside these */
    dtmin, dtmax;                            /* min and max values for data type */
  char lut[24];                              /* lut file */
  float info;                                /* slice order information */
  short intent_code;                         /* information on what kind of image e.g. zstat */
  float intent_p1, intent_p2, intent_p3;     /* parameters like dof */
  FSLIO* miscinfo;                           /* miscellaneous header information */
} image_struct;  

FDT median(double, FDT*, int);
double dmedian(double, double*, int);
FDT mean(FDT*, int);
FDT TLI(image_struct, double, double, double);
int make_isometric (image_struct, image_struct*);
int p2c(int, int, int, double, double, double, double, double, double, double*, double*, double*);
FDT getp2c(image_struct, double, double, double, double, double, double);
int find_histogram (image_struct*, int*, int);
/*void find_roi_histogram (image_struct*, int, int, int, int, int, int, int*, int);*/
void find_thresholds (image_struct*, double);
void c_of_g (image_struct, double*, double*, double*);
void invert_y (image_struct);
double find_radius (image_struct, double);
void print_image_struct(image_struct);
void init_image_struct(image_struct*);
void tempfilt(image_struct, double, double, int);
void spatgauss(image_struct, FDT*, double, double, double);
void spatfilt(image_struct, double);
void spatmedian(image_struct im, int width);
void intnorm(image_struct, int, double);
void fmrimask(image_struct, image_struct);
void resample(image_struct, image_struct*, float);
void sample(image_struct, image_struct, float);

#ifdef  __cplusplus
}
#endif

