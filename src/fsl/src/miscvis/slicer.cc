/* {{{ Copyright etc. */

/*  slicer - produce any slicings of a 3D image

    Stephen Smith & Christian Beckmann, FMRIB Image Analysis Group

    Copyright (C) 1999-2002 University of Oxford  */

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
/* {{{ defines, includes and typedefs */
 
#include "libvis/miscpic.h"
//#include "fmribmain.h"
 
using namespace NEWIMAGE;
using namespace MISCPIC;

/* }}} */
/* {{{ usage */


void usage(void)
{
  printf("\nUsage: slicer <input> [input2] [main options] [output options - any number of these]\n\n");

  printf("Main options: [-l <lut>] [-s <scale>] [-u] [-i <intensitymin> <intensitymax>] [-t] [-n]\n");
  printf("These must be before output options and must be in the above order.\n");
  printf("-l <lut> : use a different colour map from that specified in the header.\n");
  printf("-t       : produce semi-transparent (dithered) edges.\n");
  printf("-n       : use nearest-neighbour interpolation for output.\n");
  printf("-u       : do not put left-right labels in output.\n\n");

  printf("Output options:\n");
  printf("[-x/y/z <slice> <filename>]      : output sagittal, coronal or axial slice\n     (if <slice> >0 it is a fraction of image dimension, if <0, it is an absolute slice number)\n");
  printf("[-a <filename>]                  : output mid-sagittal, -coronal and -axial slices into one image\n");
  printf("[-A <width> <filename>]          : output _all_ axial slices into one image of _max_ width <width>\n");
  printf("[-S <sample> <width> <filename>] : as -A but only include every <sample>'th slice\n\n");

  exit(1);
}

/* }}} */
/* {{{ main */

//template <class T>
int fmrib_main(int argc, char* argv[])
{
  volumeinfo vol1info;
  //volume<T> vol1, vol2(1,1,1);

  volume<float> vol1, vol2(1,1,1);

  read_volume(vol1,string(argv[1]),vol1info);
  vol1.makeradiological();

  int i = 2;
  if ( (argc>i) && (argv[i][0]!='-') ){
    read_volume(vol2,string(argv[i])); i++;
    vol2.makeradiological();
  }
  
  bool dbgflag = FALSE;
  if ( (argc>i) && (string(argv[i]) == string("-d")) ){
    dbgflag = TRUE; i++;
  }

  char tmp[10000];
  sprintf(tmp," ");
  for(;i<argc;i++){
    strcat(strcat(tmp,argv[i])," ");}
  
  //miscpic<T> newpic;
  miscpic newpic;

  return newpic.slicer(vol1, vol2, tmp, &vol1info, dbgflag);
}

int main(int argc,char *argv[])
{
  if (argc<2) usage();
  
  //  return call_fmrib_main(dtype(std::string(argv[1])),argc,argv);  
  //return call_fmrib_main(DT_FLOAT,argc,argv);  
  return fmrib_main(argc,argv); 

}

/* }}} */
