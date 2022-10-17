/* {{{ Copyright etc. */

/*  avwroi.c - extract cuboid ROI from image

    Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 1999-2005 University of Oxford  */

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

#include "libss/libss.h"
#include "libss/libavw.h"

void usage(void);

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("\nUsage: avwroi <input> <output> <xmin> <xsize> <ymin> <ysize> <zmin> <zsize>\n");
  printf("       avwroi <input> <output> <tmin> <tsize>\n");
  printf("       avwroi <input> <output> <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize>\n\n");
  printf("Note: counting (in both time and space) starts with 0 not 1!\n\n");
  exit(1);
}

/* }}} */
/* {{{ main */

int main(argc,argv) 
     int argc;
     char **argv;
{
  FDT *out;
  image_struct im;
  int x,y,z,t,xmin=0,xsize=0,ymin=0,ysize=0,zmin=0,zsize=0,tmin=0,tsize=0,argindex;

  /* {{{ process args */

if (argc==1)
     usage();

avw_read(argv[1],&im);
im.t=MAX(im.t,1);

argindex=3;

if (argc==5)
{
  xmin=0;
  xsize=im.x;
  ymin=0;
  ysize=im.y;
  zmin=0;
  zsize=im.z;
  tmin=atoi(argv[argindex++]);
  tsize=atoi(argv[argindex++]);
}
else if (argc==9)
{
  xmin=atoi(argv[argindex++]);
  xsize=atoi(argv[argindex++]);
  ymin=atoi(argv[argindex++]);
  ysize=atoi(argv[argindex++]);
  zmin=atoi(argv[argindex++]);
  zsize=atoi(argv[argindex++]);
  tmin=0;
  tsize=im.t;
}
else if (argc==11)
{
  xmin=atoi(argv[argindex++]);
  xsize=atoi(argv[argindex++]);
  ymin=atoi(argv[argindex++]);
  ysize=atoi(argv[argindex++]);
  zmin=atoi(argv[argindex++]);
  zsize=atoi(argv[argindex++]);
  tmin=atoi(argv[argindex++]);
  tsize=atoi(argv[argindex++]);
}
else
  usage();

out=malloc(xsize*ysize*zsize*tsize*sizeof(FDT));  

/* }}} */
  /* {{{ process roi and output */

  for(t=0;t<tsize;t++)
    for(z=0;z<zsize;z++)
      for(y=0;y<ysize;y++)
	for(x=0;x<xsize;x++)
	  out[t*zsize*ysize*xsize + z*ysize*xsize + y*xsize + x] =
	    im.i[(t+tmin)*im.z*im.y*im.x + (z+zmin)*im.y*im.x + (y+ymin)*im.x + x+xmin];

  im.i=out;
  im.x=xsize;	  
  im.y=ysize;	  
  im.z=zsize;	  
  im.t=tsize;	  

  avw_write(argv[2],im);

/* }}} */

  exit(0);
}

/* }}} */
