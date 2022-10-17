/* {{{ copyright */

/*  avwmerge.c - concatenate AVW files into a single output

    Steve Smith, David Flitney and Stuart Clare, FMRIB Image Analysis Group

    Copyright (C) 2000-2002 University of Oxford  */

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fslio/fslio.h"

#define MAXINPUTS 10000

/* }}} */
/* {{{ usage */

void usage()
{
  printf("\nUsage: avwmerge <-x/y/z/t> <output> <file1 file2 .......>\n\n");
  printf("-t : concatenate images in time\n");
  printf("-x : concatenate images in the x direction\n");
  printf("-y : concatenate images in the y direction\n");
  printf("-z : concatenate images in the z direction\n\n");
  printf("-a : auto-choose: single slices -> volume, volumes -> 4D (time series)\n\n");
  exit(1);
}

/* }}} */

int main(int argc, char **argv)
{
  /* {{{ vars */

  FSLIO *src=NULL, *dest;
  short x[MAXINPUTS], y[MAXINPUTS], z[MAXINPUTS], v[MAXINPUTS],
    X=0, Y=0, Z=0, V=0,
    xx, yy, zz, vv,
    t;
  char *buffer[MAXINPUTS], *outbuffer;
  unsigned int i, direction, bpp=0;
  char filename[10000];

/* }}} */

  /* {{{ process initial args */

if (argc<4)
     usage();

if (!strcmp(argv[1], "-t"))
     direction=0;
else if (!strcmp(argv[1], "-x"))
     direction=1;
else if (!strcmp(argv[1], "-y"))
     direction=2;
else if (!strcmp(argv[1], "-z"))
     direction=3;
else if (!strcmp(argv[1], "-a"))
     direction=4;
else usage();

/* }}} */
  /* {{{ loop over all inputs */

for(i = 0; i < argc-3; i++)
{
  strcpy(filename,argv[i+3]);
  
  if((src=FslOpen(FslMakeBaseName(filename),"r"))==NULL) {
    perror("avwmerge");
    exit(1);
  }
  
  FslGetDim(src,&x[i],&y[i],&z[i],&v[i]);
  
  if (i==0)
    {
      bpp = FslGetDataType(src, &t) / 8;
      V=v[0]; X=x[0]; Y=y[0]; Z=z[0];
      if(direction==4){
	if((z[0]<2)&&(v[0]<2))direction=3;
	else direction=0;
      }
    }
  else
    {
      if (direction==0) V+=v[i];
      if (direction==1) X+=x[i];
      if (direction==2) Y+=y[i];
      if (direction==3) Z+=z[i];
    }

  buffer[i] = malloc(x[i]*y[i]*z[i]*v[i]*bpp);
  FslReadVolumes(src, buffer[i], v[i]);
  
  if (i<argc-4)
    FslClose(src);
}

outbuffer = malloc(X * Y * Z * V * bpp);
xx=yy=zz=vv=0;

/*printf("%d %d %d %d\n",X,Y,Z,V);*/

for(i = 0; i < argc-3; i++)
{
  if (direction==0)
    {
      memcpy(outbuffer+X*Y*Z*vv*bpp,buffer[i],X*Y*Z*v[i]*bpp);
      vv+=v[i];
    }

  if (direction==1)
    {
      for(vv=0;vv<V;vv++)
	for(zz=0;zz<Z;zz++)
	  for(yy=0;yy<Y;yy++)
	    memcpy(outbuffer+(X*Y*Z*vv+X*Y*zz+X*yy+xx)*bpp,buffer[i]+(Y*Z*vv+Y*zz+yy)*x[i]*bpp,x[i]*bpp);
      xx+=x[i];
    }    

  if (direction==2)
    {
      for(vv=0;vv<V;vv++)
	for(zz=0;zz<Z;zz++)
	  memcpy(outbuffer+(Y*Z*vv+Y*zz+yy)*X*bpp,buffer[i]+(Z*vv+zz)*y[i]*X*bpp,X*y[i]*bpp);
      yy+=y[i];
    }    

  if (direction==3)
    {
      for(vv=0;vv<V;vv++)
	memcpy(outbuffer+(vv*Z+zz)*X*Y*bpp,buffer[i]+X*Y*z[i]*vv*bpp,X*Y*z[i]*bpp);
      zz+=z[i];
    }    
}

/* }}} */
  /* {{{ write output */

  strcpy(filename,argv[2]);

  dest = FslOpen(FslMakeBaseName(filename), "w");

  FslCloneHeader(dest, src);
  FslClose(src);

  FslSetDim(dest, X, Y, Z, V);
  FslSetDimensionality(dest, 4);

  FslWriteHeader(dest);
  FslWriteVolumes(dest, outbuffer, V);

  FslClose(dest);

/* }}} */

  return(0);
}
