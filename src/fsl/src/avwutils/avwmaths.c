/* {{{ Copyright etc. */

/*  avwmaths.c - image maths

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
#include "fslio/fslio.h"

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("\nUsage: avwmaths <first_input> [operations and inputs] <output>\n\n");

  printf("\nBinary operations:\n");
  printf("  some inputs can be either an image or a number\n");
  printf("  \"current image\" may be n_volumes>1 and the second n_volumes=1\n");
  printf("-add  : add following input to current image\n");
  printf("-sub  : subtract following input from current image\n");
  printf("-mul  : multiply current image by following input\n");
  printf("-div  : divide current image by following input\n");
  printf("-mas  : use (following image>0) to mask current image\n");
  printf("-thr  : use following number to threshold current image (zero anything below the number)\n");
  printf("-uthr : use following number to upper-threshold current image (zero anything above the number)\n");
  printf("-max  : take maximum of following input and current image\n");
  printf("-min  : take minimum of following input and current image\n");

  printf("\nUnary operations:\n");
  printf("-exp  : exponential\n");
  printf("-log  : natural logarithm\n");
  printf("-sqr  : square\n");
  printf("-sqrt : square root\n");
  printf("-abs  : absolute value\n");
  printf("-bin  : use (current image>0) to binarise\n");
  printf("-dil  : use (current image>0) to dilate using 3x3x3 neighourhood (new value is average of non-zero neighbours)\n");
  printf("-dils : use (current image>0) to dilate using 3x3x3 neighourhood (simple - new value is last non-zero neighbour\n");
  printf("-dil2 : use (current image>0) to dilate using 3x3 (2D) neighourhood (new value is average of non-zero neighbours)\n");
  printf("-ero  : use (current image>0) to erode using 3x3x3 neighourhood\n");
  printf("-ero2 : use (current image>0) to erode using 3x3 (2D) neighourhood\n");
  printf("-edge : edge strength\n");
  printf("-nms  : 3D non-maximum suppression wrt immediate 3x3x3 neighbourhood\n");
  printf("-nan  : replace NaNs (improper numbers) with 0\n");
  printf("-nanm : make NaN (improper number) mask with 1 for NaN voxels, 0 otherwise\n");
  printf("-roi <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> : zero outside roi\n");

  printf("\nDimensionality reduction operations:\n");
  printf("  The \"T\" can be replaced by X, Y or Z to collapse across a different dimension\n");
  printf("-Tmean    : mean across time\n");
  printf("-Tstd     : standard deviation across time\n");
  printf("-Tar1     : temporal AR(1) coefficient (use _32R and probably demean first)\n");
  printf("-Tmax     : max across time\n");
  printf("-Tmaxn    : time index of max across time\n");
  printf("-Tmin     : min across time\n");
  printf("-Tmedian  : median across time\n");
  printf("-Tperc <percentage> : nth percentile (0-100) across time\n");

  printf("\ne.g. avwmaths grotin -add grotin2 grotout\n");
  printf("     avwmaths grotin -add 2.5 grotout\n");
  printf("     avwmaths grotin -add 2.5 -mul grotin2 grotout\n\n");

  exit(1);
}

/* }}} */
/* {{{ isavw */

int isavw(const char *filename, double *tmpd)
{
  /* check if it is an image */
  if (FslFileExists(filename)) 
  return 1;

  /* if it is not an image */
  *tmpd=atof(filename);
  return 0;
}

/* }}} */
/* {{{ main */

int main(argc, argv)
  int   argc;
  char  *argv [];
{
  /* {{{ vars */

image_struct im, im2;
FDT  *in[1000];
int  i, j, x, y, z, t, xysize, xyzsize, size;
double tmpd;

/* }}} */

  if (argc<3)
    usage();

  avw_read(argv[1],&im);
  in[0]=im.i;
  xysize=im.x*im.y;
  xyzsize=xysize*im.z;
  size=xyzsize*im.t;

  for (i = 2; i < argc-1; i++) {
    if (isupper((int)argv[i][1]))
      /* {{{ group operations (dimensionality reduction) */

#define CI FDT tmpval, *tmpdata=malloc(sizeof(FDT)*MAX(MAX(im.x,im.y),MAX(im.z,im.t))); double sum=0,sumsquares=0,sumcorr=0,min=im.dtmax,max=im.dtmin; int nmax=0, n=0

#define CU {        if (argv[i][1] == 'T')  n=t;                                                                            \
               else if (argv[i][1] == 'Z')  n=z;                                                                            \
               else if (argv[i][1] == 'Y')  n=y;                                                                            \
               else if (argv[i][1] == 'X')  n=x;                                                                            \
                                                                                                                            \
             tmpval = in[0][t*im.z*im.y*im.x+z*im.y*im.x+y*im.x+x];                                                         \
             tmpdata[n] = tmpval;                                                                                           \
                                                                                                                            \
             sum += tmpval;                                                                                                 \
                                                                                                                            \
             sumsquares += tmpval*tmpval;                                                                                   \
	                                                                                                                    \
             if (n>0)                                                                                                       \
	       sumcorr += tmpval*tmpdata[n-1];                                                                              \
                                                                                                                            \
             if (tmpval>max) {                                                                                              \
               max=tmpval;                                                                                                  \
               nmax=n;}                                                                                                     \
                                                                                                                            \
             min = MIN(min,tmpval);}

#define CO n++;                                                                                                                \
           if (!strncmp(argv[i]+2, "mean", 4))   tmpim.i[t*tmpim.z*tmpim.y*tmpim.x+z*tmpim.y*tmpim.x+y*tmpim.x+x]=sum/n;       \
           else if (!strncmp(argv[i]+2, "std", 3))    tmpim.i[t*tmpim.z*tmpim.y*tmpim.x+z*tmpim.y*tmpim.x+y*tmpim.x+x]=        \
                                                                   sqrt((sumsquares/MAX(1,n-1))-((sum*sum)/(n*MAX(1,n-1))));   \
           else if (!strncmp(argv[i]+2, "ar1", 3))    tmpim.i[t*tmpim.z*tmpim.y*tmpim.x+z*tmpim.y*tmpim.x+y*tmpim.x+x]=        \
                                                                   sumcorr / sumsquares;                                       \
	   else if (!strncmp(argv[i]+2, "maxn", 4))   tmpim.i[t*tmpim.z*tmpim.y*tmpim.x+z*tmpim.y*tmpim.x+y*tmpim.x+x]=nmax;   \
	   else if (!strncmp(argv[i]+2, "max", 3))    tmpim.i[t*tmpim.z*tmpim.y*tmpim.x+z*tmpim.y*tmpim.x+y*tmpim.x+x]=max;    \
	   else if (!strncmp(argv[i]+2, "min", 3))    tmpim.i[t*tmpim.z*tmpim.y*tmpim.x+z*tmpim.y*tmpim.x+y*tmpim.x+x]=min;    \
	   else if ( (!strncmp(argv[i]+2, "median", 6)) || (!strncmp(argv[i]+2, "perc", 4) ) )                                 \
                                                      tmpim.i[t*tmpim.z*tmpim.y*tmpim.x+z*tmpim.y*tmpim.x+y*tmpim.x+x]=        \
                                                                   median(percfrac,tmpdata,n);                                 \
           free(tmpdata);

{
  image_struct tmpim=im;
  float percfrac=0.5;
  int doperc=0;
  tmpim.i=malloc(size*sizeof(FDT));

  if (!strncmp(argv[i]+2, "perc", 4))
    {
      percfrac=MAX(MIN(atof(argv[i+1]),100),0)/100;
      doperc=1;
    }

  if (argv[i][1] == 'T')
    {
      tmpim.t=1;
      for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
	{ CI; for(t=0;t<im.t;t++) CU; t=0; CO; }
    }
  else if (argv[i][1] == 'Z')
    {
      tmpim.z=1;
      for(t=0;t<im.t;t++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
	{ CI; for(z=0;z<im.z;z++) CU; z=0; CO; }
    }
  else if (argv[i][1] == 'Y')
    {
      tmpim.y=1;
      for(z=0;z<im.z;z++) for(t=0;t<im.t;t++) for(x=0;x<im.x;x++)
	{ CI; for(y=0;y<im.y;y++) CU; y=0; CO; }
    }
  else if (argv[i][1] == 'X')
    {
      tmpim.x=1;
      for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(t=0;t<im.t;t++)
	{ CI; for(x=0;x<im.x;x++) CU; x=0; CO; }
    }

  free(in[0]);
  im=tmpim;
  in[0]=im.i;
  xysize=im.x*im.y;
  xyzsize=xysize*im.z;
  size=xyzsize*im.t;

  if (doperc)
    i++;
}

/* }}} */
    else if (!strncmp(argv[i], "-mas", 4))
      /* {{{ mask */

{
  i++;
  avw_read(argv[i],&im2);
  in[1]=im2.i;
  for(j=0;j<size;j++)
    if (in[1][j%(im2.x*im2.y*im2.z*im2.t)]<=0)
      in[0][j]=0;
  free(in[1]);
}

/* }}} */
    else if (!strncmp(argv[i], "-add", 4))
      /* {{{ add */

{
  i++;
  if ( ! isavw(argv[i],&tmpd) )
    for(j=0;j<size;j++)
      in[0][j]+=tmpd;
  else
    {
      avw_read(argv[i],&im2);
      in[1]=im2.i;
      for(j=0;j<size;j++)
	in[0][j]+=in[1][j%(im2.x*im2.y*im2.z*im2.t)];
      free(in[1]);
    }
}

/* }}} */
    else if (!strncmp(argv[i], "-sub", 4))
      /* {{{ subtract */

{
  i++;
  if ( ! isavw(argv[i],&tmpd) )
    for(j=0;j<size;j++)
      in[0][j]-=tmpd;
  else
    {
      avw_read(argv[i],&im2);
      in[1]=im2.i;
      for(j=0;j<size;j++)
	  in[0][j]-=in[1][j%(im2.x*im2.y*im2.z*im2.t)];
      free(in[1]);
    }
}

/* }}} */
    else if (!strncmp(argv[i], "-mul", 4))
      /* {{{ multiply */

{
  i++;
  if ( ! isavw(argv[i],&tmpd) )
    for(j=0;j<size;j++)
      in[0][j]*=tmpd;
  else
    {
      avw_read(argv[i],&im2);
      in[1]=im2.i;
      for(j=0;j<size;j++)
	in[0][j]*=in[1][j%(im2.x*im2.y*im2.z*im2.t)];
      free(in[1]);
    }
}

/* }}} */
    else if (!strncmp(argv[i], "-div", 4))
      /* {{{ divide */

{
  i++;
  if ( ! isavw(argv[i],&tmpd) )
    {
      if (tmpd!=0)
	for(j=0;j<size;j++)
	  in[0][j]/=tmpd;
      else
	for(j=0;j<size;j++)
	  in[0][j]=0;
    }
  else
    {
      avw_read(argv[i],&im2);
      in[1]=im2.i;
      for(j=0;j<size;j++)
	if (in[1][j%(im2.x*im2.y*im2.z*im2.t)]!=0)
	  in[0][j]/=in[1][j%(im2.x*im2.y*im2.z*im2.t)];
	else
	  in[0][j]=0;
      free(in[1]);
    }
}

/* }}} */
    else if (!strncmp(argv[i], "-thr", 4))
      /* {{{ threshold */

{
  i++;
  if ( ! isavw(argv[i],&tmpd) )
    for(j=0;j<size;j++)
      if (in[0][j]<tmpd)
	in[0][j]=0;
}

/* }}} */
    else if (!strncmp(argv[i], "-uthr", 5))
      /* {{{ upper-threshold */

{
  i++;
  if ( ! isavw(argv[i],&tmpd) )
    for(j=0;j<size;j++)
      if (in[0][j]>tmpd)
	in[0][j]=0;
}

/* }}} */
    else if (!strncmp(argv[i], "-exp", 4))
      /* {{{ exp */

{
  for(j=0;j<size;j++)
    in[0][j]=exp(in[0][j]);
}

/* }}} */
    else if (!strncmp(argv[i], "-log", 4))
      /* {{{ log */

{
  for(j=0;j<size;j++)
    if (in[0][j]>0)
      in[0][j]=log(in[0][j]);
    else
      in[0][j]=0;
}

/* }}} */
    else if (!strncmp(argv[i], "-sqrt", 5))
      /* {{{ square root */

{
  for(j=0;j<size;j++)
    if (in[0][j]>0)
      in[0][j]=sqrt(in[0][j]);
    else
      in[0][j]=0;
}

/* }}} */
    else if (!strncmp(argv[i], "-sqr", 4))
      /* {{{ square */

{
  for(j=0;j<size;j++)
    in[0][j]*=in[0][j];
}

/* }}} */
    else if (!strncmp(argv[i], "-abs", 4))
      /* {{{ abs */

{
  for(j=0;j<size;j++)
    in[0][j]=ABS(in[0][j]);
}

/* }}} */
    else if (!strncmp(argv[i], "-bin", 4))
      /* {{{ binarise */

{
  for(j=0;j<size;j++)
    if (in[0][j]>0)
      in[0][j]=1;
    else
      in[0][j]=0;
}

/* }}} */
    else if (!strncmp(argv[i], "-max", 4))
      /* {{{ max */

{
  i++;
  if ( ! isavw(argv[i],&tmpd) )
    for(j=0;j<size;j++)
      in[0][j]=MAX(tmpd,in[0][j]);
  else
    {
      avw_read(argv[i],&im2);
      in[1]=im2.i;
      for(j=0;j<size;j++)
	in[0][j]=MAX(in[1][j%(im2.x*im2.y*im2.z*im2.t)],in[0][j]);
      free(in[1]);
    }
}

/* }}} */
    else if (!strncmp(argv[i], "-min", 4))
      /* {{{ min */

{
  i++;
  if ( ! isavw(argv[i],&tmpd) )
    for(j=0;j<size;j++)
      in[0][j]=MIN(tmpd,in[0][j]);
  else
    {
      avw_read(argv[i],&im2);
      in[1]=im2.i;
      for(j=0;j<size;j++)
	in[0][j]=MIN(in[1][j%(im2.x*im2.y*im2.z*im2.t)],in[0][j]);
      free(in[1]);
    }
}

/* }}} */
    else if (!strncmp(argv[i], "-dils", 5))
      /* {{{ dilate simple */

{
  FDT *tmpim=malloc(size*sizeof(FDT));
  memcpy(tmpim,in[0],size*sizeof(FDT));

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    if (in[0][t*xyzsize+z*xysize+y*im.x+x]==0)
      {
	int xx, yy, zz;
	for(zz=MAX(z-1,0);zz<=MIN(z+1,im.z-1);zz++)
	  for(yy=MAX(y-1,0);yy<=MIN(y+1,im.y-1);yy++)
	    for(xx=MAX(x-1,0);xx<=MIN(x+1,im.x-1);xx++)
	      if (in[0][t*xyzsize+zz*xysize+yy*im.x+xx]!=0)
		tmpim[t*xyzsize+z*xysize+y*im.x+x]=in[0][t*xyzsize+zz*xysize+yy*im.x+xx];
      }

  free(in[0]);
  in[0]=tmpim;
}

/* }}} */
    else if (!strncmp(argv[i], "-dil2", 5))
      /* {{{ dilate 2D */

{
  FDT *tmpim=malloc(size*sizeof(FDT));
  memcpy(tmpim,in[0],size*sizeof(FDT));

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    if (in[0][t*xyzsize+z*xysize+y*im.x+x]==0)
      {
	int xx, yy, count=0;
	double total=0;
	for(yy=MAX(y-1,0);yy<=MIN(y+1,im.y-1);yy++)
	  for(xx=MAX(x-1,0);xx<=MIN(x+1,im.x-1);xx++)
	    if (in[0][t*xyzsize+z*xysize+yy*im.x+xx]!=0)
	      {
		total+=in[0][t*xyzsize+z*xysize+yy*im.x+xx];
		count++;
	      }
	if (count>0)
	  tmpim[t*xyzsize+z*xysize+y*im.x+x]=total/count;
      }

  free(in[0]);
  in[0]=tmpim;
}

/* }}} */
    else if (!strncmp(argv[i], "-dil", 4))
      /* {{{ dilate */

{
  FDT *tmpim=malloc(size*sizeof(FDT));
  memcpy(tmpim,in[0],size*sizeof(FDT));

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    if (in[0][t*xyzsize+z*xysize+y*im.x+x]==0)
      {
	int xx, yy, zz, count=0;
	double total=0;
	for(zz=MAX(z-1,0);zz<=MIN(z+1,im.z-1);zz++)
	  for(yy=MAX(y-1,0);yy<=MIN(y+1,im.y-1);yy++)
	    for(xx=MAX(x-1,0);xx<=MIN(x+1,im.x-1);xx++)
	      if (in[0][t*xyzsize+zz*xysize+yy*im.x+xx]!=0)
		{
		  total+=in[0][t*xyzsize+zz*xysize+yy*im.x+xx];
		  count++;
		}
	if (count>0)
	  tmpim[t*xyzsize+z*xysize+y*im.x+x]=total/count;
      }

  free(in[0]);
  in[0]=tmpim;
}

/* }}} */
    else if (!strncmp(argv[i], "-ero2", 5))
      /* {{{ erode */

{
  FDT *tmpim=malloc(size*sizeof(FDT));
  memcpy(tmpim,in[0],size*sizeof(FDT));

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    if (in[0][t*xyzsize+z*xysize+y*im.x+x]!=0)
      {
	int xx, yy;
	for(yy=MAX(y-1,0);yy<=MIN(y+1,im.y-1);yy++)
	  for(xx=MAX(x-1,0);xx<=MIN(x+1,im.x-1);xx++)
	    if (in[0][t*xyzsize+z*xysize+yy*im.x+xx]==0)
	      tmpim[t*xyzsize+z*xysize+y*im.x+x]=0;
      }

  free(in[0]);
  in[0]=tmpim;
}

/* }}} */
    else if (!strncmp(argv[i], "-ero", 4))
      /* {{{ erode */

{
  FDT *tmpim=malloc(size*sizeof(FDT));
  memcpy(tmpim,in[0],size*sizeof(FDT));

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    if (in[0][t*xyzsize+z*xysize+y*im.x+x]!=0)
      {
	int xx, yy, zz;
	for(zz=MAX(z-1,0);zz<=MIN(z+1,im.z-1);zz++)
	  for(yy=MAX(y-1,0);yy<=MIN(y+1,im.y-1);yy++)
	    for(xx=MAX(x-1,0);xx<=MIN(x+1,im.x-1);xx++)
	      if (in[0][t*xyzsize+zz*xysize+yy*im.x+xx]==0)
		tmpim[t*xyzsize+z*xysize+y*im.x+x]=0;
      }

  free(in[0]);
  in[0]=tmpim;
}

/* }}} */
    else if (!strncmp(argv[i], "-edge", 5))
      /* {{{ edge */

{
  float tmpf=2 * sqrt(1/(im.xv*im.xv) + 1/(im.yv*im.yv) + 1/(im.zv*im.zv));
  FDT *tmpim=malloc(size*sizeof(FDT));
  memcpy(tmpim,in[0],size*sizeof(FDT));

  if (im.z>2)
    {

  for(t=0;t<im.t;t++) for(z=1;z<im.z-1;z++) for(y=1;y<im.y-1;y++) for(x=1;x<im.x-1;x++)
    tmpim[t*xyzsize+z*xysize+y*im.x+x] = sqrt (

	( ((double)in[0][t*xyzsize+(z+1)*xysize+y*im.x+x]-(double)in[0][t*xyzsize+(z-1)*xysize+y*im.x+x]) *
	  ((double)in[0][t*xyzsize+(z+1)*xysize+y*im.x+x]-(double)in[0][t*xyzsize+(z-1)*xysize+y*im.x+x]) ) 
	/ (im.zv*im.zv)

	+

	( ((double)in[0][t*xyzsize+z*xysize+(y+1)*im.x+x]-(double)in[0][t*xyzsize+z*xysize+(y-1)*im.x+x]) *
	  ((double)in[0][t*xyzsize+z*xysize+(y+1)*im.x+x]-(double)in[0][t*xyzsize+z*xysize+(y-1)*im.x+x]) )
	/ (im.yv*im.yv)

	+

	( ((double)in[0][t*xyzsize+z*xysize+y*im.x+(x+1)]-(double)in[0][t*xyzsize+z*xysize+y*im.x+(x-1)]) *
	  ((double)in[0][t*xyzsize+z*xysize+y*im.x+(x+1)]-(double)in[0][t*xyzsize+z*xysize+y*im.x+(x-1)]) )
	/ (im.xv*im.xv)

	) / tmpf;

    } else {

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=1;y<im.y-1;y++) for(x=1;x<im.x-1;x++)
    tmpim[t*xyzsize+z*xysize+y*im.x+x] = sqrt (

	( ((double)in[0][t*xyzsize+z*xysize+(y+1)*im.x+x]-(double)in[0][t*xyzsize+z*xysize+(y-1)*im.x+x]) *
	  ((double)in[0][t*xyzsize+z*xysize+(y+1)*im.x+x]-(double)in[0][t*xyzsize+z*xysize+(y-1)*im.x+x]) )
	/ (im.yv*im.yv)

	+

	( ((double)in[0][t*xyzsize+z*xysize+y*im.x+(x+1)]-(double)in[0][t*xyzsize+z*xysize+y*im.x+(x-1)]) *
	  ((double)in[0][t*xyzsize+z*xysize+y*im.x+(x+1)]-(double)in[0][t*xyzsize+z*xysize+y*im.x+(x-1)]) )
	/ (im.xv*im.xv)

	) / tmpf;

    }

  free(in[0]);
  in[0]=tmpim;
}

/* }}} */
    else if (!strncmp(argv[i], "-nms", 4))
      /* {{{ non-max-suppression */

{
  int maxsearch=0;
  float origthresh=0;
  FDT *tmpdata=NULL;
  image_struct distancemap, thedata, lowercingulum;
	
  if (!strncmp(argv[i], "-nmss", 5))
    {
      maxsearch=1;
      i++;
      origthresh=atof(argv[i]);
      i++;
      avw_read(argv[i],&distancemap);
      i++;
      avw_read(argv[i],&lowercingulum);
      i++;
      avw_read(argv[i],&thedata);
      
      tmpdata=calloc(thedata.t*size,sizeof(FDT));
    }

  /* {{{ nms 3x3x3 CofG etc with flow smoothing */

{
  int MAXZ=1;
  FDT *tmpim=calloc(size,sizeof(FDT));
  short *X=calloc(size,sizeof(short)),
    *Y=calloc(size,sizeof(short)),
    *Z=calloc(size,sizeof(short)),
    *XX=calloc(size,sizeof(short)),
    *YY=calloc(size,sizeof(short)),
    *ZZ=calloc(size,sizeof(short));

  if (im.z<3)
    MAXZ=0;

  /* {{{ estimate perp from CofG and curvature, and store this in X,Y,Z */

      for(z=MAXZ;z<im.z-MAXZ;z++) for(y=1;y<im.y-1;y++) for(x=1;x<im.x-1;x++)
	{
	  FDT theval = in[0][z*xysize+y*im.x+x];

	  if ( theval != 0 )
	    {
	      float CofGx=0, CofGy=0, CofGz=0, Sum=0, CofGl;
	      int xx, yy, zz, xxx=0, yyy=0, zzz=0;

	      for(zz=-MAXZ; zz<=MAXZ; zz++)
		for(yy=-1; yy<=1; yy++)
		  for(xx=-1; xx<=1; xx++)
		    {
		      float val = in[0][(z+zz)*xysize+(y+yy)*im.x+x+xx];
		      Sum   += val;
		      CofGx += xx * val;  CofGy += yy * val;  CofGz += zz * val;
		    }	      
	      
	      CofGx /= Sum;  CofGy /= Sum;  CofGz /= Sum;
	      CofGl = sqrt(CofGx*CofGx+CofGy*CofGy+CofGz*CofGz);
	      
	      if (CofGl > .1)  /* is CofG far enough away from centre voxel? */
		{
		  xxx = FTOI(CofGx/CofGl);
		  yyy = FTOI(CofGy/CofGl);
		  zzz = FTOI(CofGz/CofGl);
		}
	      else
		/* {{{ find direction of max curvature */

	{
	  float maxcost=0, centreval=2.0 * (float)theval;

	  for(zz=0; zz<=MAXZ; zz++) /* note - starts at zero as we're only searching half the voxels */
	    for(yy=-1; yy<=1; yy++)
	      for(xx=-1; xx<=1; xx++)
		if ( (zz==1) || (yy==1) || ((yy==0)&&(xx==1)) ) /* only search half the voxels */
		  {
		    float weighting = pow( (float)(xx*xx+yy*yy+zz*zz) , -0.7 ); /* power is arbitrary: maybe test other functions here */
		    float cost = weighting * ( centreval 
					       - (float)in[0][(z+zz)*xysize+(y+yy)*im.x+x+xx]
					       - (float)in[0][(z-zz)*xysize+(y-yy)*im.x+x-xx] );

		    if (cost>maxcost)
		      {
			maxcost=cost;
			xxx=xx;
			yyy=yy;
			zzz=zz;
		      }
		  }
	}

/* }}} */
							   
	      X[z*xysize+y*im.x+x]=xxx;
	      Y[z*xysize+y*im.x+x]=yyy;
	      Z[z*xysize+y*im.x+x]=zzz;
	    }

	}

      /* {{{ COMMENT save perp image */

#ifdef FoldingComment

{
  char tmpname[10000];
  image_struct tmpim=im;

  tmpim.i=calloc(size*3,sizeof(float));
  sprintf(tmpname,"%s_flow",argv[argc-1]);

  tmpim.dt=DT_FLOAT;
  tmpim.t=3;

  for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    {
      float tmpX=X[z*xysize+y*im.x+x],
	tmpY=Y[z*xysize+y*im.x+x],
	tmpZ=Z[z*xysize+y*im.x+x],
	tmpf=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);

      if (tmpf>0)
	{
	  tmpim.i[          z*xysize+y*im.x+x]=tmpX/tmpf;
	  tmpim.i[  xyzsize+z*xysize+y*im.x+x]=tmpY/tmpf;
	  tmpim.i[2*xyzsize+z*xysize+y*im.x+x]=tmpZ/tmpf;
	}
    }

  avw_write(tmpname,tmpim);
  free(tmpim.i);
}

#endif

/* }}} */

/* }}} */
  /* {{{ smooth X,Y,Z and store in XX,YY,ZZ */

      for(z=MAXZ;z<im.z-MAXZ;z++) for(y=1;y<im.y-1;y++) for(x=1;x<im.x-1;x++)
	{
	  int xx, yy, zz, *localsum=calloc(27,sizeof(int)), localmax=0, xxx, yyy, zzz;

	  for(zz=-MAXZ; zz<=MAXZ; zz++)
	    for(yy=-1; yy<=1; yy++)
	      for(xx=-1; xx<=1; xx++)
		{
		  xxx = X[(z+zz)*xysize+(y+yy)*im.x+x+xx];
		  yyy = Y[(z+zz)*xysize+(y+yy)*im.x+x+xx];
		  zzz = Z[(z+zz)*xysize+(y+yy)*im.x+x+xx];
		  localsum[(1+zzz)*9+(1+yyy)*3+1+xxx]++;
		  localsum[(1-zzz)*9+(1-yyy)*3+1-xxx]++;
		}

	  for(zz=-MAXZ; zz<=MAXZ; zz++)
	    for(yy=-1; yy<=1; yy++)
	      for(xx=-1; xx<=1; xx++)
		{
		  if (localsum[(1+zz)*9+(1+yy)*3+1+xx]>localmax)
		    {
		      localmax=localsum[(1+zz)*9+(1+yy)*3+1+xx];
		      XX[z*xysize+y*im.x+x]=xx;
		      YY[z*xysize+y*im.x+x]=yy;
		      ZZ[z*xysize+y*im.x+x]=zz;
		    }
		}
	}

      /* {{{ COMMENT save perp image */

#ifdef FoldingComment

{
  char tmpname[10000];
  image_struct tmpim=im;

  tmpim.i=calloc(size*3,sizeof(float));
  sprintf(tmpname,"%s_flowsmooth",argv[argc-1]);

  tmpim.dt=DT_FLOAT;
  tmpim.t=3;

  for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    {
      float tmpX=XX[z*xysize+y*im.x+x],
	tmpY=YY[z*xysize+y*im.x+x],
	tmpZ=ZZ[z*xysize+y*im.x+x],
	tmpf=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);

      if (tmpf>0)
	{
	  tmpim.i[          z*xysize+y*im.x+x]=tmpX/tmpf;
	  tmpim.i[  xyzsize+z*xysize+y*im.x+x]=tmpY/tmpf;
	  tmpim.i[2*xyzsize+z*xysize+y*im.x+x]=tmpZ/tmpf;
	}
    }

  avw_write(tmpname,tmpim);
  free(tmpim.i);
}

#endif

/* }}} */

      free(X); free(Y); free(Z);

/* }}} */
  /* {{{ do non-max-suppression in the direction of perp */

      for(z=MAXZ;z<im.z-MAXZ;z++) for(y=1;y<im.y-1;y++) for(x=1;x<im.x-1;x++)
	{
	  FDT theval = in[0][z*xysize+y*im.x+x];
	  int xxx=XX[z*xysize+y*im.x+x];
	  int yyy=YY[z*xysize+y*im.x+x];
	  int zzz=ZZ[z*xysize+y*im.x+x];
	  
	  if ( ( (xxx!=0) || (yyy!=0) || (zzz!=0) ) &&
	       ( theval >= in[0][(z+zzz)*xysize+(y+yyy)*im.x+x+xxx] ) &&
	       ( theval >  in[0][(z-zzz)*xysize+(y-yyy)*im.x+x-xxx] ) &&
	       ( theval >= in[0][(z+2*zzz)*xysize+(y+2*yyy)*im.x+x+2*xxx] ) &&
	       ( theval >  in[0][(z-2*zzz)*xysize+(y-2*yyy)*im.x+x-2*xxx] ) )
	    tmpim[z*xysize+y*im.x+x] = theval;
	}

/* }}} */

  if (maxsearch==1)
    /* {{{ do search in the direction of perp */

#define SEARCHSIGMA 10 /* length in linear voxel dimensions */
#define MAXSEARCHLENGTH (3*SEARCHSIGMA)

/*#define FLOWDEBUG*/

{
  int d, T, iters;

#ifdef FLOWDEBUG
  signed short *tmpimFLOWx=calloc(size*thedata.t,sizeof(short)),
    *tmpimFLOWy=calloc(size*thedata.t,sizeof(short)),
    *tmpimFLOWz=calloc(size*thedata.t,sizeof(short));
#endif  

  for(T=0;T<thedata.t;T++) for(z=MAXZ;z<im.z-MAXZ;z++) for(y=1;y<im.y-1;y++) for(x=1;x<im.x-1;x++)
    if (tmpim[z*xysize+y*im.x+x] > origthresh)
      {
	int xxx=XX[z*xysize+y*im.x+x];
	int yyy=YY[z*xysize+y*im.x+x];
	int zzz=ZZ[z*xysize+y*im.x+x];
	float maxval=thedata.i[T*xyzsize+z*xysize+y*im.x+x];
	float maxval_weighted=maxval;
	short maxvalX=0, maxvalY=0, maxvalZ=0;
	float exponentfactor = -0.5 * (xxx*xxx+yyy*yyy+zzz*zzz) / (float)(SEARCHSIGMA*SEARCHSIGMA);

	/* need to put proper FOV bounds checks in both below searches !!! */

	if (lowercingulum.i[z*xysize+y*im.x+x] == 0)
	  /* {{{ search perp to sheet */

{
  for(iters=0;iters<2;iters++)
    {
      float distance=0;

      for(d=1;d<MAXSEARCHLENGTH;d++)
	{
	  int D=d;
	  if (iters==1) D=-d;
		    
	  if (distancemap.i[(z+zzz*D)*xysize+(y+yyy*D)*im.x+x+xxx*D]>=distance)
	    {
	      float distanceweight = exp(d * d * exponentfactor);
	      distance=distancemap.i[(z+zzz*D)*xysize+(y+yyy*D)*im.x+x+xxx*D];
	      if (distanceweight * thedata.i[T*xyzsize+(z+zzz*D)*xysize+(y+yyy*D)*im.x+x+xxx*D]>maxval_weighted)
		{
		  maxval=thedata.i[T*xyzsize+(z+zzz*D)*xysize+(y+yyy*D)*im.x+x+xxx*D];
		  maxval_weighted=maxval*distanceweight;
		  maxvalX=xxx*D;
		  maxvalY=yyy*D;
		  maxvalZ=zzz*D;
		}
	    }
	  else
	    d=MAXSEARCHLENGTH;
	}
    }
}

/* }}} */
	else
	  /* {{{ search all around tube */

{
  for(yyy=-MAXSEARCHLENGTH; yyy<=MAXSEARCHLENGTH; yyy++) for(xxx=-MAXSEARCHLENGTH; xxx<=MAXSEARCHLENGTH; xxx++) 
    {
      float distanceweight = exp(-0.5 * (xxx*xxx+yyy*yyy) / (float)(SEARCHSIGMA*SEARCHSIGMA) );

      float r=sqrt((float)(xxx*xxx+yyy*yyy));

      if (r>0)
	{
	  float rr;
	  int allok=1;
	  
	  for(rr=1; rr<=r+0.1; rr++) /* search outwards from centre to current voxel - test that distancemap always increasing */
	    {
	      int xxx1=FTOI(rr*xxx/r);
	      int yyy1=FTOI(rr*yyy/r);
	      int xxx2=FTOI((rr+1)*xxx/r);
	      int yyy2=FTOI((rr+1)*yyy/r);
	      if ( distancemap.i[z*xysize+(y+yyy1)*im.x+x+xxx1] > distancemap.i[z*xysize+(y+yyy2)*im.x+x+xxx2] )
		allok=0;
	    }

	  if ( allok &&
	       ( distanceweight * thedata.i[T*xyzsize+z*xysize+(y+yyy)*im.x+x+xxx] > maxval_weighted ) )
	    {
	      maxval=thedata.i[T*xyzsize+z*xysize+(y+yyy)*im.x+x+xxx];
	      maxval_weighted=maxval*distanceweight;
	      maxvalX=xxx;
	      maxvalY=yyy;
	      maxvalZ=0;
	    }
	}
    }

}

/* }}} */
	
	tmpdata[T*xyzsize+z*xysize+y*im.x+x]=maxval; /* output maxsearch data */

#ifdef FLOWDEBUG
	tmpimFLOWx[T*xyzsize+z*xysize+y*im.x+x]=maxvalX;
	tmpimFLOWy[T*xyzsize+z*xysize+y*im.x+x]=maxvalY;
	tmpimFLOWz[T*xyzsize+z*xysize+y*im.x+x]=maxvalZ;
#endif

      }

  thedata.i=tmpdata;
  avw_write(argv[argc-1],thedata);

#ifdef FLOWDEBUG
  {
    char tmpname[10000];
    image_struct grotim=im;
    grotim.dt=DT_SIGNED_SHORT;
    grotim.t=thedata.t;
    
    grotim.i=tmpimFLOWx;
    sprintf(tmpname,"%s_search_X",argv[argc-1]);
    avw_write(tmpname,grotim);
    grotim.i=tmpimFLOWy;
    sprintf(tmpname,"%s_search_Y",argv[argc-1]);
    avw_write(tmpname,grotim);
    grotim.i=tmpimFLOWz;
    sprintf(tmpname,"%s_search_Z",argv[argc-1]);
    avw_write(tmpname,grotim);
  }
#endif
  
  exit(0);
}

/* }}} */

  free(XX); free(YY); free(ZZ);

  free(in[0]);
  in[0]=tmpim;
}

/* }}} */
  /* {{{ COMMENT NEW2 nms 3x3x3 CofG etc with flow smoothing */

#ifdef FoldingComment

/* make thinner skeleton using NMS in only 3 main directions
   then smooth perp direction a lot
   then count how many neighbours, look at local COG, etc. */

{
  int MAXZ=1;
  FDT *tmpim=calloc(size,sizeof(FDT));

  if (im.z<3)
    MAXZ=0;

  for(t=0;t<im.t;t++) /* kill this - shouldn't do this over t */
    {
      /* {{{ setup temp images to store perp */

      short *X=calloc(size,sizeof(short)),
	*Y=calloc(size,sizeof(short)),
	*Z=calloc(size,sizeof(short)),
	*XX=calloc(size,sizeof(short)),
	*YY=calloc(size,sizeof(short)),
	*ZZ=calloc(size,sizeof(short));

/* }}} */

      /* {{{ estimate perp from CofG and curvature, and store this in X,Y,Z */

      for(z=MAXZ;z<im.z-MAXZ;z++) for(y=1;y<im.y-1;y++) for(x=1;x<im.x-1;x++)
	{
	  FDT theval = in[0][t*xyzsize+z*xysize+y*im.x+x];

	  if ( theval != 0 )
	    {
	      float CofGx=0, CofGy=0, CofGz=0, Sum=0, CofGl;
	      int xx, yy, zz, xxx=0, yyy=0, zzz=0;

	      for(zz=-MAXZ; zz<=MAXZ; zz++)
		for(yy=-1; yy<=1; yy++)
		  for(xx=-1; xx<=1; xx++)
		    {
		      float val = in[0][t*xyzsize+(z+zz)*xysize+(y+yy)*im.x+x+xx];
		      Sum   += val;
		      CofGx += xx * val;  CofGy += yy * val;  CofGz += zz * val;
		    }	      
	      
	      CofGx /= Sum;  CofGy /= Sum;  CofGz /= Sum;
	      CofGl = sqrt(CofGx*CofGx+CofGy*CofGy+CofGz*CofGz);
	      
	      if (CofGl > .1)  /* is CofG far enough away from centre voxel? */
		{
		  /*xxx = FTOI(CofGx/CofGl);
		  yyy = FTOI(CofGy/CofGl);
		  zzz = FTOI(CofGz/CofGl);*/
		  if ( (fabs(CofGx)>fabs(CofGy)) && (fabs(CofGx)>fabs(CofGz)) )
		    { xxx=1; yyy=0; zzz=0; }
		  else
		    {
		      if ( fabs(CofGy)>fabs(CofGz) )
			{ xxx=0; yyy=1; zzz=0; }
		      else
			{ xxx=0; yyy=0; zzz=1; }
		    }
		}
	      else
		/* {{{ find direction of max curvature */

	{
	  float maxcost=0, centreval=2.0 * (float)theval;

	  for(zz=0; zz<=MAXZ; zz++) /* note - starts at zero as we're only searching half the voxels */
	    for(yy=0; yy<=1; yy++)
	      for(xx=0; xx<=1; xx++)
		if ( xx+yy+zz == 1 )
		  {
		    float weighting = pow( (float)(xx*xx+yy*yy+zz*zz) , -0.7 ); /* power is arbitrary: maybe test other functions here */
		    float cost = weighting * ( centreval 
					       - (float)in[0][t*xyzsize+(z+zz)*xysize+(y+yy)*im.x+x+xx]
					       - (float)in[0][t*xyzsize+(z-zz)*xysize+(y-yy)*im.x+x-xx] );

		    if (cost>maxcost)
		      {
			maxcost=cost;
			xxx=xx;
			yyy=yy;
			zzz=zz;
		      }
		  }
	}

/* }}} */
							   
	      X[z*xysize+y*im.x+x]=xxx;
	      Y[z*xysize+y*im.x+x]=yyy;
	      Z[z*xysize+y*im.x+x]=zzz;
	    }

	}

      /* {{{ COMMENT save perp image */

#ifdef FoldingComment

{
  char tmpname[10000];
  image_struct tmpim=im;

  tmpim.i=calloc(size*3,sizeof(float));
  sprintf(tmpname,"%s_flow",argv[argc-1]);

  tmpim.dt=DT_FLOAT;
  tmpim.t=3;

  for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    {
      float tmpX=X[z*xysize+y*im.x+x],
	tmpY=Y[z*xysize+y*im.x+x],
	tmpZ=Z[z*xysize+y*im.x+x],
	tmpf=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);

      if (tmpf>0)
	{
	  tmpim.i[          z*xysize+y*im.x+x]=tmpX/tmpf;
	  tmpim.i[  xyzsize+z*xysize+y*im.x+x]=tmpY/tmpf;
	  tmpim.i[2*xyzsize+z*xysize+y*im.x+x]=tmpZ/tmpf;
	}
    }

  avw_write(tmpname,tmpim);
  free(tmpim.i);
}

#endif

/* }}} */

/* }}} */
      /* {{{ smooth X,Y,Z and store in XX,YY,ZZ */

#define SMOOTHSEARCH 3

      for(z=MAXZ+SMOOTHSEARCH;z<im.z-MAXZ-SMOOTHSEARCH;z++) for(y=SMOOTHSEARCH;y<im.y-SMOOTHSEARCH;y++) for(x=SMOOTHSEARCH;x<im.x-SMOOTHSEARCH;x++)
	{
	  int xx, yy, zz, xxx, yyy, zzz;
	  float *localsum=calloc(27,sizeof(float)), localmax=0;

	  for(zz=-SMOOTHSEARCH; zz<=SMOOTHSEARCH; zz++)
	    for(yy=-SMOOTHSEARCH; yy<=SMOOTHSEARCH; yy++)
	      for(xx=-SMOOTHSEARCH; xx<=SMOOTHSEARCH; xx++)
		{
		  xxx = X[(z+zz)*xysize+(y+yy)*im.x+x+xx];
		  yyy = Y[(z+zz)*xysize+(y+yy)*im.x+x+xx];
		  zzz = Z[(z+zz)*xysize+(y+yy)*im.x+x+xx];
		  localsum[(1+zzz)*9+(1+yyy)*3+1+xxx] += in[0][t*xyzsize+(z+zz)*xysize+(y+yy)*im.x+x+xx];
		  localsum[(1-zzz)*9+(1-yyy)*3+1-xxx] += in[0][t*xyzsize+(z+zz)*xysize+(y+yy)*im.x+x+xx];
		}

  	  for(zz=-1; zz<=1; zz++)
	    for(yy=-1; yy<=1; yy++)
	      for(xx=-1; xx<=1; xx++)
		{
		  if (localsum[(1+zz)*9+(1+yy)*3+1+xx]>localmax)
		    {
		      localmax=localsum[(1+zz)*9+(1+yy)*3+1+xx];
		      XX[z*xysize+y*im.x+x]=xx;
		      YY[z*xysize+y*im.x+x]=yy;
		      ZZ[z*xysize+y*im.x+x]=zz;
		    }
		}
	}

      /* {{{ COMMENT save perp image */

#ifdef FoldingComment

{
  char tmpname[10000];
  image_struct tmpim=im;

  tmpim.i=calloc(size*3,sizeof(float));
  sprintf(tmpname,"%s_flowsmooth",argv[argc-1]);

  tmpim.dt=DT_FLOAT;
  tmpim.t=3;

  for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    {
      float tmpX=XX[z*xysize+y*im.x+x],
	tmpY=YY[z*xysize+y*im.x+x],
	tmpZ=ZZ[z*xysize+y*im.x+x],
	tmpf=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);

      if (tmpf>0)
	{
	  tmpim.i[          z*xysize+y*im.x+x]=tmpX/tmpf;
	  tmpim.i[  xyzsize+z*xysize+y*im.x+x]=tmpY/tmpf;
	  tmpim.i[2*xyzsize+z*xysize+y*im.x+x]=tmpZ/tmpf;
	}
    }

  avw_write(tmpname,tmpim);
  free(tmpim.i);
}

#endif

/* }}} */

/* }}} */
      /* {{{ do non-max-suppression in the direction of perp */

      for(z=MAXZ;z<im.z-MAXZ;z++) for(y=1;y<im.y-1;y++) for(x=1;x<im.x-1;x++)
	{
	  FDT theval = in[0][t*xyzsize+z*xysize+y*im.x+x];
	  int xxx=XX[z*xysize+y*im.x+x];
	  int yyy=YY[z*xysize+y*im.x+x];
	  int zzz=ZZ[z*xysize+y*im.x+x];
	  
	  if ( ( (xxx!=0) || (yyy!=0) || (zzz!=0) ) &&
	       ( theval >= in[0][t*xyzsize+(z+zzz)*xysize+(y+yyy)*im.x+x+xxx] ) &&
	       ( theval >  in[0][t*xyzsize+(z-zzz)*xysize+(y-yyy)*im.x+x-xxx] ) &&
	       ( theval >= in[0][t*xyzsize+(z+2*zzz)*xysize+(y+2*yyy)*im.x+x+2*xxx] ) &&
	       ( theval >  in[0][t*xyzsize+(z-2*zzz)*xysize+(y-2*yyy)*im.x+x-2*xxx] ) )
	    tmpim[t*xyzsize+z*xysize+y*im.x+x] = theval;
	}

/* }}} */

      /* {{{ COMMENT count number of 3x3x3 neighbours */

#ifdef FoldingComment

#define TUBESEARCH 3
#define TUBETHRESH 2000

{
  FDT *tmpim2=calloc(size,sizeof(FDT)), *tmpim3=calloc(size,sizeof(FDT));

  for(z=TUBESEARCH;z<im.z-TUBESEARCH;z++) for(y=TUBESEARCH;y<im.y-TUBESEARCH;y++) for(x=TUBESEARCH;x<im.x-TUBESEARCH;x++)
    {
      int xx, yy, zz, sum=0, xxx=0, yyy=0, zzz=0;

      if (tmpim[t*xyzsize+z*xysize+y*im.x+x]>TUBETHRESH)
	{
	  for(zz=-TUBESEARCH; zz<=TUBESEARCH; zz++) for(yy=-TUBESEARCH; yy<=TUBESEARCH; yy++) for(xx=-TUBESEARCH; xx<=TUBESEARCH; xx++)
	    if (tmpim[t*xyzsize+(z+zz)*xysize+(y+yy)*im.x+x+xx]>TUBETHRESH)
	      {
		sum++; xxx+=xx; yyy+=yy; zzz+=zz;
	      }

	  if (sum>0)
	    {
	      float sumfactor, cogfactor;

	      cogfactor = sqrt( ((float)(xxx*xxx+yyy*yyy+zzz*zzz)) / (sum*sum) ) / TUBESEARCH; /* between 0 and 1 */
	      cogfactor = MAX(0, 1-cogfactor);

	      sumfactor = ((float)sum) / ((TUBESEARCH*2+1)*(TUBESEARCH*2+1)); /* generally between 0 and 2; the lower the more tubelike */
	      sumfactor = 0.5 * MAX(0, 2-sumfactor);
	      
	      tmpim2[t*xyzsize+z*xysize+y*im.x+x] = (FDT)(100 * sumfactor );
	    }
	}
    }

  free(tmpim);
  tmpim=tmpim2;
}

#endif

/* }}} */
      /* {{{ count number of 3x3x3 neighbours */

#define TUBESEARCH 2
#define TUBETHRESH 2000

{
  FDT *tmpim2=calloc(size,sizeof(FDT)), *tmpim3=calloc(size,sizeof(FDT));

  for(z=TUBESEARCH;z<im.z-TUBESEARCH;z++) for(y=TUBESEARCH;y<im.y-TUBESEARCH;y++) for(x=TUBESEARCH;x<im.x-TUBESEARCH;x++)
    {
      int xx, yy, zz, sum=0;

      if (tmpim[t*xyzsize+z*xysize+y*im.x+x]>TUBETHRESH)
	{
	  tmpim2[t*xyzsize+z*xysize+y*im.x+x] = 1;

	  for(zz=-TUBESEARCH; zz<=TUBESEARCH; zz++) for(yy=-TUBESEARCH; yy<=TUBESEARCH; yy++) for(xx=-TUBESEARCH; xx<=TUBESEARCH; xx++)
	    if (tmpim[t*xyzsize+(z+zz)*xysize+(y+yy)*im.x+x+xx]>TUBETHRESH)
	      sum++;

	  if (sum>(TUBESEARCH*2+1)*(TUBESEARCH*2+1)*2/3)
	    tmpim2[t*xyzsize+z*xysize+y*im.x+x] = 2;
	}
    }

  for(z=TUBESEARCH;z<im.z-TUBESEARCH;z++) for(y=TUBESEARCH;y<im.y-TUBESEARCH;y++) for(x=TUBESEARCH;x<im.x-TUBESEARCH;x++)
    {
      int xx, yy, zz, sum=0;

      tmpim3[t*xyzsize+z*xysize+y*im.x+x] = tmpim2[t*xyzsize+z*xysize+y*im.x+x];

      if (tmpim2[t*xyzsize+z*xysize+y*im.x+x]>1)
	{
	  for(zz=-TUBESEARCH; zz<=TUBESEARCH; zz++) for(yy=-TUBESEARCH; yy<=TUBESEARCH; yy++) for(xx=-TUBESEARCH; xx<=TUBESEARCH; xx++)
	    if (tmpim2[t*xyzsize+(z+zz)*xysize+(y+yy)*im.x+x+xx]>1)
	      sum++;

	  if (sum>(TUBESEARCH*2+1)*(TUBESEARCH*2+1)*2/3)
	    tmpim3[t*xyzsize+z*xysize+y*im.x+x] = 3;
	}
    }

  free(tmpim);
  tmpim=tmpim3;
}

/* }}} */

      free(X); free(Y); free(Z); free(XX); free(YY); free(ZZ);
    }

  free(in[0]);
  in[0]=tmpim;
}

#endif

/* }}} */
}

/* }}} */
    else if (!strncmp(argv[i], "-nanm", 5))
      /* {{{ NaN mask */

{
  for(j=0;j<size;j++)
    if ( finite((double)in[0][j]) )
      in[0][j]=0;
    else
      in[0][j]=1;
}

/* }}} */
    else if (!strncmp(argv[i], "-nan", 4))
      /* {{{ delete NaNs */

{
  for(j=0;j<size;j++)
    if ( ! finite((double)in[0][j]) )
      in[0][j]=0;
}

/* }}} */
    else if (!strncmp(argv[i], "-roi", 4))
      /* {{{ ROI masking */

{
  int xmin,xsize,ymin,ysize,zmin,zsize,tmin,tsize;

  i++; xmin  = atoi(argv[i]);
  i++; xsize = atoi(argv[i]);
  i++; ymin  = atoi(argv[i]);
  i++; ysize = atoi(argv[i]);
  i++; zmin  = atoi(argv[i]);
  i++; zsize = atoi(argv[i]);
  i++; tmin  = atoi(argv[i]);
  i++; tsize = atoi(argv[i]);

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    if ( (t<tmin) || (t>=tmin+tsize) ||
	 (z<zmin) || (z>=zmin+zsize) ||
	 (y<ymin) || (y>=ymin+ysize) ||
	 (x<xmin) || (x>=xmin+xsize) )
      in[0][t*xyzsize+z*xysize+y*im.x+x]=0;
}

/* }}} */
    else
      usage();
  }

  im.i=in[0];
  find_thresholds(&im,0.1);
  /* if any maths was done then who knows what the stats would be now... */
  if (argc > 3) { im.intent_code = NIFTI_INTENT_NONE; }
  avw_write(argv[i],im);
  return(0);
}

/* }}} */

