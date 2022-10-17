/* {{{ Copyright etc. */

/*  avwstats.c - basic image stats measurements

    Stephen Smith, FMRIB Image Analysis Group

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

#include "libss/libss.h"
#include "libss/libavw.h"

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("\nUsage: avwstats <input> [options]\n\n");
  printf("-l <lthresh> : set lower threshold\n");
  printf("-u <uthresh> : set upper threshold\n");
  printf("-r           : output <robust min intensity> <robust max intensity>\n");
  printf("-R           : output <min intensity> <max intensity>\n");
  printf("-e           : output mean entropy ; mean(-i*ln(i))\n");
  printf("-v           : output <voxels> <volume>\n");
  printf("-V           : output <voxels> <volume> (for nonzero voxels)\n");
  printf("-m           : output mean\n");
  printf("-M           : output mean (for nonzero voxels)\n");
  printf("-s           : output standard deviation\n");
  printf("-S           : output standard deviation (for nonzero voxels)\n");
  printf("-w           : output smallest ROI <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> containing nonzero voxels\n");
  printf("-x           : output co-ordinates of maximum voxel\n\n");
  printf("Note - thresholds are not inclusive ie lthresh<allowed<uthresh\n\n");
  exit(1);
}

/* }}} */
/* {{{ main */

int main(argc, argv)
  int   argc;
  char  *argv [];
{
  /* {{{ vars */

image_struct im;
int i;

/* }}} */

  if (argc<3)
    usage();

  avw_read(argv[1],&im);

  im.thresh=im.min=im.max=0; /* flag up that the thresholds need finding */

  for (i = 2; i < argc; i++) {
    if (!strcmp(argv[i], "-l"))
      /* {{{ set lower threshold */

{
  i++;
  if (argc<i+1)
    {
      printf("Error: no value given following -l\n");
      usage();
    }

  im.lthresh=atof(argv[i]);
  im.thresh=im.min=im.max=0; /* flag up that the thresholds need re-finding */
}

/* }}} */
    else if (!strcmp(argv[i], "-u"))
      /* {{{ set upper threshold */

{
  i++;
  if (argc<i+1)
    {
      printf("Error: no value given following -u\n");
      usage();
    }

  im.uthresh=atof(argv[i]);
  im.thresh=im.min=im.max=0; /* flag up that the thresholds need re-finding */
}

/* }}} */
    else if (!strcmp(argv[i], "-r"))
      /* {{{ robust range */

{
  if (im.thresh==im.min)
    find_thresholds(&im,0.1); /*set all these back to 0.1*/

  printf("%f %f ",(double)im.thresh2,(double)im.thresh98);
}

/* }}} */
    else if (!strcmp(argv[i], "-R"))
      /* {{{ complete range */

{
  if (im.thresh==im.min)
    find_thresholds(&im,0.1);

  printf("%f %f ",(double)im.min,(double)im.max);
}

/* }}} */
    else if (!strcmp(argv[i], "-e"))
      /* {{{ mean entropy */

/* warning - hist uses >=lthresh etc not >lthresh */

#define HISTOGRAM_BINS 1000

{
  int hist[HISTOGRAM_BINS], j, count;
  double entropy=0;

  if (im.thresh==im.min)
    find_thresholds(&im,0.1);

  count=find_histogram(&im,hist,HISTOGRAM_BINS);

  for(j=0;j<HISTOGRAM_BINS;j++)
    entropy -= ((double)hist[j]/count) * log((double)hist[j]/count);

  printf("%f ",entropy/log((double)HISTOGRAM_BINS));
}

/* }}} */
    else if (!strcmp(argv[i], "-v"))
      /* {{{ voxels/volume between thresholds */

{
  int j, count=0;

  for(j=0;j<im.x*im.y*im.z*im.t;j++)
    if ( (im.i[j]>im.lthresh) && (im.i[j]<im.uthresh) )
      count++;

  printf("%d %f ",count,count*im.xv*im.yv*im.zv);
}

/* }}} */
    else if (!strcmp(argv[i], "-V"))
      /* {{{ nonzero voxels/volume between thresholds */

{
  int j, count=0;

  for(j=0;j<im.x*im.y*im.z*im.t;j++)
    if ( (im.i[j]>im.lthresh) && (im.i[j]<im.uthresh) && (im.i[j]!=0) )
      count++;

  printf("%d %f ",count,count*im.xv*im.yv*im.zv);
}

/* }}} */
    else if (!strcmp(argv[i], "-m"))
      /* {{{ mean */

{
  int j, count=0;
  float sum=0;

  for(j=0;j<im.x*im.y*im.z*im.t;j++)
    if ( (im.i[j]>im.lthresh) && (im.i[j]<im.uthresh) )
      {
	sum+=im.i[j];
	count++;
      }

  if (count>0)
    printf("%f ",sum/count);
  else
    printf("0 ");
}

/* }}} */
    else if (!strcmp(argv[i], "-M"))
      /* {{{ nonzero mean */

{
  int j, count=0;
  float sum=0;

  for(j=0;j<im.x*im.y*im.z*im.t;j++)
    if ( (im.i[j]>im.lthresh) && (im.i[j]<im.uthresh) && (im.i[j]!=0) )
      {
	sum+=im.i[j];
	count++;
      }

  if (count>0)
    printf("%f ",sum/count);
  else
    printf("0 ");
}

/* }}} */
    else if (!strcmp(argv[i], "-s"))
      /* {{{ standard deviation */

{
  int j, count=0;
  float sum=0, sumsq=0;

  for(j=0;j<im.x*im.y*im.z*im.t;j++)
    if ( (im.i[j]>im.lthresh) && (im.i[j]<im.uthresh) )
      {
	sum+=im.i[j];
	sumsq+=im.i[j]*im.i[j];
	count++;
      }

  if (count>0)
    {
      float tmpf = (sumsq - sum * sum / count) / count;
      if (tmpf>0) tmpf=sqrt(tmpf);
      else        tmpf=0;
      printf("%f ",tmpf);
    }
  else
    printf("0 ");
}

/* }}} */
    else if (!strcmp(argv[i], "-S"))
      /* {{{ standard deviation */

{
  int j, count=0;
  float sum=0, sumsq=0;

  for(j=0;j<im.x*im.y*im.z*im.t;j++)
    if ( (im.i[j]>im.lthresh) && (im.i[j]<im.uthresh) && (im.i[j]!=0) )
      {
	sum+=im.i[j];
	sumsq+=im.i[j]*im.i[j];
	count++;
      }

  if (count>0)
    {
      float tmpf = (sumsq - sum * sum / count) / count;
      if (tmpf>0) tmpf=sqrt(tmpf);
      else        tmpf=0;
      printf("%f ",tmpf);
    }
  else
    printf("0 ");
}

/* }}} */
    else if (!strcmp(argv[i], "-x"))
      /* {{{ co-ordinates of max voxel */

{
  int x,y,z,t, xm=0,ym=0,zm=0,tm=0;
  float maxi=im.i[0];

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    if ( (im.i[t*im.x*im.y*im.z+z*im.x*im.y+y*im.x+x]>im.lthresh) &&
	 (im.i[t*im.x*im.y*im.z+z*im.x*im.y+y*im.x+x]<im.uthresh) &&
	 (im.i[t*im.x*im.y*im.z+z*im.x*im.y+y*im.x+x]>maxi) )
      {
	maxi=im.i[t*im.x*im.y*im.z+z*im.x*im.y+y*im.x+x];
	xm=x; ym=y; zm=z; tm=t;
      }

  if (t<2)
    printf("%d %d %d ",xm,ym,zm);
  else
    printf("%d %d %d %d ",xm,ym,zm,tm);
}

/* }}} */
    else if (!strcmp(argv[i], "-w"))
      /* {{{ nonzero ROI */

{
  int x,y,z,t,xmin=im.x-1,xmax=0,ymin=im.y-1,ymax=0,zmin=im.z-1,zmax=0,tmin=im.t-1,tmax=0;

  for(t=0;t<im.t;t++) for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
    if ( (im.i[t*im.x*im.y*im.z+z*im.x*im.y+y*im.x+x]>im.lthresh) &&
	 (im.i[t*im.x*im.y*im.z+z*im.x*im.y+y*im.x+x]<im.uthresh) &&
	 (im.i[t*im.x*im.y*im.z+z*im.x*im.y+y*im.x+x]!=0) )
      {
	if (x<xmin) xmin=x;
	if (x>xmax) xmax=x;
	if (y<ymin) ymin=y;
	if (y>ymax) ymax=y;
	if (z<zmin) zmin=z;
	if (z>zmax) zmax=z;
	if (t<tmin) tmin=t;
	if (t>tmax) tmax=t;
      }

  printf("%d %d %d %d %d %d %d %d ",xmin,1+xmax-xmin,ymin,1+ymax-ymin,zmin,1+zmax-zmin,tmin,1+tmax-tmin);
}

/* }}} */
    else
      usage();
  }

  printf("\n");
  return(0);
}

/* }}} */
