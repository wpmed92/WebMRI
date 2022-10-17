/* {{{ Copyright etc. */

/*  avwfill - fill non-brain parts of image with brain data

    Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

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

#define LOOP for(z=MAX(cgz-bz/im.zv,0);z<=MIN(cgz+bz/im.zv,im.z-1);z++) \
               for(y=MAX(cgy-by/im.yv,0);y<=MIN(cgy+by/im.yv,im.y-1);y++) \
                 for(x=MAX(cgx-bx/im.xv,0);x<=MIN(cgx+bx/im.xv,im.x-1);x++) \
                   if (mask[z*im.y*im.x+y*im.x+x]==0) \
                     ok=0

/* }}} */
/* {{{ usage */

void usage(char *progname)
{
  printf("\nUsage: %s <input> <half paradigm period (scans)> <output>\n\n",progname);
  exit(1);
}

/* }}} */
/* {{{ main(argc, argv) */

int main(argc, argv)
  int   argc;
  char  *argv [];
{
/* {{{ vars */

unsigned char *mask, *mask2;
int          x,y,z,t,P,HP,ok;
double       bx,by,bz,cgx,cgy,cgz,scale;
image_struct im,im2;

/* }}} */

  /* {{{ process initial arguments */

if (argc<4)
     usage(argv[0]);

avw_read(argv[1],&im);

im2=im;
im2.i=malloc(sizeof(FDT)*im.x*im.y*im.z*im.t);
mask=malloc(im.x*im.y*im.z);
mask2=malloc(im.x*im.y*im.z);

HP=atoi(argv[2]);
P=2*HP;

/* }}} */
  /* {{{ find initial threshold and CofG */

find_thresholds (&im, 0.2);
printf("thresholds: 0=%f 2=%f t=%f 98=%f 100=%f\n",im.min,im.thresh2,im.thresh,im.thresh98,im.max);

c_of_g (im,&cgx,&cgy,&cgz);
printf("CofG (%f,%f,%f) voxels\n",cgx,cgy,cgz);

/* }}} */
  /* {{{ create binary mask and "close" */

for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
  if (im.i[im.t/2*im.z*im.y*im.x+z*im.y*im.x+y*im.x+x]<im.thresh)
     mask[z*im.y*im.x+y*im.x+x]=0;
  else
     mask[z*im.y*im.x+y*im.x+x]=1;

for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
{
  int xx, yy, zz;
  mask2[z*im.y*im.x+y*im.x+x]=0;
  for(zz=MAX(z-1,0);zz<=MIN(z+1,im.z-1);zz++)
    for(yy=MAX(y-1,0);yy<=MIN(y+1,im.y-1);yy++)
      for(xx=MAX(x-1,0);xx<=MIN(x+1,im.x-1);xx++)
	if (mask[zz*im.y*im.x+yy*im.x+xx]==1)
	  mask2[z*im.y*im.x+y*im.x+x]=1;
}

for(z=0;z<im.z;z++) for(y=0;y<im.y;y++) for(x=0;x<im.x;x++)
{
  int xx, yy, zz;
  mask[z*im.y*im.x+y*im.x+x]=1;
  for(zz=MAX(z-1,0);zz<=MIN(z+1,im.z-1);zz++)
    for(yy=MAX(y-1,0);yy<=MIN(y+1,im.y-1);yy++)
      for(xx=MAX(x-1,0);xx<=MIN(x+1,im.x-1);xx++)
	if (mask2[zz*im.y*im.x+yy*im.x+xx]==0)
	  mask[z*im.y*im.x+y*im.x+x]=0;
}

/*for(x=0;x<im.x*im.y*im.z;x++) im.i[x]=mask[x]; im.t=1; avw_write("MASK",im);exit(0);*/

/* }}} */
  /* {{{ find biggest box */

for(bx=by=bz=10, scale=0.75; scale<0.9999; scale=sqrt(sqrt(scale)))
{
  for(ok=1; bx<1000 && ok; bx++) LOOP; bx=MAX(1,scale*bx);
  for(ok=1; by<1000 && ok; by++) LOOP; by=MAX(1,scale*by);
  for(ok=1; bz<1000 && ok; bz++) LOOP; bz=MAX(1,scale*bz);
  printf("scale=%f bx=%f by=%f bz=%f volume=%f\n",scale,bx,by,bz,bx*by*bz*8);
}

/* }}} */
  /* {{{ use box */

for(t=0;t<im.t/2;t++)
{
  int donedebug=0;
  
  for(z=0;z<im.z;z++)
    for(y=0;y<im.y;y++)
      for(x=0;x<im.x;x++)
	{
	  int xi,yi,zi,toldB,toldNB,inbox;

	  inbox = (z>=cgz-bz/im.zv) && (z<cgz+bz/im.zv) &&
	    (y>=cgy-by/im.yv) && (y<cgy+by/im.yv) &&
	    (x>=cgx-bx/im.xv) && (x<cgx+bx/im.xv);

	  toldB  = (t/P)*2*P+t%P;
	  if (t%P<HP)
	    toldNB = toldB;
	  else
	    toldNB = toldB+HP;

	  if(!donedebug)
	    {
	      donedebug=1;
	      printf("[%d %d %d] ",t,toldB,toldNB);
	    }

	  /* treated as within-brain */
	  if ( (im.i[z*im.y*im.x+y*im.x+x]>=im.thresh) || (inbox) )
	    im2.i[t*im.z*im.y*im.x+z*im.y*im.x+y*im.x+x]=
	      im.i[toldB*im.z*im.y*im.x+z*im.y*im.x+y*im.x+x];
	  
	  /* treated as not-within-brain */
	  if ( inbox )
	    {
	      for(zi=-4;zi<=4;zi++)
		for(yi=-4;yi<=4;yi++)
		  for(xi=-4;xi<=4;xi++)

		    if ( !((xi==0)&&(yi==0)&&(zi==0)) &&
			 (z+zi*(int)(2*bz/im.zv)>=0) && (z+zi*(int)(2*bz/im.zv)<im.z) &&
			 (y+yi*(int)(2*by/im.yv)>=0) && (y+yi*(int)(2*by/im.yv)<im.y) &&
			 (x+xi*(int)(2*bx/im.xv)>=0) && (x+xi*(int)(2*bx/im.xv)<im.x) &&
			 (im.i[(z+zi*(int)(2*bz/im.zv))*im.y*im.x+(y+yi*(int)(2*by/im.yv))*im.x+(x+xi*(int)(2*bx/im.xv))]<im.thresh) )
		      im2.i[t*im.z*im.y*im.x+(z+zi*(int)(2*bz/im.zv))*im.y*im.x+(y+yi*(int)(2*by/im.yv))*im.x+(x+xi*(int)(2*bx/im.xv))]=
			im.i[toldNB*im.z*im.y*im.x+z*im.y*im.x+y*im.x+x];
	    }
	}
}

im.t/=2;
im.i=im2.i;

printf("\n");

/* }}} */
  /* {{{ output and exit */

avw_write(argv[3],im);

return(0);

/* }}} */
}

/* }}} */
