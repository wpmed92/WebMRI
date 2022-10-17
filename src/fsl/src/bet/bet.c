/* {{{ Copyright etc. */

/*  bet.c - Brain Extraction Tool

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
/* {{{ defines, includes and typedefs */

#include "libss/libss.h"
#include "libss/libavw.h"
#include "libss/libtessa.h"

void usage(void);

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("\nBET (Brain Extraction Tool) v1.2 - FMRIB Analysis Group, Oxford\n\n");
  printf("Usage: bet <input fileroot> <output fileroot> [options]\n\n");

  printf("-o : generate brain surface outline overlaid onto original image\n");
  printf("-m : generate binary brain mask\n");
  printf("-s : generate approximate skull image\n");
  printf("-n : don't generate segmented brain image output\n");
  printf("-f <fractional_threshold> : fractional intensity threshold (0->1); default=0.5; smaller values give larger brain outline estimates\n");
  printf("-g <fractional_threshold> : vertical gradient in fractional intensity threshold (-1->1); default=0; positive values give larger brain outline at bottom, smaller at top\n");
  printf("-t : apply thresholding to segmented brain image and mask\n");
  printf("-c <x y z> : co-ordinates (mm not voxels) for centre of initial brain surface sphere\n");
  printf("-r <r>: head radius (mm not voxels); initial surface sphere is set to half of this\n");
  printf("-v : verbose text output\n\n");

  /*printf("-X : generate xtopol format output; file extensions=coo,dat\n");
  printf("-C : generate approximate image cost function output; fileroot=<output filroot>_cost\n");
  printf("-S : as above but \"colour code\" brain-skull distance\n");*/

  exit(1);
}

/* }}} */
/* {{{ main */

#define TESSELATE_ORDER             5       /* 5 */
#define ITERATIONS                  1000    /* 1000 */
#define BRAIN_THRESHOLD_DEFAULT     0.5     /* 0.5 */
#define COST_SEARCH                 7       /* 7 mm */
#define NORMAL_MAX_UPDATE_FRACTION  0.5     /* 0.5 */
#define LAMBDA_TANGENT              0.5     /* 0.5 */
#define LAMBDA_FIT                  0.1     /* 0.1 */
#define RADIUSMAX                   10.0    /* 10.0 mm */
#define RADIUSMIN                   3.33    /* 3.33 mm */
#define SELF_INTERSECTION_THRESHOLD 4000    /* 4000 */
#define SKULL_SEARCH                30      /* 30 mm */
#define SKULL_START                 -3      /* -3 mm */
#define MIN_FOV                     20      /* 20 mm */

/*#define DEBUG_NORMALS*/
/*#define DEBUG_MOVEMENT*/
/*#define DEBUG_EVOLVE*/

int main(argc, argv)
  int   argc;
  char  *argv [];
{
  /* {{{ vars */

FDT *in, *mask=NULL, *raw=NULL, threshold, thresh2, thresh98,
  hist_min=0, hist_max=0, medianval;
int x_size, y_size, z_size, t_size, x, y, z, i, pc=0, iters, pass=1,
  output_brain=1, output_xtopol=0, output_cost=0, output_mask=0,
  output_overlay=0, output_skull=0, apply_thresholding=0, code_skull=0, verbose=0;
double cgx=-1000000, cgy, cgz, radius=-1000000, scale, ml=0, ml0=0, tmpf,
  brain_threshold=BRAIN_THRESHOLD_DEFAULT, gradthresh=0, incfactor=0,
  rE = 0.5 * (1/RADIUSMIN + 1/RADIUSMAX), rF = 6 / (1/RADIUSMIN - 1/RADIUSMAX);
char filename[1000];
image_struct im;
points_struc *v;

/* }}} */

  /* {{{ process arguments */

if (argc<3)
     usage();

avw_read(argv[1],&im);

in=im.i;
x_size=im.x; y_size=im.y; z_size=im.z;
t_size=MAX(im.t,1); im.t=1;
scale=MIN(im.xv,MIN(im.yv,im.zv)); /* minimum voxel side length */

for (i = 3; i < argc; i++) {
  if (!strcmp(argv[i], "-X"))
    output_xtopol=1;
  else if (!strcmp(argv[i], "-C"))
    output_cost=1;
  else if (!strcmp(argv[i], "-s"))
    output_skull=1;
  else if (!strcmp(argv[i], "-S"))
    /* {{{ colour-coded skull image */

    {
      output_skull=1;
      code_skull=1;
    }

/* }}} */
  else if (!strcmp(argv[i], "-m"))
    output_mask=1;
  else if (!strcmp(argv[i], "-o"))
    output_overlay=1;
  else if (!strcmp(argv[i], "-n"))
    output_brain=0;
  else if (!strcmp(argv[i], "-t"))
    apply_thresholding=1;
  else if (!strcmp(argv[i], "-v"))
    verbose=1;
  else if (!strcmp(argv[i], "-f"))
    /* {{{ fractional brain threshold */

    {
      i++;
      if (argc<i+1) /* option following f hasn't been given */
      {
	printf("Error: no value given following -f\n");
	usage();
      }
      brain_threshold=atof(argv[i]);
      if ( (brain_threshold<=0) || (brain_threshold>=1) )
      {
	printf("Error: value following -f must lie between 0 and 1\n");
	usage();
      }
    }

/* }}} */
  else if (!strcmp(argv[i], "-g"))
    /* {{{ gradient fractional brain threshold */

    {
      i++;
      if (argc<i+1) /* option following g hasn't been given */
      {
	printf("Error: no value given following -g\n");
	usage();
      }
      gradthresh=atof(argv[i]);
      if ( (gradthresh<-1) || (gradthresh>1) )
      {
	printf("Error: value following -g must lie between -1 and 1\n");
	usage();
      }
    }

/* }}} */
  else if (!strcmp(argv[i], "-c"))
    /* {{{ initial sphere centre */

    {
      i++;
      if (argc<i+3) /* options following -c hasn't been given */
      {
	printf("Error: need three values after -c\n");
	usage();
      }
      cgx=atof(argv[i]);
      i++;
      cgy=atof(argv[i]);
      i++;
      cgz=atof(argv[i]);
    }

/* }}} */
  else if (!strcmp(argv[i], "-r"))
    /* {{{ initial sphere radius */

    {
      i++;
      if (argc<i+1) /* option following -r hasn't been given */
      {
	printf("Error: no value given following -r\n");
	usage();
      }
      radius=atof(argv[i]);
      if (radius<0)
      {
	printf("Error: value following -r must be greater than 0\n");
	usage();
      }
    }

/* }}} */
  else
    usage();
}

brain_threshold=pow(brain_threshold,0.275);

if ( !output_xtopol && !output_cost && !output_skull && !output_mask && !output_overlay && !output_brain )
{
  printf("No outputs requested!\n");
  usage();
}

/* }}} */
  /* {{{ image preprocessing */

/* {{{ setup for single-slice data if present */

if (im.x*im.xv<MIN_FOV) im.xv=200;
if (im.y*im.yv<MIN_FOV) im.yv=200;
if (im.z*im.zv<MIN_FOV) im.zv=200;

/* }}} */

if (verbose) printf("\nBET (Brain Extraction Tool) v1.2 - FMRIB Analysis Group, Oxford\n\n");

im.min=im.max=0;
find_thresholds(&im,0.1);
hist_min=im.min; hist_max=im.max; thresh2=im.thresh2; thresh98=im.thresh98; threshold=im.thresh;
if (verbose) printf("hist_min=%.2f thresh2=%.2f thresh=%.2f thresh98=%.2f hist_max=%.2f\n",
       im.min,(double)thresh2,(double)threshold,(double)thresh98,im.max);
if (verbose) printf("THRESHOLD %.1f\n",(double)threshold);

if (cgx<-999999)
{
  c_of_g (im,&cgx,&cgy,&cgz);
  cgx*=im.xv; cgy*=im.yv; cgz*=im.zv;
}
if (verbose) printf("CofG (%.1f,%.1f,%.1f) mm\n",cgx,cgy,cgz);

if (radius<-999999)
     radius = find_radius (im,im.xv*im.yv*im.zv);
if (verbose) printf("RADIUS %.1f\n",radius);

if ( im.thresh98 - im.thresh2 > 1.5 )
{
  FDT *tmpimage = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size);

  i=0;
  for(z=0; z<z_size; z++)
    for(y=0; y<y_size; y++)
      for(x=0; x<x_size; x++)
	{
	  FDT tmp=IA(in,x,y,z);
	  
	  if ( (tmp>thresh2) && (tmp<thresh98) &&
	       ( (x*im.xv-cgx)*(x*im.xv-cgx) + (y*im.yv-cgy)*(y*im.yv-cgy) + (z*im.zv-cgz)*(z*im.zv-cgz) < radius*radius ) )
	    tmpimage[i++]=tmp;
	}
  medianval = median(0.5,tmpimage,i);
  if (verbose)   printf("MEDIANVAL %.1f\n",(double)medianval);
  free(tmpimage);
} else {
  medianval=im.thresh98;
  if (verbose)   printf("MEDIANVAL %.1f\n",(double)medianval);
}

if (output_cost) /* prepare cost function image for writing into */
{
  raw = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size);
  memset(raw,(unsigned char)0,sizeof(FDT)*x_size*y_size*z_size);
}

/* }}} */

  while (pass>0)
    {
      /* {{{ initialize tessellation */

if (pass==1)
{
  object *old = &ico;                /* start with icosohedral tessellation */
  
  tessa((int)TESSELATE_ORDER,&old);  /* create tessellated triangles; level 5 gives 2562 points */
  pc=points_list(old,&v);            /* convert triangles to vertex list */
  free(old);                         /* don't need triangular tessellation any more */

  /*printf("VERTICES %d\n",pc);*/

  /* measure initial spacing for use later in self-intersection calculations */
  ml0 = sqrt( (v[0].xorig-v[v[0].n[0]].xorig)*(v[0].xorig-v[v[0].n[0]].xorig) +
	      (v[0].yorig-v[v[0].n[0]].yorig)*(v[0].yorig-v[v[0].n[0]].yorig) +
	      (v[0].zorig-v[v[0].n[0]].zorig)*(v[0].zorig-v[v[0].n[0]].zorig) );
  /*printf("ml0=%f\n",ml0);*/

#ifdef DEBUG_NORMALS
  v = (points_struc *)realloc((void *)v,sizeof(points_struc)*2*(pc+10)); /* make space for storing normals */
#endif

}

/* scale vertex positions for this image, set surface small and allow to grow */
for(i=0; i<pc; i++)
{
  v[i].x = v[i].xorig * radius * 0.5 + cgx;
  v[i].y = v[i].yorig * radius * 0.5 + cgy;
  v[i].z = v[i].zorig * radius * 0.5 + cgz;
}

/* }}} */
      /* {{{ find brain surface */

for(iters=0; iters<ITERATIONS; iters++)
{
  /* {{{ find local surface normals */

for(i=0; i<pc; i++)
{
  double nx, ny, nz, tmpf;
  int k, l;

  nx=ny=nz=0.0;

  for(k=0; v[i].n[k]>-1; k++); /* find number of connections */

  for(l=0; l<k; l++) /* for each pair of consecutive neighbours form a vector product to get normal */
    {
      double adx = v[v[i].n[l]].x - v[i].x,
	ady = v[v[i].n[l]].y - v[i].y,
	adz = v[v[i].n[l]].z - v[i].z,
	bdx = v[v[i].n[(l+1)%k]].x - v[i].x,
	bdy = v[v[i].n[(l+1)%k]].y - v[i].y,
	bdz = v[v[i].n[(l+1)%k]].z - v[i].z;
  
      nx += ady*bdz - adz*bdy;
      ny += adz*bdx - adx*bdz;
      nz += adx*bdy - ady*bdx;
    }

  /* make the normal vector of length 1 */
  tmpf = sqrt(nx*nx+ny*ny+nz*nz);
  nx/=(double)tmpf; ny/=(double)tmpf; nz/=(double)tmpf;

  /* {{{ debug normals */

#ifdef DEBUG_NORMALS

if (iters==ITERATIONS-1) /* final iteration */
{
     v[pc+i].x=v[i].x+nx*20.0;
     v[pc+i].y=v[i].y+ny*20.0;
     v[pc+i].z=v[i].z+nz*20.0;
     v[pc+i].n[0]=i;
     v[pc+i].n[1]=-1;
}

#endif

/* }}} */

  v[i].nx=nx;
  v[i].ny=ny;
  v[i].nz=nz;
}

/* }}} */
  /* {{{ find mean connection length every now and then */

if ( (iters==50) || (iters%100==0) ) /* add the 50 as the rate of change is highest at start; thus do 0,50,100,200,..... */
{
  int l;

  ml=0;

  for(i=0; i<pc; i++)
    {
      double mml=0;
      
      for(l=0; v[i].n[l]>-1; l++)
	mml += sqrt( (v[i].x-v[v[i].n[l]].x)*(v[i].x-v[v[i].n[l]].x) +
		     (v[i].y-v[v[i].n[l]].y)*(v[i].y-v[v[i].n[l]].y) +
		     (v[i].z-v[v[i].n[l]].z)*(v[i].z-v[v[i].n[l]].z) );

      ml += mml/l;
    }

  ml /= pc;
  /*printf("iteration=%d, ml=%f\n",iters,ml);*/
}

/* }}} */
  /* {{{ increased smoothing for pass>1 */

if (pass>1)
{
  incfactor=pow((double)10.0,(double)pass);

  if (iters>ITERATIONS*0.75)
    incfactor=4.0*(1.0-((double)iters)/ITERATIONS)*(incfactor-1.0) + 1.0;
}

/* }}} */

  for(i=0; i<pc; i++) /* calculate tessellation update */
    {
      /* {{{ variables, and setup k and normal */

FDT lmin, lmax;
int d, k, l;
double nx=v[i].nx, ny=v[i].ny, nz=v[i].nz, sx, sy, sz, fit, sn, stx, sty, stz;

for(k=0; v[i].n[k]>-1; k++); /* find number of connections */

/* }}} */
      /* {{{ average position of neighbours: smoothing vector */

/* s is vector from current vertex to the mean position of its neighbours */

sx=sy=sz=0;

for(l=0; l<k; l++)
{
  sx += v[v[i].n[l]].x;
  sy += v[v[i].n[l]].y;
  sz += v[v[i].n[l]].z;
}

sx = sx/k - v[i].x;
sy = sy/k - v[i].y;
sz = sz/k - v[i].z;

/* part of s normal to surface, sn = n * (s.n)
   part of s tangent to surface, st = s - sn */

sn = sx*nx + sy*ny + sz*nz; /* this is just the s.n part - will multiply by n later */

stx = sx - nx * sn;
sty = sy - ny * sn;
stz = sz - nz * sn;

/* }}} */
      /* {{{ COMMENT ORIG find intensity-based part of cost function - local max */

#ifdef FoldingComment

{
  FDT lnew;
  int all_inside_image=1;
  double lthresh, local_brain_threshold;

  fit=0;
  lmin=hist_max;
  lmax=thresh98/2+thresh2/2;

  tmpf = brain_threshold + gradthresh * ( (v[i].z-cgz) / radius );
  local_brain_threshold = MIN( 1.0 , MAX( tmpf , 0.0 ) );  

  d=1;   /* start d at 1 not 0 so that boundary is just outside brain not on the actual edge */
  x=FTOI((v[i].x-((double)d)*nx)/im.xv); y=FTOI((v[i].y-((double)d)*ny)/im.yv); z=FTOI((v[i].z-((double)d)*nz)/im.zv);
  if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) && (z>=0) && (z<z_size) )
    {
      lnew=IA(in,x,y,z);
      lmin = MIN(lmin,lnew);
      lmax = MAX(lmax,lnew);
    }
  else all_inside_image=0;

  d=COST_SEARCH;
  x=FTOI((v[i].x-((double)d)*nx)/im.xv); y=FTOI((v[i].y-((double)d)*ny)/im.yv); z=FTOI((v[i].z-((double)d)*nz)/im.zv);
  if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) && (z>=0) && (z<z_size) )
    {
      lnew=IA(in,x,y,z);
      lmin = MIN(lmin,lnew);
      /*lmax = MAX(lmax,lnew);*/
    }
  else all_inside_image=0;

  if (all_inside_image)
    {
      for(d=2;d<COST_SEARCH;d++)
	{
	  x=FTOI((v[i].x-((double)d)*nx)/im.xv); y=FTOI((v[i].y-((double)d)*ny)/im.yv); z=FTOI((v[i].z-((double)d)*nz)/im.zv);
	  lnew=IA(in,x,y,z);
	  lmin = MIN(lmin,lnew);
	  if (d<COST_SEARCH/2) /* only look relatively locally for maximum intensity */
	    lmax = MAX(lmax,lnew);
	}

      lmin = MAX(lmin,thresh2);   /* so that extreme values don't screw this up */
      lmax = MIN(lmax,thresh98);

      lthresh = (lmax - thresh2)*local_brain_threshold + thresh2;
      tmpf = lmin - lthresh;
      fit = tmpf / ((lmax - thresh2)*0.5); /* scale range to around -1:1 */
  
      if (output_cost)
	{
	  x=(v[i].x/im.xv); y=(v[i].y/im.yv); z=(v[i].z/im.zv);
	  if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) &&(z>=0) && (z<z_size) )
	    IA(raw,x,y,z)=thresh98/2 + thresh2/2 + (0.5*fit*((double)thresh98 - (double)thresh2));
	}  
    }
}

#endif

/* }}} */
      /* {{{ find intensity-based part of cost function - local max */

{
  FDT lnew;
  int all_inside_image=1;
  double lthresh, local_brain_threshold=brain_threshold;

  fit=0;
  lmin=medianval;
  lmax=threshold;

  /* {{{ change local threshold if gradient threshold used */

  if (gradthresh!=0)
    {
      tmpf = brain_threshold + gradthresh * ( (v[i].z-cgz) / radius );
      local_brain_threshold = MIN( 1.0 , MAX( tmpf , 0.0 ) );
    }

/* }}} */

  d=1;   /* start d at 1 not 0 so that boundary is just outside brain not on the actual edge */
  x=FTOI((v[i].x-((double)d)*nx)/im.xv); y=FTOI((v[i].y-((double)d)*ny)/im.yv); z=FTOI((v[i].z-((double)d)*nz)/im.zv);
  if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) && (z>=0) && (z<z_size) )
    {
      lnew=IA(in,x,y,z);
      lmin = MIN(lmin,lnew);
      lmax = MAX(lmax,lnew);
    }
  else all_inside_image=0;

  d=COST_SEARCH;
  x=FTOI((v[i].x-((double)d)*nx)/im.xv); y=FTOI((v[i].y-((double)d)*ny)/im.yv); z=FTOI((v[i].z-((double)d)*nz)/im.zv);
  if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) && (z>=0) && (z<z_size) )
    {
      lnew=IA(in,x,y,z);
      lmin = MIN(lmin,lnew);
    }
  else all_inside_image=0;

  if (all_inside_image)
    {
      for(d=2;d<COST_SEARCH;d++)
	{
	  x=FTOI((v[i].x-((double)d)*nx)/im.xv); y=FTOI((v[i].y-((double)d)*ny)/im.yv); z=FTOI((v[i].z-((double)d)*nz)/im.zv);
	  lnew=IA(in,x,y,z);
	  lmin = MIN(lmin,lnew);
	  if (d<COST_SEARCH/2) /* only look relatively locally for maximum intensity */
	    lmax = MAX(lmax,lnew);
	}

      lmin = MAX(lmin,thresh2);
      lmax = MIN(lmax,medianval);

      lthresh = (lmax - thresh2)*local_brain_threshold + thresh2;

      tmpf = lmin - lthresh;
      if (lmax - thresh2>0)
	fit = tmpf / ((lmax - thresh2)*0.5); /* scale range to around -1:1 */
      else
	fit = tmpf / 0.5;
  
      if (output_cost)
	{
	  x=(v[i].x/im.xv); y=(v[i].y/im.yv); z=(v[i].z/im.zv);
	  if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) &&(z>=0) && (z<z_size) )
	    IA(raw,x,y,z)=thresh98/2 + thresh2/2 + (0.5*fit*((double)thresh98 - (double)thresh2));
	}  
    }
}

/* }}} */
      /* {{{ estimate the update */

/* normal component of smoothing */
tmpf = 2 * ABS(sn) / (ml*ml);              /* 1/r ; r=local radius of curvature */
tmpf = 0.5 * ( 1 + tanh((tmpf-rE)*rF) );   /* f(1/r) ; this makes big r have low correction and vice versa */

if ( (pass>1) && (sn > 0) )  /* if need increased smoothing, to avoid self-intersection; */
{                            /* only apply to concave curvature (as seen from outside surface) */
  tmpf *= incfactor;
  tmpf = MIN(tmpf,1);
}

tmpf = sn*tmpf;                                              /* multiply normal component by correction factor */

tmpf += NORMAL_MAX_UPDATE_FRACTION * ml * LAMBDA_FIT * fit ; /* combine normal smoothing and image fit */

v[i].xnew = v[i].x + stx*LAMBDA_TANGENT + tmpf*nx;
v[i].ynew = v[i].y + sty*LAMBDA_TANGENT + tmpf*ny;
v[i].znew = v[i].z + stz*LAMBDA_TANGENT + tmpf*nz;

/* }}} */
    }

  /* {{{ debug evolve: output single slice to show how shape evolves */

#ifdef DEBUG_EVOLVE

#define SPACING 5

if (iters%SPACING==0)
{
  image_struct imd=im;
  FILE *ofp;
  FDT *tmpimage=(FDT *)calloc(x_size*y_size*z_size,sizeof(FDT));

  if (iters==0)
    {
      imd.x=im.y; imd.y=im.z; imd.z=ITERATIONS/SPACING; imd.t=1;
      imd.xv0=im.yv; imd.yv0=im.zv; imd.zv0=1;
      imd.min=thresh2; imd.max=thresh98;
      imd.i=tmpimage;

      sprintf(filename,"%s_debug",argv[2]);
      avw_write(filename,imd);

      sprintf(filename,"%s_debug.img",argv[2]);
      ofp=fopen(filename,"wb");
    }
  else
    {
      sprintf(filename,"%s_debug.img",argv[2]);
      ofp=fopen(filename,"ab");
    }

  memcpy(tmpimage,in,sizeof(FDT)*x_size*y_size*z_size);
  im.i=tmpimage;
  draw_surface(&im,hist_max,thresh2,v,pc);
  im.i=in;

  x=x_size*0.52;
  for(z=0;z<z_size;z++)
    for(y=0;y<y_size;y++)
      fwrite(&tmpimage[z*y_size*x_size+y*x_size+x],sizeof(FDT),1,ofp);

  free(tmpimage);
  fclose(ofp);
}

#endif

/* }}} */
  /* {{{ debug movement */

#ifdef DEBUG_MOVEMENT

if ( (iters==50) || (iters%100==0) )
{
  for(mm=0, i=0; i<pc; i++)
    mm += sqrt ( (v[i].x-v[i].xnew)*(v[i].x-v[i].xnew) + (v[i].y-v[i].ynew)*(v[i].y-v[i].ynew) + (v[i].z-v[i].znew)*(v[i].z-v[i].znew) );

  mm /= pc;

  printf("mean movement=%f\n",mm);
}

#endif

/* }}} */
  /* {{{ update tessellation */

for(i=0; i<pc; i++)
{
  v[i].x = v[i].xnew;
  v[i].y = v[i].ynew;
  v[i].z = v[i].znew;
}

/* }}} */
}

/* }}} */
      /* {{{ test surface for non-spherecity */

if (pass<10)
{
  int j, l;
  double intersection=0;

  for(i=0; i<pc; i++) /* loop round all points */
      for(j=0; j<pc; j++) /* inner loop round all points */
	{
	  int do_it=1;
	  
	  if (j==i) /* other point is same as current one - don't use */
	    do_it=0;
      
	  for(l=0; v[i].n[l]>-1; l++) /* other point is connected to current one - don't use */
	    if (j==v[i].n[l])
	      do_it=0;

	  if ( (do_it) &&
	       ( (v[i].x-v[j].x)*(v[i].x-v[j].x) + (v[i].y-v[j].y)*(v[i].y-v[j].y) + (v[i].z-v[j].z)*(v[i].z-v[j].z) < ml*ml ) )
	    {
	      double dist = sqrt ( (v[i].x-v[j].x)*(v[i].x-v[j].x) +
				   (v[i].y-v[j].y)*(v[i].y-v[j].y) +
				   (v[i].z-v[j].z)*(v[i].z-v[j].z) ),
		distorig = sqrt ( (v[i].xorig-v[j].xorig)*(v[i].xorig-v[j].xorig) +
				  (v[i].yorig-v[j].yorig)*(v[i].yorig-v[j].yorig) +
				  (v[i].zorig-v[j].zorig)*(v[i].zorig-v[j].zorig) );

	      tmpf = (distorig/ml0) - (dist/ml);    /* orig distance (in units of mean length) - current distance */
	      tmpf *= tmpf;                         /* weight higher values more highly */
	      /*printf("self-intersection value %f at: %f %f %f\n",tmpf,v[i].x/im.xv,v[i].y/im.yv,v[i].z/im.zv);*/
	      intersection += tmpf;
	    }
	}

  if (verbose) printf("self-intersection total = %.1f (threshold=%.1f)\n",intersection,(double)SELF_INTERSECTION_THRESHOLD);

  if (intersection>SELF_INTERSECTION_THRESHOLD)
    {
      if (verbose) printf("thus will rerun with higher smoothness constraint\n");
      pass++;
    }
  else
    pass=0;
}
else
  pass=0;

/* }}} */
    }

  /* {{{ write brain image and mask */

if ( (output_mask) || (output_brain) )
{
  mask = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size);
  im.i=mask;
  fill_surface(&im,v,pc);

  if (apply_thresholding)
    for(i=0; i<z_size*y_size*x_size; i++)
      if (in[i]<threshold)
	mask[i]=0;
}

if (output_mask)
{
  im.min=0;
  im.max=1;
  sprintf(filename,"%s_mask",argv[2]);
  avw_write(filename,im);
}

if (output_brain)
{
  FDT *tmpin = malloc(sizeof(FDT)*x_size*y_size*z_size*t_size);
  int goesneg=0;

  im.min=thresh2;
  im.max=thresh98;

  if ( (im.dtmin<0) && ((double)hist_min<0) )
    goesneg=1;

  for(i=0; i<z_size*y_size*x_size*t_size; i++)
    {
      if ( goesneg )
	{
	  if (mask[i%(z_size*y_size*x_size)]<0.5)
	    tmpin[i]=0;
	  else
	    tmpin[i]=(FDT)( ((double)(MIN(MAX(in[i],thresh2),thresh98))-(double)thresh2) * (im.dtmax-1) / ((double)thresh98-(double)thresh2) + 1 );

	  im.min=0;
	  im.max=(FDT)im.dtmax;
	}
      else
	tmpin[i]=mask[i%(z_size*y_size*x_size)]*in[i];
    }

  im.i=tmpin;
  im.t=t_size;
  avw_write(argv[2],im);
  im.t=1;
  free(tmpin);
}

if ( (output_mask) || (output_brain) )
  free(mask);

/* }}} */
  /* {{{ output cost */

if (output_cost)
{
  im.i=raw;
  im.min=raw[0];
  im.max=raw[0];

  for(i=0; i<z_size*y_size*x_size; i++)
    {
      if (raw[i]>im.max) im.max=raw[i];
      if (raw[i]<im.min) im.min=raw[i];
    }

  sprintf(filename,"%s_cost",argv[2]);
  avw_write(filename,im);
  free(raw);
}

/* }}} */
  /* {{{ find skull and output */

if (output_skull)
{
  FDT *skull = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size);

  im.i=skull;
  memset(skull,(unsigned char)0,sizeof(FDT)*x_size*y_size*z_size);
  draw_surface(&im,1,1,v,pc); /* tell skull searching where to start */

  /* {{{ create skull image */

im.i=in;

for(z=0; z<z_size; z++)
  for(y=0; y<y_size; y++)
    for(x=0; x<x_size; x++)
      if (IA(skull,x,y,z)==1)
        {
	  FDT val, val2, minval, maxval;
	  double nx, ny, nz, d, d_max=0, d_min, grad, maxgrad, X, Y, Z, maxJ, lastJ;
	  int xx, yy, zz, j=0;

	  /* {{{ zero this point */

IA(skull,x,y,z)=0;

/* }}} */
	  /* {{{ find nearest node and setup normal */

{
  double mind=10000000;

  for(i=0; i<pc; i++)
  { 
    tmpf = (x*im.xv-v[i].x)*(x*im.xv-v[i].x) +
      (y*im.yv-v[i].y)*(y*im.yv-v[i].y) +
      (z*im.zv-v[i].z)*(z*im.zv-v[i].z);
    if (tmpf<mind)
      {
	j=i;
	mind=tmpf;
      }
  }

  nx=v[j].nx;
  ny=v[j].ny;
  nz=v[j].nz;
}

/* }}} */
	  /* {{{ find minval, d_max and maxval up to SKULL_SEARCH distance */

maxval=threshold;
minval=IA(in,x,y,z);

	  for(d=0; d<SKULL_SEARCH; d+=scale*0.5)
	    {
	      xx=x+FTOI(d*nx/im.xv); yy=y+FTOI(d*ny/im.yv); zz=z+FTOI(d*nz/im.zv);

	      if ( (xx>=0) && (xx<x_size) && (yy>=0) && (yy<y_size) && (zz>=0) && (zz<z_size) )
		{
		  val=IA(in,xx,yy,zz);

		  if (val>maxval)
		    {
		      maxval=val;
		      d_max=d;
		    }

		  if (val<minval)
		      minval=val;
		}
	    }

/* }}} */
	  if (maxval>threshold) /* can we see an intensity peak on the far side of the skull? */
	    {
	      /* {{{ find furthest out point that has small intensity */

maxgrad=1;
d_min=SKULL_START;
maxJ =-1000000;
lastJ=-2000000;

	  for(d=SKULL_START; d<d_max; d+=scale*0.5)
	    {
	      xx=x+FTOI(d*nx/im.xv); yy=y+FTOI(d*ny/im.yv); zz=z+FTOI(d*nz/im.zv);

	      if ( (xx>=0) && (xx<x_size) && (yy>=0) && (yy<y_size) && (zz>=0) && (zz<z_size) )
		{
		  tmpf=d/30 - IA(in,xx,yy,zz)/((double)(thresh98-thresh2));

		  /*		  if ( (tmpf>maxJ) && (lastJ+1.0/30>tmpf*0.5) )*/
		  if (tmpf>maxJ)
		  {
		    maxJ=tmpf;
		    d_min=d;
		  }

		  lastJ=tmpf;
		}
	    }

/* }}} */
	      /* {{{ find _first_ max gradient out from d_min */

maxgrad=0;

X=x+d_min*nx/im.xv; Y=y+d_min*ny/im.yv; Z=z+d_min*nz/im.zv;

if ( (X>0) && (X<x_size-1) && (Y>0) && (Y<y_size-1) && (Z>0) && (Z<z_size-1) )
{
  val2 = TLI(im,X,Y,Z);

	  for(d=d_min+scale; d<d_max; d+=scale*0.5)
	    {
	      X=x+d*nx/im.xv; Y=y+d*ny/im.yv; Z=z+d*nz/im.zv;

	      if ( (X>0) && (X<x_size-1) && (Y>0) && (Y<y_size-1) && (Z>0) && (Z<z_size-1) )
		{
		  val = TLI(im,X,Y,Z);
		  /*printf("%d %d %d   %f %f %f   %d %d\n",x,y,z,X,Y,Z,(int)val,(int)val2);*/

		  grad=val-val2;
		  val2=val;		 

		  if (grad>0) /* this so that we don't do anything if we're still in the same voxel */
		    {
		      if (grad > maxgrad)
			{
			  maxgrad=grad;
			  d_min=d;
			}
		      else
			d=d_max;
		    }
		}
	    }
}

/*printf("\n");*/

/* }}} */
	      /* {{{ mark this point as skull */

if (maxgrad>0)
{
  xx=x+FTOI(d_min*nx/im.xv);
  yy=y+FTOI(d_min*ny/im.yv);
  zz=z+FTOI(d_min*nz/im.zv);

  if ( (xx>=0) && (xx<x_size) && (yy>=0) && (yy<y_size) && (zz>=0) && (zz<z_size) )
    {
      if (code_skull)
	IA(skull,xx,yy,zz)=(FDT)(d_min*100.0);
      else
	IA(skull,xx,yy,zz)=100;
    }
}

/* }}} */
	    }
	}

im.i=skull;

/* }}} */
  /* {{{ output */

im.min=0;
im.max=100;
sprintf(filename,"%s_skull",argv[2]);
avw_write(filename,im);

/* }}} */

  free(skull);
}

/* }}} */
  /* {{{ output overlay (corrupts input image) */

if (output_overlay)
{
  im.min=thresh2;
  im.max=hist_max;
  im.i=in;
  draw_surface(&im,im.max,im.max,v,pc);
  sprintf(filename,"%s_overlay",argv[2]);
  avw_write(filename,im);
}

/* }}} */
  /* {{{ output xtopol surface (corrupts tessellation) */

  /* {{{ output xtopol surface (corrupts tessellation) */

if (output_xtopol)
{
  int xtopol_pc=pc;
  double ax=0, ay=0, az=0, tmpf=0;

#ifdef DEBUG_NORMALS
  xtopol_pc *= 2;
#endif

  for(i=0; i<xtopol_pc; i++)
    {
      ax += v[i].x;
      ay += v[i].y;
      az += v[i].z;
    }
  ax/=(double)xtopol_pc; ay/=(double)xtopol_pc; az/=(double)xtopol_pc; 

  for(i=0; i<xtopol_pc; i++)
    tmpf += sqrt( (v[i].x-ax)*(v[i].x-ax) + (v[i].y-ay)*(v[i].y-ay) + (v[i].z-az)*(v[i].z-az) );
  tmpf/=(double)xtopol_pc;

  for(i=0; i<xtopol_pc; i++)
    {
      v[i].x=(v[i].x-ax)/tmpf;;
      v[i].y=(v[i].y-ay)/tmpf;;
      v[i].z=(v[i].z-az)/tmpf;;
    }

  xtopol_output(v,xtopol_pc,argv[2]);
}

/* }}} */
  /* {{{ COMMENT output xtopol surface (x>0) */

#ifdef FoldingComment

if (output_xtopol)
{
  int xtopol_pc=pc, j, k, l, *mapping;
  double ax=0, ay=0, az=0, tmpf=0;

#ifdef DEBUG_NORMALS
  xtopol_pc *= 2;
#endif

  mapping = (int*)malloc(sizeof(int)*xtopol_pc);

  for(i=0; i<xtopol_pc; i++)
    {
      ax += v[i].x;
      ay += v[i].y;
      az += v[i].z;
    }
  ax/=(double)xtopol_pc; ay/=(double)xtopol_pc; az/=(double)xtopol_pc; 

  for(i=0; i<xtopol_pc; i++)
    tmpf += sqrt( (v[i].x-ax)*(v[i].x-ax) + (v[i].y-ay)*(v[i].y-ay) + (v[i].z-az)*(v[i].z-az) );
  tmpf/=(double)xtopol_pc;

  j=0;
  for(i=0; i<xtopol_pc; i++)
    {
      if (v[i].x-ax>0)
	{
	  mapping[i]=j;
	  v[j].x=v[i].x;
	  v[j].y=v[i].y;
	  v[j].z=v[i].z;
	  l=0;
	  for(k=0; v[i].n[k]!=-1; k++)
	    {
	      if (v[v[i].n[k]].x-ax>0)
		{
		  v[j].n[l]=v[i].n[k];
		  l++;
		}
	    }
	  v[j].n[l]=-1;
	  j++;
	}
    }

  for(i=0; i<j; i++)
    {
      v[i].x=(v[i].x-ax)/tmpf;
      v[i].y=(v[i].y-ay)/tmpf;
      v[i].z=(v[i].z-az)/tmpf;
      for(k=0; v[i].n[k]!=-1; k++)
	v[i].n[k] = mapping[v[i].n[k]];
    }

  xtopol_output(v,j,argv[2]);
}

#endif

/* }}} */
  /* {{{ COMMENT output xtopol surface (y>0) */

#ifdef FoldingComment

if (output_xtopol)
{
  int xtopol_pc=pc, j, k, l, *mapping;
  double ax=0, ay=0, az=0, tmpf=0;

#ifdef DEBUG_NORMALS
  xtopol_pc *= 2;
#endif

  mapping = (int*)malloc(sizeof(int)*xtopol_pc);

  for(i=0; i<xtopol_pc; i++)
    {
      ax += v[i].x;
      ay += v[i].y;
      az += v[i].z;
    }
  ax/=(double)xtopol_pc; ay/=(double)xtopol_pc; az/=(double)xtopol_pc; 

  for(i=0; i<xtopol_pc; i++)
    tmpf += sqrt( (v[i].x-ax)*(v[i].x-ax) + (v[i].y-ay)*(v[i].y-ay) + (v[i].z-az)*(v[i].z-az) );
  tmpf/=(double)xtopol_pc;

  j=0;
  for(i=0; i<xtopol_pc; i++)
    {
      if (v[i].y-ay>0)
	{
	  mapping[i]=j;
	  v[j].x=v[i].x;
	  v[j].y=v[i].y;
	  v[j].z=v[i].z;
	  l=0;
	  for(k=0; v[i].n[k]!=-1; k++)
	    {
	      if (v[v[i].n[k]].y-ay>0)
		{
		  v[j].n[l]=v[i].n[k];
		  l++;
		}
	    }
	  v[j].n[l]=-1;
	  j++;
	}
    }

  for(i=0; i<j; i++)
    {
      v[i].x=(v[i].x-ax)/tmpf;
      v[i].y=(v[i].y-ay)/tmpf;
      v[i].z=(v[i].z-az)/tmpf;
      for(k=0; v[i].n[k]!=-1; k++)
	v[i].n[k] = mapping[v[i].n[k]];
    }

  xtopol_output(v,j,argv[2]);
}

#endif

/* }}} */
  /* {{{ COMMENT output xtopol surface (z>0) */

#ifdef FoldingComment

if (output_xtopol)
{
  int xtopol_pc=pc, j, k, l, *mapping;
  double ax=0, ay=0, az=0, tmpf=0;

#ifdef DEBUG_NORMALS
  xtopol_pc *= 2;
#endif

  mapping = (int*)malloc(sizeof(int)*xtopol_pc);

  for(i=0; i<xtopol_pc; i++)
    {
      ax += v[i].x;
      ay += v[i].y;
      az += v[i].z;
    }
  ax/=(double)xtopol_pc; ay/=(double)xtopol_pc; az/=(double)xtopol_pc; 

  for(i=0; i<xtopol_pc; i++)
    tmpf += sqrt( (v[i].x-ax)*(v[i].x-ax) + (v[i].y-ay)*(v[i].y-ay) + (v[i].z-az)*(v[i].z-az) );
  tmpf/=(double)xtopol_pc;

  j=0;
  for(i=0; i<xtopol_pc; i++)
    {
      if (v[i].z-az>0)
	{
	  mapping[i]=j;
	  v[j].x=v[i].x;
	  v[j].y=v[i].y;
	  v[j].z=v[i].z;
	  l=0;
	  for(k=0; v[i].n[k]!=-1; k++)
	    {
	      if (v[v[i].n[k]].z-az>0)
		{
		  v[j].n[l]=v[i].n[k];
		  l++;
		}
	    }
	  v[j].n[l]=-1;
	  j++;
	}
    }

  for(i=0; i<j; i++)
    {
      v[i].x=(v[i].x-ax)/tmpf;
      v[i].y=(v[i].y-ay)/tmpf;
      v[i].z=(v[i].z-az)/tmpf;
      for(k=0; v[i].n[k]!=-1; k++)
	v[i].n[k] = mapping[v[i].n[k]];
    }

  xtopol_output(v,j,argv[2]);
}

#endif

/* }}} */

/* }}} */

  if (verbose) printf("\n");
  return(0);
}

/* }}} */
