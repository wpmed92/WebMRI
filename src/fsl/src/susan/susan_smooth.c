/* {{{ Copyright */

/*

  SUSAN Version 2i_medx - nonlinear 2D/3D smoothing

  Oxford Centre for Functional Magnetic Resonance Imaging of the Brain,
  Department of Clinical Neurology, Oxford University, Oxford, UK
  (Previously in Computer Vision and Image Processing Group - now
  Computer Vision and Electro Optics Group - DERA Chertsey, UK)
  Email:    steve@fmrib.ox.ac.uk
  WWW:      http://www.fmrib.ox.ac.uk/~steve

  (C) Crown Copyright (1995-1999), Defence Evaluation and Research Agency,
  Farnborough, Hampshire, GU14 6TD, UK
  DERA WWW site:
  http://www.dera.gov.uk/
  DERA Computer Vision and Electro Optics Group WWW site:
  http://www.dera.gov.uk/imageprocessing/dera/group_home.html
  DERA Computer Vision and Electro Optics Group point of contact:
  Dr. John Savage, jtsavage@dera.gov.uk, +44 1344 633203

  A UK patent has been granted: "Method for digitally processing
  images to determine the position of edges and/or corners therein for
  guidance of unmanned vehicle", UK Patent 2272285. Proprietor:
  Secretary of State for Defence, UK. 15 January 1997

  This code is issued for research purposes only and remains the
  property of the UK Secretary of State for Defence. This code must
  not be passed on without this header information being kept
  intact. This code must not be sold.

*/

/* }}} */
/* {{{ Readme First */

/**********************************************************************\

  This 3D version derived from 2d version SUSAN Version 2i

  SUSAN = Smallest Univalue Segment Assimilating Nucleus

  Email:    steve@fmrib.ox.ac.uk
  WWW:      http://www.fmrib.ox.ac.uk/~steve

  Related paper:
  @article{Smith97,
        author = "Smith, S.M. and Brady, J.M.",
        title = "{SUSAN} - A New Approach to Low Level Image Processing",
        journal = "Int. Journal of Computer Vision",
        pages = "45--78",
        volume = "23",
        number = "1",
        month = "May",
        year = 1997}

  To be registered for automatic (bug) updates of SUSAN, send an email.

  Note: FDT is the image data type

  Note: edge and corner finding have been taken out of the 3D version
  of SUSAN - it is not thought that they would be of great value.

  See following section for different machine information. Please
  report any bugs (and fixes). There are a few optional changes that
  can be made in the "defines" section which follows shortly.

  This code is written using an emacs folding mode, making moving
  around the different sections very easy. This is why there are
  various marks within comments and why comments are indented.

  SPATIAL CONTROL: d

  In SUSAN smoothing d controls the size of the Gaussian mask; its
  default is 4.0. Increasing d gives more smoothing. In edge finding,
  a fixed flat mask is used, either 37 pixels arranged in a "circle"
  (default), or a 3 by 3 mask which gives finer detail. In corner
  finding, only the larger 37 pixel mask is used; d is not
  variable. In smoothing, the flat 3 by 3 mask can be used instead of
  a larger Gaussian mask; this gives low smoothing and fast operation.

  BRIGHTNESS CONTROL: t

  In all three algorithms, t can be varied (default=20); this is the
  main threshold to be varied. It determines the maximum difference in
  greylevels between two pixels which allows them to be considered
  part of the same "region" in the image. Thus it can be reduced to
  give more edges or corners, i.e. to be more sensitive, and vice
  versa. In smoothing, reducing t gives less smoothing, and vice
  versa. Set t=10 for the test image available from the SUSAN web
  page.

  ITERATIONS:

  With SUSAN smoothing, more smoothing can also be obtained by
  iterating the algorithm several times. This has a different effect
  from varying d or t.

  BRIGHTNESS FUNCTION LUT IMPLEMENTATION:
  (Only read this if you are interested in the C implementation)

  The SUSAN brightness function is implemented as a LUT
  (Look-Up-Table) for speed. The resulting pointer-based code is a
  little hard to follow, so here is a brief explanation. In
  setup_brightness_lut() the LUT is setup. This mallocs enough space
  for *bp and then repositions the pointer to the centre of the
  malloced space. The SUSAN function e^-(x^6) or e^-(x^2) is
  calculated and converted to a unsigned char in the range 0-100, for all
  possible image brightness differences (including negative
  ones). Thus bp[23] is the output for a brightness difference of 23
  greylevels. In the SUSAN algorithms this LUT is used as follows:

  p=in + (i-3)*x_size + j - 1;
  p points to the first image pixel in the circular mask surrounding
  point (x,y).

  cp=bp + in[i*x_size+j];
  cp points to a position in the LUT corresponding to the brightness
  of the centre pixel (x,y).

  now for every pixel within the mask surrounding (x,y),
  n+=*(cp-*p++);
  the brightness difference function is found by moving the cp pointer
  down by an amount equal to the value of the pixel pointed to by p,
  thus subtracting the two brightness values and performing the
  exponential function. This value is added to n, the running USAN
  area.

  in SUSAN smoothing, the variable height mask is implemented by
  multiplying the above by the moving mask pointer, reset for each new
  centre pixel.
  tmp = *dpt++ * *(cp-brightness);

\**********************************************************************/

/* }}} */
/* {{{ defines, includes and typedefs */

/* ********** Optional settings */

typedef float TOTAL_TYPE;

/* ********** Leave the rest - but you may need to remove one or both of sys/file.h and malloc.h lines */

#include "libss/libss.h"
#include "libss/libavw.h"

void usage();
void setup_brightness_lut(unsigned char**, int, int);
FDT susan_median(FDT*, int, int, int);
void senlarge(FDT**, FDT*, int*, int*, int, int);
void susan_smoothing(int, FDT*, float, float, int, int, int, unsigned char*, int);
void susan_smoothing_usan(int, FDT*, FDT*, int*, float, float, int, int, int, unsigned char*, int);
void susan_smoothing_usan2(int, FDT*, FDT*, FDT*, int*, float, float, int, int, int, unsigned char*, unsigned char*, int);
FDT susan_median3d(FDT*, int, int, int, int, int);
void senlarge3d(FDT**, FDT *, int*, int*, int*, int);
void susan_smoothing3d(int, FDT*, float, float, float, int, int, int, unsigned char*, int);
void susan_smoothing3d_usan(int, FDT*, FDT*, int*, float, float, float, int, int, int, unsigned char*, int);
void susan_smoothing3d_usan2(int, FDT*, FDT*, FDT*, int*, float, float, float, int, int, int, unsigned char*, unsigned char*, int);

/* }}} */
/* {{{ usage() */

void usage()
{
  printf("\nUsage: susan_smooth <input> <bt> <output> <dt> <dim> <use_median> <n_usans> [<usan1> <bt1> [<usan2> <bt2>] <usanoutput>]\n\n");
  printf("<bt> is brightness threshold and should be greater than noise level and less than contrast of edges to be preserved.\n");
  printf("<dt> is spatial size (sigma) of smoothing, in mm.\n");
  printf("<dim> is dimensionality (2 or 3), depending on whether smoothing is to be within-plane (2) or fully 3D (3).\n");
  printf("<use_median> determines whether to use a local median filter in the cases where single-point noise is detected (0 or 1).\n");
  printf("<n_usans> determines whether the smoothing area (USAN) is to be found from secondary images (0, 1 or 2).\n\n");
  exit(0);
}

/* }}} */
/* {{{ setup_brightness_lut(bp,thresh,form) */

void setup_brightness_lut(unsigned char **bp, int thresh, int form)
{
  int   k, bp_size;
  float temp;

  if (thresh<=0) thresh=1; 

  bp_size = (int) pow( ((float)2) , ((float)(8*sizeof(FDT))) );
  /*fprintf(stderr,"bp_size=%d\n",bp_size);*/

  /* {{{ COMMENT test for total's type */

#ifdef FoldingComment

  /* {{{ test for total's type */

  tmp=0.1;
  if (tmp==0)
    fprintf(stderr,"FDT is non-float\n");
  else
  {
    fprintf(stderr,"FDT is float\n");
    bp_size=100;
  }

/* }}} */

#endif

/* }}} */

  *bp=(unsigned char *)malloc(bp_size*2+1);
  *bp=*bp+bp_size;

  for(k=-bp_size;k<=bp_size;k++)
  {
    temp=((float)k)/((float)thresh);
    temp=temp*temp;
    if (form==6)
      temp=temp*temp*temp;
    temp=100.0*exp(-temp);
    *(*bp+k)= (unsigned char)temp;
  }
}

/* }}} */
/* {{{ smoothing */

/* {{{ susan_median */

FDT susan_median(FDT *in, int i, int j, int x_size)
{
  FDT p[8],tmp;
  int k,l;

  p[0]=in[(i-1)*x_size+j-1];
  p[1]=in[(i-1)*x_size+j  ];
  p[2]=in[(i-1)*x_size+j+1];
  p[3]=in[(i  )*x_size+j-1];
  p[4]=in[(i  )*x_size+j+1];
  p[5]=in[(i+1)*x_size+j-1];
  p[6]=in[(i+1)*x_size+j  ];
  p[7]=in[(i+1)*x_size+j+1];

  for(k=0; k<7; k++)
    for(l=0; l<(7-k); l++)
      if (p[l]>p[l+1])
      {
        tmp=p[l]; p[l]=p[l+1]; p[l+1]=tmp;
      }

  return( (FDT) ( (((float)p[3])+((float)p[4])) / 2.0 ) );
}

/* }}} */
/* {{{ senlarge */

/* this enlarges "in" so that borders can be dealt with easily */

void senlarge(FDT **in, FDT *tmp_image, int *x_size, int *y_size, int z, int border)
{
int   i, j;

  for(i=0; i<*y_size; i++)   /* copy *in into tmp_image */
    memcpy(tmp_image+(i+border)*(*x_size+2*border)+border, *in+i* *x_size, *x_size *sizeof(FDT));

  for(i=0; i<border; i++) /* copy top and bottom rows; invert as many as necessary */
  {
    memcpy(tmp_image+(border-1-i)*(*x_size+2*border)+border,*in+i* *x_size,*x_size*sizeof(FDT));
    memcpy(tmp_image+(*y_size+border+i)*(*x_size+2*border)+border,*in+(*y_size-i-1)* *x_size,*x_size*sizeof(FDT));
  }

  for(i=0; i<border; i++) /* copy left and right columns */
    for(j=0; j<*y_size+2*border; j++)
    {
      tmp_image[j*(*x_size+2*border)+border-1-i]=tmp_image[j*(*x_size+2*border)+border+i];
      tmp_image[j*(*x_size+2*border)+ *x_size+border+i]=tmp_image[j*(*x_size+2*border)+ *x_size+border-1-i];
    }

  *x_size+=2*border;  /* alter image size */
  *y_size+=2*border;
  *in=tmp_image;      /* repoint in */

  /* {{{ COMMENT output tmp_image to file for debugging */

#ifdef FoldingComment

  /* {{{ output tmp_image to file for debugging */

  {
  FILE *ofp;
  char ff[100];

    sprintf(ff,"/tmp/Enlarge%d.raw",z);
    ofp=fopen(ff,"wb");
    fwrite(tmp_image,*x_size * *y_size*sizeof(FDT),1,ofp);
    fprintf(stderr,"%d %d %d",*x_size,*y_size,z);
  }

/* }}} */

#endif

/* }}} */
}

/* }}} */
/* {{{ void susan_smoothing */

void susan_smoothing(int three_by_three, FDT *in, float x_dt, float y_dt, int x_size, int y_size, int z, unsigned char *bp, int use_median)
{
/* {{{ vars */

int   increment, mask_size, x_mask_size=1, y_mask_size=1,
      i,j,x,y,area,tmp;
FDT   *ip, *tmp_image, *out=in, centre,brightness;
unsigned char *dp, *dpt, *cp;
TOTAL_TYPE total;

/* }}} */

  /* {{{ setup larger image and border sizes */

  if (three_by_three==0)
  {
    x_mask_size = ((int)(2.0 * x_dt)) + 1;
    y_mask_size = ((int)(2.0 * y_dt)) + 1;
  }
  else
    y_mask_size = x_mask_size = 1;

  mask_size=MAX(x_mask_size,y_mask_size);

  tmp_image = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * sizeof(FDT) );
  senlarge(&in,tmp_image,&x_size,&y_size,z,mask_size);

/* }}} */

  if (three_by_three==0)
  {     /* large Gaussian masks */
    /* {{{ setup distance lut */

  increment = x_size - ((x_mask_size*2)+1);

  dp     = (unsigned char *)malloc( ((x_mask_size*2)+1) * ((y_mask_size*2)+1) );
  dpt    = dp;
  x_dt   = -(2*x_dt*x_dt);
  y_dt   = -(2*y_dt*y_dt);

  for(i=-y_mask_size; i<=y_mask_size; i++)
  {
    for(j=-x_mask_size; j<=x_mask_size; j++)
    {
      x = (int) (100.0 * exp( ((float)(i*i))/y_dt + ((float)(j*j))/x_dt ));
      /*fprintf(stderr,"%d ",x);*/
      *dpt++ = (unsigned char)x;
    }
    /*fprintf(stderr,"\n");*/
  }

/* }}} */
    /* {{{ main section */

  for (i=mask_size;i<y_size-mask_size;i++) /* use mask_size as senlarge was isotropic */
  {
    for (j=mask_size;j<x_size-mask_size;j++)
    {
      area = 0;
      total = 0;
      dpt = dp;
      ip = in + ((i-y_mask_size)*x_size) + j - x_mask_size;
      centre = in[i*x_size+j];
      cp = bp + centre;
      for(y=-y_mask_size; y<=y_mask_size; y++)
      {
        for(x=-x_mask_size; x<=x_mask_size; x++)
	{
          brightness = *ip++;
          tmp = *dpt++ * *(cp-brightness);
          area += tmp;
          total += tmp * brightness;
        }
        ip += increment;
      }

      if (use_median)
	{
	  tmp = area-10000;
	  if (tmp<10000)
	    *out++=susan_median(in,i,j,x_size);
	  else
	    *out++=((total-(centre*10000))/tmp);
	}
      else
	*out++=total/area;
    }
  }

/* }}} */
  }
  else
  {     /* 3x3 constant mask */
    /* {{{ main section */

  for (i=1;i<y_size-1;i++)
  {
    for (j=1;j<x_size-1;j++)
    {
      area = 0;
      total = 0;
      ip = in + ((i-1)*x_size) + j - 1;
      centre = in[i*x_size+j];
      cp = bp + centre;

      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size-2;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size-2;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;

      if (use_median)
	{
	  tmp = area-100;
	  if (tmp<100)
	    *out++=susan_median(in,i,j,x_size);
	  else
	    *out++=(total-(centre*100))/tmp;
	}
      else
	*out++=total/area;
    }
  }

/* }}} */
  }
}

/* }}} */
/* {{{ void susan_smoothing_usan */

void susan_smoothing_usan(int three_by_three, FDT *in, FDT *usan, int *usansize, float x_dt, float y_dt, int x_size, int y_size, int z, unsigned char *bp, int use_median)
{
/* {{{ vars */

int   increment, mask_size, x_mask_size=1, y_mask_size=1,
      i,j,x,y;
FDT   *ip, *up, *tmp_image, *tmp_usan, *out=in, centre, ucentre, brightness, ubrightness;
unsigned char *dp, *dpt, *cp;
float area, tmp, total;

/* }}} */

  /* {{{ setup larger image and border sizes */

  if (three_by_three==0)
  {
    x_mask_size = ((int)(2.0 * x_dt)) + 1;
    y_mask_size = ((int)(2.0 * y_dt)) + 1;
  }
  else
    y_mask_size = x_mask_size = 1;

  mask_size=MAX(x_mask_size,y_mask_size);

  tmp_image = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * sizeof(FDT) );
  senlarge(&in,tmp_image,&x_size,&y_size,z,mask_size);

  x_size-=2*mask_size;  /* reset image size */
  y_size-=2*mask_size;

  tmp_usan = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * sizeof(FDT) );
  senlarge(&usan,tmp_usan,&x_size,&y_size,z,mask_size);

/* }}} */

  if (three_by_three==0)
  {     /* large Gaussian masks */
    /* {{{ setup distance lut */

  increment = x_size - ((x_mask_size*2)+1);

  dp     = (unsigned char *)malloc( ((x_mask_size*2)+1) * ((y_mask_size*2)+1) );
  dpt    = dp;
  x_dt   = -(2*x_dt*x_dt);
  y_dt   = -(2*y_dt*y_dt);

  for(i=-y_mask_size; i<=y_mask_size; i++)
  {
    for(j=-x_mask_size; j<=x_mask_size; j++)
    {
      x = (int) (100.0 * exp( ((float)(i*i))/y_dt + ((float)(j*j))/x_dt ));
      /*fprintf(stderr,"%d ",x);*/
      *dpt++ = (unsigned char)x;
    }
    /*fprintf(stderr,"\n");*/
  }

/* }}} */
    /* {{{ main section */

  for (i=mask_size;i<y_size-mask_size;i++) /* use mask_size as senlarge was isotropic */
  {
    for (j=mask_size;j<x_size-mask_size;j++)
    {
      area = 0;
      total = 0;
      dpt = dp;
      ip = in + ((i-y_mask_size)*x_size) + j - x_mask_size;
      up = usan + ((i-y_mask_size)*x_size) + j - x_mask_size;
      centre = in[i*x_size+j];
      ucentre = usan[i*x_size+j];
      cp = bp + ucentre;
      for(y=-y_mask_size; y<=y_mask_size; y++)
      {
        for(x=-x_mask_size; x<=x_mask_size; x++)
	{
          brightness = *ip++;
          ubrightness = *up++;
          tmp = *dpt++ * *(cp-ubrightness);
          area += tmp;
          total += tmp * brightness;
        }
        ip += increment;
        up += increment;
      }

      *usansize++ = area;

      if (use_median)
	{
	  tmp = area-10000.0;
	  if (tmp<10000.0)
	    *out++=susan_median(in,i,j,x_size);
	  else
	    *out++=((total-(centre*10000.0))/tmp);
	}
      else
	*out++=total/area;
    }
  }

/* }}} */
  }
  else
  {     /* 3x3 constant mask */
    /* {{{ main section */

  for (i=1;i<y_size-1;i++)
  {
    for (j=1;j<x_size-1;j++)
    {
      area = 0;
      total = 0;
      ip = in + ((i-1)*x_size) + j - 1;
      up = usan + ((i-1)*x_size) + j - 1;
      centre = in[i*x_size+j];
      ucentre = usan[i*x_size+j];
      cp = bp + ucentre;

      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip; ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size-2; up += x_size-2;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip; ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size-2; up += x_size-2;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip; ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;

      *usansize++ = area;

      if (use_median)
	{
	  tmp = area-100;
	  if (tmp<100)
	    *out++=susan_median(in,i,j,x_size);
	  else
	    *out++=(total-(centre*100))/tmp;
	}
      else
	*out++=total/area;
    }
  }

/* }}} */
  }
}

/* }}} */
/* {{{ void susan_smoothing_usan2 */

void susan_smoothing_usan2(int three_by_three, FDT *in, FDT *usan1, FDT *usan2, int *usansize, float x_dt, float y_dt, int x_size, int y_size, int z, unsigned char *bp1, unsigned char *bp2, int use_median)
{
/* {{{ vars */

int   increment, mask_size, x_mask_size=1, y_mask_size=1,
      i,j,x,y;
FDT   *ip, *up1, *up2, *tmp_image, *tmp_usan1, *tmp_usan2, *out=in, centre, ucentre1, ucentre2, brightness, ubrightness1, ubrightness2;
unsigned char *dp, *dpt, *cp1, *cp2;
TOTAL_TYPE area, total, tmp;

/* }}} */

  /* {{{ setup larger image and border sizes */

  if (three_by_three==0)
  {
    x_mask_size = ((int)(2.0 * x_dt)) + 1;
    y_mask_size = ((int)(2.0 * y_dt)) + 1;
  }
  else
    y_mask_size = x_mask_size = 1;

  mask_size=MAX(x_mask_size,y_mask_size);

  tmp_image = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * sizeof(FDT) );
  senlarge(&in,tmp_image,&x_size,&y_size,z,mask_size);

  x_size-=2*mask_size;  /* reset image size */
  y_size-=2*mask_size;

  tmp_usan1 = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * sizeof(FDT) );
  senlarge(&usan1,tmp_usan1,&x_size,&y_size,z,mask_size);

  x_size-=2*mask_size;  /* reset image size */
  y_size-=2*mask_size;

  tmp_usan2 = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * sizeof(FDT) );
  senlarge(&usan2,tmp_usan2,&x_size,&y_size,z,mask_size);

/* }}} */

  if (three_by_three==0)
  {     /* large Gaussian masks */
    /* {{{ setup distance lut */

  increment = x_size - ((x_mask_size*2)+1);

  dp     = (unsigned char *)malloc( ((x_mask_size*2)+1) * ((y_mask_size*2)+1) );
  dpt    = dp;
  x_dt   = -(2*x_dt*x_dt);
  y_dt   = -(2*y_dt*y_dt);

  for(i=-y_mask_size; i<=y_mask_size; i++)
  {
    for(j=-x_mask_size; j<=x_mask_size; j++)
    {
      x = (int) (100.0 * exp( ((float)(i*i))/y_dt + ((float)(j*j))/x_dt ));
      /*fprintf(stderr,"%d ",x);*/
      *dpt++ = (unsigned char)x;
    }
    /*fprintf(stderr,"\n");*/
  }

/* }}} */
    /* {{{ main section */

  for (i=mask_size;i<y_size-mask_size;i++) /* use mask_size as enlarge was isotropic */
  {
    for (j=mask_size;j<x_size-mask_size;j++)
    {
      area = 0;
      total = 0;
      dpt = dp;
      ip = in + ((i-y_mask_size)*x_size) + j - x_mask_size;
      up1 = usan1 + ((i-y_mask_size)*x_size) + j - x_mask_size;
      up2 = usan2 + ((i-y_mask_size)*x_size) + j - x_mask_size;
      centre = in[i*x_size+j];
      ucentre1 = usan1[i*x_size+j];
      ucentre2 = usan2[i*x_size+j];
      cp1 = bp1 + ucentre1;
      cp2 = bp2 + ucentre2;
      for(y=-y_mask_size; y<=y_mask_size; y++)
      {
        for(x=-x_mask_size; x<=x_mask_size; x++)
	{
          brightness = *ip++;
          ubrightness1 = *up1++;
          ubrightness2 = *up2++;
          tmp = (int)*dpt++ * (int)*(cp1-ubrightness1) * (int)*(cp2-ubrightness2);
          area += tmp;
          total += tmp * brightness;
        }
        ip += increment;
        up1 += increment;
        up2 += increment;
      }

      *usansize++ = area;

      if (use_median)
	{
	  tmp = area-1000000.0;
	  if (tmp<1000000.0)
	    *out++=susan_median(in,i,j,x_size);
	  else
	    *out++=((total-(centre*1000000.0))/tmp);
	}
      else
	*out++=total/area;
    }
  }

/* }}} */
  }
  else
  {     /* 3x3 constant mask */
    /* {{{ main section */

  for (i=1;i<y_size-1;i++)
  {
    for (j=1;j<x_size-1;j++)
    {
      area = 0;
      total = 0;
      ip = in + ((i-1)*x_size) + j - 1;
      up1 = usan1 + ((i-1)*x_size) + j - 1;
      up2 = usan2 + ((i-1)*x_size) + j - 1;
      centre = in[i*x_size+j];
      ucentre1 = usan1[i*x_size+j];
      ucentre2 = usan2[i*x_size+j];
      cp1 = bp1 + ucentre1;
      cp2 = bp2 + ucentre2;

      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;
      ip += x_size-2; up1 += x_size-2; up2 += x_size-2;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;
      ip += x_size-2; up1 += x_size-2; up2 += x_size-2;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area += tmp; total += tmp * brightness;

      *usansize++ = area;

      if (use_median)
	{
	  tmp = area-10000.0;
	  if (tmp<10000.0)
	    *out++=susan_median(in,i,j,x_size);
	  else
	    *out++=(total-(centre*10000.0))/tmp;
	}
      else
	*out++=total/area;
    }
  }

/* }}} */
  }
}

/* }}} */

/* {{{ susan_median3d */

FDT susan_median3d(FDT *in, int i, int j, int k, int x_size, int y_size)
{
FDT p[26],tmp;
int kk,ll;

  p[0] =in[(k-1)*x_size*y_size+(i-1)*x_size+j-1];
  p[1] =in[(k-1)*x_size*y_size+(i-1)*x_size+j  ];
  p[2] =in[(k-1)*x_size*y_size+(i-1)*x_size+j+1];
  p[3] =in[(k-1)*x_size*y_size+(i  )*x_size+j-1];
  p[4] =in[(k-1)*x_size*y_size+(i  )*x_size+j  ];
  p[5] =in[(k-1)*x_size*y_size+(i  )*x_size+j+1];
  p[6] =in[(k-1)*x_size*y_size+(i+1)*x_size+j-1];
  p[7] =in[(k-1)*x_size*y_size+(i+1)*x_size+j  ];
  p[8] =in[(k-1)*x_size*y_size+(i+1)*x_size+j+1];

  p[9] =in[k*x_size*y_size+(i-1)*x_size+j-1];
  p[10]=in[k*x_size*y_size+(i-1)*x_size+j  ];
  p[11]=in[k*x_size*y_size+(i-1)*x_size+j+1];
  p[12]=in[k*x_size*y_size+(i  )*x_size+j-1];
  p[13]=in[k*x_size*y_size+(i  )*x_size+j+1];
  p[14]=in[k*x_size*y_size+(i+1)*x_size+j-1];
  p[15]=in[k*x_size*y_size+(i+1)*x_size+j  ];
  p[16]=in[k*x_size*y_size+(i+1)*x_size+j+1];

  p[17]=in[(k+1)*x_size*y_size+(i-1)*x_size+j-1];
  p[18]=in[(k+1)*x_size*y_size+(i-1)*x_size+j  ];
  p[19]=in[(k+1)*x_size*y_size+(i-1)*x_size+j+1];
  p[20]=in[(k+1)*x_size*y_size+(i  )*x_size+j-1];
  p[21]=in[(k+1)*x_size*y_size+(i  )*x_size+j  ];
  p[22]=in[(k+1)*x_size*y_size+(i  )*x_size+j+1];
  p[23]=in[(k+1)*x_size*y_size+(i+1)*x_size+j-1];
  p[24]=in[(k+1)*x_size*y_size+(i+1)*x_size+j  ];
  p[25]=in[(k+1)*x_size*y_size+(i+1)*x_size+j+1];

  for(kk=0; kk<25; kk++)
    for(ll=0; ll<(25-kk); ll++)
      if (p[ll]>p[ll+1])
      {
        tmp=p[ll]; p[ll]=p[ll+1]; p[ll+1]=tmp;
      }

  return( (FDT)(((double)p[12]+(double)p[13]) / 2.0 ));
}

/* }}} */
/* {{{ senlarge3d */

/* this enlarges "in" so that borders can be dealt with easily */

void senlarge3d(FDT **in, FDT *tmp_image, int *x_size, int *y_size, int *z_size, int border)
{
int   i, j, k;

  for(k=0; k<*z_size; k++)   /* copy *in into tmp_image */
    for(i=0; i<*y_size; i++)
      memcpy(tmp_image+(k+border)*(*x_size+2*border)*(*y_size+2*border)+(i+border)*(*x_size+2*border)+border,
             *in+k* *x_size* *y_size+i* *x_size,   *x_size*sizeof(FDT));

  for(k=0; k<*z_size; k++)
  {
    for(i=0; i<border; i++) /* copy top and bottom rows; invert as many as necessary */
    {
      memcpy(tmp_image+(k+border)*(*x_size+2*border)*(*y_size+2*border)+(border-1-i)*(*x_size+2*border)+border,
             *in+ k* *y_size * *x_size + i* *x_size,   *x_size*sizeof(FDT));
      memcpy(tmp_image+(k+border)*(*x_size+2*border)*(*y_size+2*border)+(*y_size+border+i)*(*x_size+2*border)+border,
             *in+ k* *y_size * *x_size + (*y_size-i-1)* *x_size,   *x_size*sizeof(FDT));
    }
    for(i=0; i<border; i++) /* copy left and right columns */
      for(j=0; j<*y_size+2*border; j++)
      {
        tmp_image[(k+border)*(*x_size+2*border)*(*y_size+2*border)+j*(*x_size+2*border)+border-1-i]=
	  tmp_image[(k+border)*(*x_size+2*border)*(*y_size+2*border)+j*(*x_size+2*border)+border+i];
        tmp_image[(k+border)*(*x_size+2*border)*(*y_size+2*border)+j*(*x_size+2*border)+ *x_size+border+i]=
	  tmp_image[(k+border)*(*x_size+2*border)*(*y_size+2*border)+j*(*x_size+2*border)+ *x_size+border-1-i];
      }
  }

  for(k=0; k<border; k++)
  {
    memcpy(tmp_image+(border-1-k)*(*x_size+2*border)*(*y_size+2*border),
	   tmp_image+(border+k)*(*x_size+2*border)*(*y_size+2*border),  (*x_size+2*border)*(*y_size+2*border)*sizeof(FDT));
    memcpy(tmp_image+(*z_size+border+k)*(*x_size+2*border)*(*y_size+2*border),
	   tmp_image+(*z_size+border-1-k)*(*x_size+2*border)*(*y_size+2*border),  (*x_size+2*border)*(*y_size+2*border)*sizeof(FDT));
  }

  *x_size+=2*border;  /* alter image size */
  *y_size+=2*border;
  *z_size+=2*border;
  *in=tmp_image;      /* repoint in */

  /* {{{ COMMENT output tmp_image to file for debugging */

#ifdef FoldingComment

  {
    FILE *ofp;
    char filename[100];
    
    sprintf(filename,"/tmp/Enlarge.raw%d",in);
    ofp=fopen(filename,"wb");
    fwrite(tmp_image,*x_size * *y_size * *z_size*sizeof(FDT),1,ofp);
    fprintf(stderr,"%d %d %d",*x_size,*y_size,*z_size);
  }

#endif

/* }}} */
}

/* }}} */
/* {{{ void susan_smoothing3d(three_by_three,in,x_dt,y_dt,x_size,y_size,bp) */

void susan_smoothing3d(int three_by_three, FDT *in, float x_dt, float y_dt, float z_dt, int x_size, int y_size, int z_size, unsigned char *bp, int use_median)
{
/* {{{ vars */

int   increment, mask_size, x_mask_size=1, y_mask_size=1, z_mask_size=1,
      i,j,k,x,y,z,area,tmp;
FDT   *ip, *tmp_image, *out=in, centre, brightness;
unsigned char *dp, *dpt, *cp;
TOTAL_TYPE total;

/* }}} */

  /* {{{ setup larger image and border sizes */

  if (three_by_three==0)
  {
    x_mask_size = ((int)(2.0 * x_dt)) + 1;
    y_mask_size = ((int)(2.0 * y_dt)) + 1;
    z_mask_size = ((int)(2.0 * z_dt)) + 1;
  }
  else
    z_mask_size = y_mask_size = x_mask_size = 1;

  mask_size=MAX(MAX(x_mask_size,y_mask_size),z_mask_size);

  tmp_image = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * (z_size+mask_size*2)  * sizeof(FDT) );
  senlarge3d(&in,tmp_image,&x_size,&y_size,&z_size,mask_size);

/* }}} */

  if (three_by_three==0)
  {     /* large Gaussian masks */
    /* {{{ setup distance lut */

  increment = x_size - (x_mask_size*2+1);

  dp     = (unsigned char *)malloc( ((x_mask_size*2)+1) * ((y_mask_size*2)+1) * ((z_mask_size*2)+1) );
  dpt    = dp;
  x_dt   = -(2*x_dt*x_dt);
  y_dt   = -(2*y_dt*y_dt);
  z_dt   = -(2*z_dt*z_dt);

for(k=-z_mask_size; k<=z_mask_size; k++)
  for(i=-y_mask_size; i<=y_mask_size; i++)
  {
    for(j=-x_mask_size; j<=x_mask_size; j++)
    {
      x = (int) (100.0 * exp( ((float)(i*i))/y_dt + ((float)(j*j))/x_dt + ((float)(k*k))/z_dt ));
      /*fprintf(stderr,"%.4d ",x);*/
      *dpt++ = (unsigned char)x;
    }
    /*fprintf(stderr,"\n");*/
  }

/* }}} */
    /* {{{ main section */

for (k=mask_size;k<z_size-mask_size;k++) /* use mask_size as senlarge was isotropic */
{
  for (i=mask_size;i<y_size-mask_size;i++)
  {
    for (j=mask_size;j<x_size-mask_size;j++)
    {
      area = 0;
      total = 0;
      dpt = dp;
      ip = in + (k-z_mask_size)*x_size*y_size + (i-y_mask_size)*x_size + j-x_mask_size;
      centre = in[k*x_size*y_size+i*x_size+j];
      cp = bp + centre;
      for(z=-z_mask_size; z<=z_mask_size; z++)
      {
        for(y=-y_mask_size; y<=y_mask_size; y++)
        {
          for(x=-x_mask_size; x<=x_mask_size; x++)
	  {
            brightness = *ip++;
            tmp = *dpt++ * *(cp-brightness);
            area += tmp;
            total += tmp * brightness;
          }
          ip += increment;
        }
        ip += x_size*(y_size-(y_mask_size*2+1));
      }

      if (use_median)
	{
	  tmp = area-10000;
	  if (tmp<20000)
	    *out++=susan_median3d(in,i,j,k,x_size,y_size);
	  else
	    *out++=((total-(centre*10000))/tmp);
	}
      else
	*out++=total/area;
    }
  }
}

/* }}} */
  }
  else
  {     /* 3x3 constant mask */
    /* {{{ main section */

for (k=1;k<z_size-1;k++)
{
  for (i=1;i<y_size-1;i++)
  {
    for (j=1;j<x_size-1;j++)
    {
      area = 0;
      total = 0;
      ip = in + (k-1)*x_size*y_size + (i-1)*x_size + j-1;
      centre = in[k*x_size*y_size+i*x_size+j];
      cp = bp + centre;

      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size-2;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size-2;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size*(y_size-2) - 2;

      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size-2;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size-2;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size*(y_size-2) - 2;

      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size-2;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      ip += x_size-2;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;
      brightness=*ip; tmp=*(cp-brightness); area += tmp; total += tmp * brightness;

      if (use_median)
	{
          tmp = area-100;
	  if (tmp<200)
	    *out++=susan_median3d(in,i,j,k,x_size,y_size);
	  else
	    *out++=(total-(centre*100))/tmp;
	}
      else
	*out++=total/area;
    }
  }
}

/* }}} */
  }
}

/* }}} */
/* {{{ void susan_smoothing3d_usan(three_by_three,in,usan,x_dt,y_dt,x_size,y_size,bp) */

void susan_smoothing3d_usan(int three_by_three, FDT *in, FDT *usan, int *usansize, float x_dt, float y_dt, float z_dt, int x_size, int y_size, int z_size, unsigned char *bp, int use_median)
{
/* {{{ vars */

int   increment, mask_size, x_mask_size=1, y_mask_size=1, z_mask_size=1,
      i,j,k,x,y,z;
FDT   *ip, *up, *tmp_image, *tmp_usan, *out=in, centre, ucentre, brightness, ubrightness;
unsigned char *dp, *dpt, *cp;
float area, tmp, total;

/* }}} */

  /* {{{ setup larger image and border sizes */

  if (three_by_three==0)
  {
    x_mask_size = ((int)(2.0 * x_dt)) + 1;
    y_mask_size = ((int)(2.0 * y_dt)) + 1;
    z_mask_size = ((int)(2.0 * z_dt)) + 1;
  }
  else
    z_mask_size = y_mask_size = x_mask_size = 1;

  mask_size=MAX(MAX(x_mask_size,y_mask_size),z_mask_size);

  tmp_image = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * (z_size+mask_size*2)  * sizeof(FDT) );
  senlarge3d(&in,tmp_image,&x_size,&y_size,&z_size,mask_size);

  x_size-=2*mask_size;  /* reset image size */
  y_size-=2*mask_size;
  z_size-=2*mask_size;

  tmp_usan = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * (z_size+mask_size*2)  * sizeof(FDT) );
  senlarge3d(&usan,tmp_usan,&x_size,&y_size,&z_size,mask_size);

/* }}} */

  if (three_by_three==0)
  {     /* large Gaussian masks */
    /* {{{ setup distance lut */

  increment = x_size - (x_mask_size*2+1);

  dp     = (unsigned char *)malloc( ((x_mask_size*2)+1) * ((y_mask_size*2)+1) * ((z_mask_size*2)+1) );
  dpt    = dp;
  x_dt   = -(2*x_dt*x_dt);
  y_dt   = -(2*y_dt*y_dt);
  z_dt   = -(2*z_dt*z_dt);

for(k=-z_mask_size; k<=z_mask_size; k++)
  for(i=-y_mask_size; i<=y_mask_size; i++)
  {
    for(j=-x_mask_size; j<=x_mask_size; j++)
    {
      x = (int) (100.0 * exp( ((float)(i*i))/y_dt + ((float)(j*j))/x_dt + ((float)(k*k))/z_dt ));
      /*fprintf(stderr,"%.4d ",x);*/
      *dpt++ = (unsigned char)x;
    }
    /*fprintf(stderr,"\n");*/
  }

/* }}} */
    /* {{{ main section */

for (k=mask_size;k<z_size-mask_size;k++) /* use mask_size as senlarge was isotropic */
{
  for (i=mask_size;i<y_size-mask_size;i++)
  {
    for (j=mask_size;j<x_size-mask_size;j++)
    {
      area = 0;
      total = 0;
      dpt = dp;
      ip = in + (k-z_mask_size)*x_size*y_size + (i-y_mask_size)*x_size + j-x_mask_size;
      up = usan + (k-z_mask_size)*x_size*y_size + (i-y_mask_size)*x_size + j-x_mask_size;
      centre = in[k*x_size*y_size+i*x_size+j];
      ucentre = usan[k*x_size*y_size+i*x_size+j];
      cp = bp + ucentre;
      for(z=-z_mask_size; z<=z_mask_size; z++)
      {
        for(y=-y_mask_size; y<=y_mask_size; y++)
        {
          for(x=-x_mask_size; x<=x_mask_size; x++)
	  {
            brightness = *ip++;
            ubrightness = *up++;
            tmp = *dpt++ * *(cp-ubrightness);
            area += tmp;
            total += tmp * brightness;
          }
          ip += increment;
          up += increment;
        }
        ip += x_size*(y_size-(y_mask_size*2+1));
        up += x_size*(y_size-(y_mask_size*2+1));
      }

      *usansize++ = area;

      if (use_median)
	{
	  tmp = area-10000.0;
	  if (tmp<20000.0)
	    *out++=susan_median3d(in,i,j,k,x_size,y_size);
	  else
	    *out++=((total-(centre*10000.0))/tmp);
	}
      else
	*out++=total/area;
    }
  }
}

/* }}} */
  }
  else
  {     /* 3x3 constant mask */
    /* {{{ main section */

for (k=1;k<z_size-1;k++)
{
  for (i=1;i<y_size-1;i++)
  {
    for (j=1;j<x_size-1;j++)
    {
      area = 0;
      total = 0;
      ip = in + (k-1)*x_size*y_size + (i-1)*x_size + j-1;
      up = usan + (k-1)*x_size*y_size + (i-1)*x_size + j-1;
      centre = in[k*x_size*y_size+i*x_size+j];
      ucentre = usan[k*x_size*y_size+i*x_size+j];
      cp = bp + ucentre;

      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size-2; up += x_size-2;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size-2; up += x_size-2;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size*(y_size-2) - 2; up += x_size*(y_size-2) - 2;

      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size-2; up += x_size-2;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size-2; up += x_size-2;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size*(y_size-2) - 2; up += x_size*(y_size-2) - 2;

      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size-2; up += x_size-2;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      ip += x_size-2; up += x_size-2;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip++; ubrightness=*up++; tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;
      brightness=*ip;   ubrightness=*up;   tmp=*(cp-ubrightness); area += tmp; total += tmp * brightness;

      *usansize++ = area;

      if (use_median)
	{
          tmp = area-100.0;
	  if (tmp<200)
	    *out++=susan_median3d(in,i,j,k,x_size,y_size);
	  else
	    *out++=(total-(centre*100.0))/tmp;
	}
      else
	*out++=total/area;
    }
  }
}

/* }}} */
  }
}

/* }}} */
/* {{{ void susan_smoothing3d_usan2 */

void susan_smoothing3d_usan2(int three_by_three, FDT *in, FDT *usan1, FDT *usan2, int *usansize, float x_dt, float y_dt, float z_dt, int x_size, int y_size, int z_size, unsigned char *bp1, unsigned char *bp2, int use_median)
{
/* {{{ vars */

int   increment, mask_size, x_mask_size=1, y_mask_size=1, z_mask_size=1,
      i,j,k,x,y,z;
FDT   *ip, *up1, *up2, *tmp_image, *tmp_usan1, *tmp_usan2, *out=in, centre, ucentre1, ucentre2, brightness, ubrightness1, ubrightness2;
unsigned char *dp, *dpt, *cp1, *cp2;
TOTAL_TYPE area, total, tmp;

/* }}} */

  /* {{{ setup larger image and border sizes */

  if (three_by_three==0)
  {
    x_mask_size = ((int)(2.0 * x_dt)) + 1;
    y_mask_size = ((int)(2.0 * y_dt)) + 1;
    z_mask_size = ((int)(2.0 * z_dt)) + 1;
  }
  else
    z_mask_size = y_mask_size = x_mask_size = 1;

  mask_size=MAX(MAX(x_mask_size,y_mask_size),z_mask_size);

  tmp_image = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * (z_size+mask_size*2)  * sizeof(FDT) );
  senlarge3d(&in,tmp_image,&x_size,&y_size,&z_size,mask_size);

  x_size-=2*mask_size;  /* reset image size */
  y_size-=2*mask_size;
  z_size-=2*mask_size;

  tmp_usan1 = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * (z_size+mask_size*2)  * sizeof(FDT) );
  senlarge3d(&usan1,tmp_usan1,&x_size,&y_size,&z_size,mask_size);

  x_size-=2*mask_size;  /* reset image size */
  y_size-=2*mask_size;
  z_size-=2*mask_size;

  tmp_usan2 = (FDT *) malloc( (x_size+mask_size*2) * (y_size+mask_size*2) * (z_size+mask_size*2)  * sizeof(FDT) );
  senlarge3d(&usan2,tmp_usan2,&x_size,&y_size,&z_size,mask_size);

/* }}} */

  if (three_by_three==0)
  {     /* large Gaussian masks */
    /* {{{ setup distance lut */

  increment = x_size - (x_mask_size*2+1);

  dp     = (unsigned char *)malloc( ((x_mask_size*2)+1) * ((y_mask_size*2)+1) * ((z_mask_size*2)+1) );
  dpt    = dp;
  x_dt   = -(2*x_dt*x_dt);
  y_dt   = -(2*y_dt*y_dt);
  z_dt   = -(2*z_dt*z_dt);

for(k=-z_mask_size; k<=z_mask_size; k++)
  for(i=-y_mask_size; i<=y_mask_size; i++)
  {
    for(j=-x_mask_size; j<=x_mask_size; j++)
    {
      x = (int) (100.0 * exp( ((float)(i*i))/y_dt + ((float)(j*j))/x_dt + ((float)(k*k))/z_dt ));
      /*fprintf(stderr,"%.4d ",x);*/
      *dpt++ = (unsigned char)x;
    }
    /*fprintf(stderr,"\n");*/
  }

/* }}} */
    /* {{{ main section */

for (k=mask_size;k<z_size-mask_size;k++) /* use mask_size as senlarge was isotropic */
{
  for (i=mask_size;i<y_size-mask_size;i++)
  {
    for (j=mask_size;j<x_size-mask_size;j++)
    {
      area = 0;
      total = 0;
      dpt = dp;
      ip = in + (k-z_mask_size)*x_size*y_size + (i-y_mask_size)*x_size + j-x_mask_size;
      up1 = usan1 + (k-z_mask_size)*x_size*y_size + (i-y_mask_size)*x_size + j-x_mask_size;
      up2 = usan2 + (k-z_mask_size)*x_size*y_size + (i-y_mask_size)*x_size + j-x_mask_size;
      centre = in[k*x_size*y_size+i*x_size+j];
      ucentre1 = usan1[k*x_size*y_size+i*x_size+j];
      ucentre2 = usan2[k*x_size*y_size+i*x_size+j];
      cp1 = bp1 + ucentre1;
      cp2 = bp2 + ucentre2;
      for(z=-z_mask_size; z<=z_mask_size; z++)
      {
        for(y=-y_mask_size; y<=y_mask_size; y++)
        {
          for(x=-x_mask_size; x<=x_mask_size; x++)
	  {
            brightness = *ip++;
            ubrightness1 = *up1++;
            ubrightness2 = *up2++;
            tmp = (int)*dpt++ * (int)*(cp1-ubrightness1) * (int)*(cp2-ubrightness2);
            area += tmp;
            total += tmp * brightness;
          }
          ip += increment;
          up1 += increment;
          up2 += increment;
        }
        ip += x_size*(y_size-(y_mask_size*2+1));
        up1 += x_size*(y_size-(y_mask_size*2+1));
        up2 += x_size*(y_size-(y_mask_size*2+1));
      }

      *usansize++ = area;

      if (use_median)
	{
	  tmp = area-1000000.0;
	  if (tmp<2000000.0)
	    *out++=susan_median3d(in,i,j,k,x_size,y_size);
	  else
	    *out++=((total-(centre*1000000.0))/tmp);
	}
      else
	*out++=total/area;
    }
  }
}

/* }}} */
  }
  else
  {     /* 3x3 constant mask */
    /* {{{ main section */

for (k=1;k<z_size-1;k++)
{
  for (i=1;i<y_size-1;i++)
  {
    for (j=1;j<x_size-1;j++)
    {
      area = 0;
      total = 0;
      ip = in + (k-1)*x_size*y_size + (i-1)*x_size + j-1;
      up1 = usan1 + (k-1)*x_size*y_size + (i-1)*x_size + j-1;
      up2 = usan2 + (k-1)*x_size*y_size + (i-1)*x_size + j-1;
      centre = in[k*x_size*y_size+i*x_size+j];
      ucentre1 = usan1[k*x_size*y_size+i*x_size+j];
      ucentre2 = usan2[k*x_size*y_size+i*x_size+j];
      cp1 = bp1 + ucentre1;
      cp2 = bp2 + ucentre2;

      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      ip += x_size-2; up1 += x_size-2; up2 += x_size-2;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      ip += x_size-2; up1 += x_size-2; up2 += x_size-2;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      ip += x_size*(y_size-2) - 2; up1 += x_size*(y_size-2) - 2; up2 += x_size*(y_size-2) - 2;

      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      ip += x_size-2; up1 += x_size-2; up2 += x_size-2;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      ip += x_size-2; up1 += x_size-2; up2 += x_size-2;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      ip += x_size*(y_size-2) - 2; up1 += x_size*(y_size-2) - 2; up2 += x_size*(y_size-2) - 2;

      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      ip += x_size-2; up1 += x_size-2; up2 += x_size-2;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      ip += x_size-2; up1 += x_size-2; up2 += x_size-2;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip++; ubrightness1=*up1++; ubrightness2=*up2++; tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;
      brightness=*ip;   ubrightness1=*up1;   ubrightness2=*up2;   tmp=*(cp1-ubrightness1) * *(cp2-ubrightness2); area+=tmp; total+=tmp * brightness;

      *usansize++ = area;

      if (use_median)
	{
          tmp = area-10000.0;
	  if (tmp<20000.0)
	    *out++=susan_median3d(in,i,j,k,x_size,y_size);
	  else
	    *out++=(total-(centre*10000.0))/tmp;
	}
      else
	*out++=total/area;
    }
  }
}

/* }}} */
  }
}

/* }}} */

/* }}} */
/* {{{ main(argc, argv) */

int main(int argc, char *argv[])
{
  /* {{{ vars */

char   tmpc[1000];
unsigned char  *bp[2];
float  dt, x_dt, y_dt, z_dt;
int    bt, i, z, dimension,
       three_by_three=0,
       use_median=1,
       usans, x_size, y_size, z_size, *usansize=NULL;
FDT    *InVolume, *UsanVolume[2], *UsanSizeVolume=NULL;
image_struct im, im2, im3, im4;

/* }}} */

  /* {{{ setup */

if (argc < 8)
{
  usage();
  return(1);
}

avw_read(argv[1],&im);
InVolume = im.i;
bt=atoi(argv[2]); /* if usans=0 bt is from input */

/* {{{ setup dimensions, dt, dimensionality and use-median */

x_size = im.x;
y_size = im.y;
z_size = im.z;

dt = atof(argv[4]);
x_dt = dt / im.xv;
y_dt = dt / im.yv;
z_dt = dt / im.zv;
if ((x_dt<0.01)||(y_dt<0.01)||(z_dt<0.01))
     three_by_three=1;

strncpy(tmpc,argv[5],1);
tmpc[1]=0;
dimension = atoi(tmpc);
if (z_size<2)
     dimension=2;

strncpy(tmpc,argv[6],1);
tmpc[1]=0;
use_median = atoi(tmpc);

/* }}} */

usans = atoi(argv[7]);

if (usans>0)
{
  avw_read(argv[8],&im2);
  UsanVolume[0] = im2.i;
  bt=atoi(argv[9]);  /* if usans=1, bt is from usan1 */
  avw_read(argv[8],&im4); /* dummy input - just sets up im4 really, we're not interested in the data */
  UsanSizeVolume = im4.i;
  usansize=(int*)malloc((sizeof(int))*x_size*y_size*z_size);
}

setup_brightness_lut(&bp[0],bt,2);

if (usans==2)
{
  avw_read(argv[10],&im3);
  UsanVolume[1] = im3.i;
  bt=atoi(argv[11]);
  setup_brightness_lut(&bp[1],bt,2);
}

/* }}} */

  if (usans==0)
    /* {{{ main stuff */

{
  if (dimension==2)
    for (z=0; z<z_size; z++)
      susan_smoothing(three_by_three,InVolume+z*x_size*y_size,x_dt,y_dt,x_size,y_size,z,bp[0],use_median);
  else
    susan_smoothing3d(three_by_three,InVolume,x_dt,y_dt,z_dt,x_size,y_size,z_size,bp[0],use_median);
}

/* }}} */
  if (usans==1)
    /* {{{ main stuff */

{
  if (dimension==2)
    for (z=0; z<z_size; z++)
      susan_smoothing_usan(three_by_three,InVolume+z*x_size*y_size,UsanVolume[0]+z*x_size*y_size,
			   usansize+z*x_size*y_size,x_dt,y_dt,x_size,y_size,z,bp[0],use_median);
  else
    susan_smoothing3d_usan(three_by_three,InVolume,UsanVolume[0],
			   usansize,x_dt,y_dt,z_dt,x_size,y_size,z_size,bp[0],use_median);
}

/* }}} */
  if (usans==2)
    /* {{{ main stuff */

{
  if (dimension==2)
    for (z=0; z<z_size; z++)
      susan_smoothing_usan2(three_by_three,InVolume+z*x_size*y_size,UsanVolume[0]+z*x_size*y_size,UsanVolume[1]+z*x_size*y_size,
			    usansize+z*x_size*y_size,x_dt,y_dt,x_size,y_size,z,bp[0],bp[1],use_median);
  else
    susan_smoothing3d_usan2(three_by_three,InVolume,UsanVolume[0],UsanVolume[1],
			    usansize,x_dt,y_dt,z_dt,x_size,y_size,z_size,bp[0],bp[1],use_median);
}

/* }}} */

  avw_write(argv[3],im);

  /* {{{ rescale and write usansize */

/* note that this bit pretty much relies on FDT>=16SI as the output may be too big for 8UI/8SI */

if (usans>0)
{
  char filename[1000];
  int divide=10000;

  if (usans==2)
    divide*=100;

  if (three_by_three)
    divide/=100;
  
  for(i=0;i<x_size*y_size*z_size;i++)
    UsanSizeVolume[i]=(usansize[i]/divide);

  sprintf(filename,"%s_usan_size",argv[3]);
  avw_write(filename,im4);
}

/* }}} */

  return 0;
}

/* }}} */
