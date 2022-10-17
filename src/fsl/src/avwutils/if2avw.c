/* {{{ Copyright etc. */

/*  if2avw - convert from interfile to AVW

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

/* }}} */
/* {{{ defines, includes and typedefs */

#include "libss/libssbase.h"
#include "image/avwimage.h"
#include <assert.h>

void usage(void);
char *find_key(FILE *, char *, const char *);

/* }}} */
/* {{{ usage */

void usage(void)
{
  fprintf(stderr,"Usage: if2avw <header> <data> [-s (for SPET)]\n");
  exit(1);
}

/* }}} */
/* {{{ find_key */

char *find_key(FILE *fd, char *control_line, const char *key)
{
int  i;

  while (fgets(control_line, 1000, fd)!=NULL)
    if ( (strncmp(control_line,key,strlen(key))==0) ||
	 (strncmp(control_line+1,key,strlen(key))==0) )
      {
	for (i=0; (i<strlen(control_line)-3)&&(strncmp(":=",control_line+i,2)!=0); i++);
	if (i<strlen(control_line)-3)
	  {
	    rewind(fd);
	    if (strncmp(" ",control_line+i+2,1)==0) i++; /* gobble space */
	    if (control_line[strlen(control_line)-1]==10) control_line[strlen(control_line)-1]=0; /* gobble RTN */
	    return(control_line+i+2);
	  }
      }
	  
  printf("Error: key \"%s\" not found.\n",key);
  exit(1);
}

/* }}} */
/* {{{ main */

int main(int argc, char **argv)
{
  /* {{{ vars */

  FILE *fd;
  AVWImage *im=NewAVW();
  char control_line[1000], outfile[1000], dob[1000], dos[1000], name[1000], tmps[1000];
  int x, y, z, t, spet=0;

/* }}} */

  if (argc<3)
    usage();

  if (argc>3)
    spet=1;

  /* {{{ read header and setup AVW header */

  fd=fopen(argv[1],"rb");

  im->xdim=atoi(find_key(fd,control_line,"matrix size [1]"));
  im->ydim=atoi(find_key(fd,control_line,"matrix size [2]"));
  im->zdim=atoi(find_key(fd,control_line,"total number of images"));
  im->vdim=1;

  im->xvox=atof(find_key(fd,control_line,"scaling factor (mm/pixel) [1]"));
  im->yvox=atof(find_key(fd,control_line,"scaling factor (mm/pixel) [2]"));
  im->zvox=atof(find_key(fd,control_line,"centre-centre slice separation (pixels)"));
  im->zvox*=im->xvox;

  if (spet==0)
     im->zvox=2; /* kludge because headers got screwed up somehow */

  im->dt=4;

  if (spet==0)
  {  
    strncpy(im->study,find_key(fd,control_line,"original institution"),18); im->study[17]=0;
  }
  else
    im->study[0]=0;

  strncpy(im->descrip,find_key(fd,control_line,"patient dob"),80);        im->descrip[79]=0;
  strncpy(im->scannum,find_key(fd,control_line,"study ID"),10);           im->scannum[9]=0;
  strncpy(im->date,find_key(fd,control_line,"study date"),10);            im->date[9]=0;
  strncpy(im->patient,find_key(fd,control_line,"patient name"),10);       im->patient[9]=0;

  strncpy(name,find_key(fd,control_line,"patient name"),10); name[1]=0; /* only use 1 character */
  strncpy(dob,find_key(fd,control_line,"patient dob"),100);
  strncpy(dos,find_key(fd,control_line,"study date"),100);

  for(x=0, y=0; x<=strlen(dob); x++) /* remove ":" from dob */
    if(dob[x]!=':')
    {
      dob[y]=dob[x];
      y++;
    }

  sprintf(outfile,"%s_%s.%s",name,dob,dos);

  fclose(fd);

/* }}} */
  /* {{{ process image data */

  fd=fopen(argv[2],"rb");
  im->data=malloc(2*im->xdim*im->ydim*im->zdim*im->vdim);
  fread(im->data,2,im->xdim*im->ydim*im->zdim*im->vdim,fd);
  fclose(fd);

  for(t=0;t<im->vdim;t++)
    for(z=0;z<im->zdim/2;z++)
      for(y=0;y<im->ydim;y++)
	for(x=0;x<im->xdim;x++)
	  {
	   signed short tmp=
	     ((int)((unsigned short*)im->data)[t*im->zdim*im->ydim*im->xdim+z*im->ydim*im->xdim+y*im->xdim+x])-SHRT_MAX-1;
	    ((signed short*)im->data)[t*im->zdim*im->ydim*im->xdim+z*im->ydim*im->xdim+y*im->xdim+x]=
	      ((int)((unsigned short*)im->data)[t*im->zdim*im->ydim*im->xdim+(im->zdim-1-z)*im->ydim*im->xdim+y*im->xdim+x])-SHRT_MAX-1;
	    ((signed short*)im->data)[t*im->zdim*im->ydim*im->xdim+(im->zdim-1-z)*im->ydim*im->xdim+y*im->xdim+x]=
	      tmp;
	  }

  if (im->zdim%2==1) /* do middle slice if odd number of slices */
    {
      z=im->zdim/2+1;
      for(t=0;t<im->vdim;t++)
	for(y=0;y<im->ydim;y++)
	  for(x=0;x<im->xdim;x++)
	    ((signed short*)im->data)[t*im->zdim*im->ydim*im->xdim+z*im->ydim*im->xdim+y*im->xdim+x]=
	      ((int)((unsigned short*)im->data)[t*im->zdim*im->ydim*im->xdim+z*im->ydim*im->xdim+y*im->xdim+x])-SHRT_MAX-1;
    }

/* }}} */
  /* {{{ test if already exists, otherwise write AVW image */

sprintf(tmps,"%s.hdr",outfile);
if((fd=fopen(tmps,"rb"))==NULL)
{
  assert(im->header->hk.sizeof_hdr == 348);
  WriteAVW(outfile,im);
  return 0;
}
else
{
  printf("Problem - %s already exists\n",outfile);
  return 1;
}

/* }}} */
}

/* }}} */
