/* {{{ Copyright */

/*  libpicio - collection of 2D image io routines

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
/* {{{ includes */

#include "libss.h"

/* }}} */
/* {{{ write_pgm */

int write_pgm ( char *filename, int x_size, int y_size,	unsigned char *i )
{
  FILE *ofp;
  int  x, y;

  if ((ofp=fopen(filename,"wb"))==NULL)
    {
      printf("Can't open %s for writing\n",filename);
      return(1);
    }

  fprintf(ofp,"P5\n");
  fprintf(ofp,"%d %d\n",x_size,y_size);
  fprintf(ofp,"255\n");

  for(y=0; y<y_size; y++)
    for(x=0; x<x_size; x++)
      fwrite(&i[y*x_size+x],1,1,ofp);

  fclose(ofp);
  return(0);
}

/* }}} */
/* {{{ write_ppm */

int write_ppm ( char *filename, int x_size, int y_size,
		unsigned char *r, unsigned char *g, unsigned char *b )
{
  FILE *ofp;
  int  x, y;

  if ((ofp=fopen(filename,"wb"))==NULL)
    {
      printf("Can't open %s for writing\n",filename);
      return(1);
    }

  fprintf(ofp,"P6\n");
  fprintf(ofp,"%d %d\n",x_size,y_size);
  fprintf(ofp,"255\n");

  for(y=0; y<y_size; y++)
    for(x=0; x<x_size; x++)
      {
	fwrite(&r[y*x_size+x],1,1,ofp);
	fwrite(&g[y*x_size+x],1,1,ofp);
	fwrite(&b[y*x_size+x],1,1,ofp);
      }

  fclose(ofp);
  return(0);
}

/* }}} */
