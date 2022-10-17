/*  costfns.cc

    Mark Jenkinson, FMRIB Image Analysis Group

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

// Interpolation functions
//  Written by Mark Jenkinson  18/2/99

#include <iostream>
#include <cassert>
#include <values.h>
#include "newimage/costfns.h"
#include "miscmaths/miscmaths.h"
#include "miscimfns.h"
#include "interpolation.h"
#include "newmatio.h"
#include "Globaloptions.h"

#ifndef NO_NAMESPACE
 using namespace MISCIMFNS;
 using namespace MISCMATHS;
 using namespace INTERPOLATION;

 namespace COSTFNS {
#endif

   float sinc(float x)
     {
       if (fabs(x)<1e-7) { return 1.0-fabs(x); }
       float y=M_PI*x;
       return sin(y)/y;
     }

   float hanning(float x, int w)
     {
       if (fabs(x)>w) 
	 return 0.0;
       else
	 return (0.5 + 0.5 *cos(M_PI*x/w));
     }

   void setupkernel()
     {
       Globaloptions::getInstance().kernelwidth=3;
       // set x between +/- kernelwidth
       for (int n=0; n<=200; n++) {
	 float x=(n-100)/100.0*Globaloptions::getInstance().kernelwidth;
	 Globaloptions::getInstance().sinckernel[n] = 
	   sinc(x)*hanning(x,Globaloptions::getInstance().kernelwidth);
       }
     }

   float kernelval(float x, int w)
     {
       // effectively returns  sinc(x)*hanning(x,w);
       if (fabs(x)>w) return 0.0;
       float dn = x/w*100.0 + 100;
       int n = (int) floor(dn);
       dn -= n;
       if (n>=200) return 0.0;
       if (n<0) return 0.0;

       return Globaloptions::getInstance().sinckernel[n]*(1.0-dn) 
	        + Globaloptions::getInstance().sinckernel[n+1]*dn;
     }


   float sinc_interpolation(const volume& v, const float x, const float y,
			    const float z)
     {
          // kernel half-width  (i.e. range is +/- w)
       int w=Globaloptions::getInstance().kernelwidth;  
       if (w<1) { 
	 setupkernel(); 
	 w=Globaloptions::getInstance().kernelwidth;
       }
       
       int ix0, iy0, iz0;
       ix0 = (int) floor(x);
       iy0 = (int) floor(y);
       iz0 = (int) floor(z);

       float convsum=0.0, interpval=0.0, kersum=0.0;
       float *sincz, *sincy, *sincx;
       sincz = new float[2*w+1];
       sincy = new float[2*w+1];
       sincx = new float[2*w+1];

       for (int d=-w; d<=w; d++) {
	 sincz[d+w] = kernelval((z-iz0+d),w);
	 sincy[d+w] = kernelval((y-iy0+d),w);
	 sincx[d+w] = kernelval((x-ix0+d),w);
       }

       int xj, yj, zj;
       for (int z1=iz0-w; z1<=iz0+w; z1++) {
	 zj=iz0-z1+w;
	 for (int y1=iy0-w; y1<=iy0+w; y1++) {
	   yj=iy0-y1+w;
	   for (int x1=ix0-w; x1<=ix0+w; x1++) {
	     if (v.in_bounds(x1,y1,z1)) {
	       xj=ix0-x1+w;
	       float sincfac = sincx[xj] * sincy[yj] * sincz[zj];
	       convsum += v(x1,y1,z1) * sincfac;
	       kersum += sincfac;
	     } 
	   }
	 }
       }

       delete [] sincx; delete [] sincy; delete [] sincz;

       if (fabs(kersum)>1e-9) {
	 interpval = convsum / kersum;
       } else {
	 return v.backgroundval();
       }
       return interpval;

     }
   


   void findrangex(unsigned int &xmin1 , unsigned int &xmax1,
		   float o1, float o2, float o3,
		   float a11, float a21, float a31,
		   unsigned int xb1, unsigned int yb1, unsigned int zb1,
		   float xb2, float yb2, float zb2) {
     
     float x1, x2, xmin, xmax, xmin0, xmax0;
     
     xmin0 = 0;
     xmax0 = xb1;
      
     if (fabs(a11)<1.0e-8) {
       if ((0.0<=o1) && (o1<=xb2)) {
	 x1 = -1.0e8; x2 = 1.0e8;
       } else {
	 x1 = -1.0e8; x2 = -1.0e8;
       }
     } else {
       x1 = -o1/a11;
       x2 = (xb2-o1)/a11;
     }
     xmin = Min(x1,x2);
     xmax = Max(x1,x2);
     // intersect ranges
     xmin0 = Max(xmin0,xmin);
     xmax0 = Min(xmax0,xmax);
	  
     if (fabs(a21)<1.0e-8) {
       if ((0.0<=o2) && (o2<=yb2)) {
	 x1 = -1.0e8; x2 = 1.0e8;
       } else {
	 x1 = -1.0e8; x2 = -1.0e8;
       }
     } else {
       x1 = -o2/a21;
       x2 = (yb2-o2)/a21;
     }
     xmin = Min(x1,x2);
     xmax = Max(x1,x2);
     // intersect ranges
     xmin0 = Max(xmin0,xmin);
     xmax0 = Min(xmax0,xmax);

     if (fabs(a31)<1.0e-8) {
       if ((0.0<=o3) && (o3<=zb2)) {
	 x1 = -1.0e8; x2 = 1.0e8;
       } else {
	 x1 = -1.0e8; x2 = -1.0e8;
       }
     } else {
       x1 = -o3/a31;
       x2 = (zb2-o3)/a31;
     }
     xmin = Min(x1,x2);
     xmax = Max(x1,x2);
     // intersect ranges
     xmin0 = Max(xmin0,xmin);
     xmax0 = Min(xmax0,xmax);
    
     //assert(xmin0>=0.0);
     //assert(xmax0<=xb1);

     if (xmax0<xmin0) {
       xmax1=0;
       xmin1=1;
     } else {
       xmin1 = (unsigned int) ceil(xmin0);
       xmax1 = (unsigned int) floor(xmax0);
     }

   }

   //--------------------------------------------------------------------//

   float corr_ratio_smoothed(const volume& vref, const volume& vtest,
		    int *bindex, const Matrix& aff,
		    const int no_bins, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      float *numy;
      numy = new float[no_bins+1];
      int b=0;
 
      for (int i=0; i<=no_bins; i++) {
	numy[i]=0.0; sumy[i]=0.0;  sumy2[i]=0.0;
      }

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4);
      float val,o1,o2,o3;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.getx();
      smoothy = smoothsize / vtest.gety();
      smoothz = smoothsize / vtest.getz();

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      unsigned int xmin, xmax;
      int *bptr;

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    b=*bptr;
	    weight=1.0;
	    if (o1<smoothx)  weight*=o1/smoothx;
	    else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
	    if (o2<smoothy)  weight*=o2/smoothy;
	    else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
	    if (o3<smoothz)  weight*=o3/smoothz;
	    else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
	    if (weight<0.0)  weight=0.0;
	    numy[b]+=weight;
	    sumy[b]+=weight*val;
	    sumy2[b]+=weight*val*val;

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }


      float corr_ratio=0.0, var=0.0, totsumy=0.0, totsumy2=0.0;
      float numtoty=0.0;

      // correct for occasion lapses into the last bin
      numy[no_bins-1] += numy[no_bins];
      sumy[no_bins-1] += sumy[no_bins];
      sumy2[no_bins-1] += sumy2[no_bins];
      numy[no_bins]=0.0;
      sumy[no_bins]=0.0;
      sumy2[no_bins]=0.0;

      // now calculate the individual variances for each iso-set
      //  weighting them by the number of pixels from Image x that contribute
      for (b=0; b<no_bins; b++) {
	if (numy[b]>2.0) {
	  numtoty += numy[b];
	  totsumy += sumy[b];
	  totsumy2 += sumy2[b];
	  // the following should be the variance of the bth iso-subset
	  var = (sumy2[b] - sumy[b]*sumy[b]/numy[b] ) / ( numy[b]-1);
	  // cerr << "Set #" << b << " has " << numy[b] << " elements and " 
	  //   << var << " variance" << endl;
	  corr_ratio += var * ((float) numy[b]);
	}
      }
      delete [] numy; delete [] sumy; delete [] sumy2;

      // normalise the weighting of numy[]
      if (numtoty>0)  corr_ratio/=((float) numtoty);
      // calculate the total variance of Image y and then normalise by this
      if (numtoty>1)
	var = ( totsumy2 - totsumy*totsumy/numtoty ) / (numtoty - 1);
      //cerr << "TOTALS are:" << endl 
      //   << " numerator variance is : " << corr_ratio << endl
      //   << " and denominator variance is: " << var << " from " << numtoty 
      //   << " valid elements" << endl;
      if (var>0.0)  corr_ratio/=var;
      // the above is actually 1 - correlation ratio, so correct this now
      if ( (numtoty<=1) || (var<=0.0) )
	return 0.0;   // the totally uncorrelated condition
      else
	return (1.0 - corr_ratio);

      // an alternative is to return 1.0/corr_ratio (=1/(1-correlation ratio))
      //  which may be better at rewarding gains near the best solution

      return 0;

    }

  ///////////////////////////////////////////////////////////////////////

   float corr_ratio(const volume& vref, const volume& vtest,
		    int *bindex, const Matrix& aff,
		    const int no_bins)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float *sumy, *sumy2;
      sumy = new float[no_bins+1];
      sumy2 = new float[no_bins+1];
      int *numy;
      numy = new int[no_bins+1];
      int b=0;
 
      for (int i=0; i<=no_bins; i++) {
	numy[i]=0; sumy[i]=0.0;  sumy2[i]=0.0;
      }

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4);
      float val,o1,o2,o3;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      unsigned int xmin, xmax;
      int *bptr;

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    b=*bptr;
	    numy[b]++;
	    sumy[b]+=val;
	    sumy2[b]+=val*val;

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }


      float corr_ratio=0.0, var=0.0, totsumy=0.0, totsumy2=0.0;
      int numtoty=0;

      // correct for occasion lapses into the last bin
      numy[no_bins-1] += numy[no_bins];
      sumy[no_bins-1] += sumy[no_bins];
      sumy2[no_bins-1] += sumy2[no_bins];
      numy[no_bins]=0;
      sumy[no_bins]=0.0;
      sumy2[no_bins]=0.0;

      // now calculate the individual variances for each iso-set
      //  weighting them by the number of pixels from Image x that contribute
      for (b=0; b<no_bins; b++) {
	if (numy[b]>2) {
	  numtoty += numy[b];
	  totsumy += sumy[b];
	  totsumy2 += sumy2[b];
	  // the following should be the variance of the bth iso-subset
	  var = (sumy2[b] - sumy[b]*sumy[b]/((float) numy[b]) ) /
	    ((float) (numy[b]-1));
	  // cerr << "Set #" << b << " has " << numy[b] << " elements and " 
	  //   << var << " variance" << endl;
	  corr_ratio += var * ((float) numy[b]);
	}
      }
      delete [] numy; delete [] sumy; delete [] sumy2;

      // normalise the weighting of numy[]
      if (numtoty>0)  corr_ratio/=((float) numtoty);
      // calculate the total variance of Image y and then normalise by this
      if (numtoty>1)
	var = ( totsumy2 - totsumy*totsumy/((float) numtoty) ) /
	  ((float) (numtoty - 1));
      //cerr << "TOTALS are:" << endl 
      //   << " numerator variance is : " << corr_ratio << endl
      //   << " and denominator variance is: " << var << " from " << numtoty 
      //   << " valid elements" << endl;
      if (var>0.0)  corr_ratio/=var;
      // the above is actually 1 - correlation ratio, so correct this now
      if ( (numtoty<=1) || (var<=0.0) )
	return 0.0;   // the totally uncorrelated condition
      else
	return (1.0 - corr_ratio);

      // an alternative is to return 1.0/corr_ratio (=1/(1-correlation ratio))
      //  which may be better at rewarding gains near the best solution

      return 0;

    }

  ///////////////////////////////////////////////////////////////////////


  float woods_fn(const volume& vref, const volume& vtest, int *bindex, 
		 const Matrix& aff, const int no_bins)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float *sum, *sum2;
      sum = new float[no_bins+1];
      sum2 = new float[no_bins+1];
      int *num;
      num = new int[no_bins+1];
      int b=0;

      for (int i=0; i<=no_bins; i++) {
	num[i]=0; sum[i]=0.0;  sum2[i]=0.0;
      }
  
      float val=0.0;
      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    b=*bptr;
	    num[b]++;
	    sum[b]+=val;
	    sum2[b]+=val*val;

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      // now calculate  W = sum_j (n_j/N)*(sigma_j / mu_j)
      //  where n_j = num[j], N = sum_j n_j, mu_j = sum[j]/num[j]
      //        sigma_j^2 = sum_j 1/(n_j - 1) * (sum2[j] - mu_j^2 * n_j)
      float woods=0.0, stdev=0.0, var=0.0;
      int numtot=0;
      for (b=0; b<=no_bins; b++) {
	if (num[b]>2) {
	  numtot += num[b];
	  // the following should be the variance of the bth subset
	  var = (sum2[b] - sum[b]*sum[b]/((float) num[b]) ) /
	    ((float) (num[b]-1));
	  if (var>0.0)
	    stdev = sqrt(var);
	  else
	    stdev = 0.0;
	  if (sum[b]>0)
	    woods += Sqr((float) num[b])*stdev/sum[b];
	  else
	    woods += Sqr((float) num[b])*stdev;
	}
      }
      delete [] num; delete [] sum; delete [] sum2;
      if (numtot>0) {
	woods/=((float) numtot);
	return woods;
      } else {
	return MAXFLOAT;
      }
    }


  ///////////////////////////////////////////////////////////////////////


  float woods_fn_smoothed(const volume& vref, const volume& vtest, int *bindex, 
			  const Matrix& aff, const int no_bins, 
			  const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float *sum, *sum2;
      sum = new float[no_bins+1];
      sum2 = new float[no_bins+1];
      float *num;
      num = new float[no_bins+1];
      int b=0;

      for (int i=0; i<=no_bins; i++) {
	num[i]=0.0; sum[i]=0.0;  sum2[i]=0.0;
      }
  
      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.getx();
      smoothy = smoothsize / vtest.gety();
      smoothz = smoothsize / vtest.getz();

      float val=0.0;
      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    b=*bptr;
	    weight=1.0;
	    if (o1<smoothx)  weight*=o1/smoothx;
	    else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
	    if (o2<smoothy)  weight*=o2/smoothy;
	    else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
	    if (o3<smoothz)  weight*=o3/smoothz;
	    else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
	    if (weight<0.0)  weight=0.0;
	    num[b]+=weight;
	    sum[b]+=weight*val;
	    sum2[b]+=weight*val*val;

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      // now calculate  W = sum_j (n_j/N)*(sigma_j / mu_j)
      //  where n_j = num[j], N = sum_j n_j, mu_j = sum[j]/num[j]
      //        sigma_j^2 = sum_j 1/(n_j - 1) * (sum2[j] - mu_j^2 * n_j)
      float woods=0.0, stdev=0.0, var=0.0;
      float numtot=0.0;
      for (b=0; b<=no_bins; b++) {
	if (num[b]>2.0) {
	  numtot += num[b];
	  // the following should be the variance of the bth subset
	  var = (sum2[b] - sum[b]*sum[b]/(num[b])) / (num[b]-1.0);
	  if (var>0.0)
	    stdev = sqrt(var);
	  else
	    stdev = 0.0;
	  if (sum[b]>0)
	    woods += Sqr(num[b])*stdev/sum[b];
	  else
	    woods += Sqr(num[b])*stdev;
	}
      }
      delete [] num; delete [] sum; delete [] sum2;
      if (numtot>0) {
	woods/=numtot;
	return woods;
      } else {
	return MAXFLOAT;
      }
    }


  ///////////////////////////////////////////////////////////////////////


  float normcorr(const volume& vref, const volume& vtest,
		 const Matrix& aff)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr=0.0;
      float sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
      float sumxA=0.0, sumyA=0.0, sumx2A=0.0, sumy2A=0.0, sumxyA=0.0;
      float sumxB=0.0, sumyB=0.0, sumx2B=0.0, sumy2B=0.0, sumxyB=0.0;
      float varx=0.0, vary=0.0, varxy=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      long int num=0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    num++;
	    valx = vref(x,y,z);
	    valy = val;
	    sumx += valx;
	    sumx2 += valx*valx;
	    sumy += valy;
	    sumy2 += valy*valy;
	    sumxy += valx*valy;

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumxA+=sumx; sumyA+=sumy; 
	  sumx2A+=sumx2; sumy2A+=sumy2; sumxyA+=sumxy; 
	  sumx=0.0; sumy=0.0; sumx2=0.0; sumy2=0.0; sumxy=0.0;
	}
	sumxB+=sumxA; sumyB+=sumyA; 
	sumx2B+=sumx2A; sumy2B+=sumy2A; sumxyB+=sumxyA; 
	sumxA=0.0; sumyA=0.0; sumx2A=0.0; sumy2A=0.0; sumxyA=0.0;
      }
      assert(fabs(sumxA+sumx)<1e-9);
      sumx=sumxB; sumy=sumyB; sumx2=sumx2B; sumy2=sumy2B; sumxy=sumxyB;

      corr = 0.0;  // uncorrelated (worst) case
      if (num>2) {
	float numsq = ((float) num)*((float) num);
	varxy = sumxy/((float) num-1) - (sumx*sumy)/numsq;
	varx = sumx2/((float) num-1) - (sumx*sumx)/numsq;
	vary = sumy2/((float) num-1) - (sumy*sumy)/numsq;
	if ((varx>0.0) && (vary>0.0)) {
	  corr = varxy/sqrt(varx)/sqrt(vary);
	} 
      }
      return corr;
    }

  
  ///////////////////////////////////////////////////////////////////////


  float normcorr_smoothed(const volume& vref, const volume& vtest,
			  const Matrix& aff, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr=0.0;
      float sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
      float sumxA=0.0, sumyA=0.0, sumx2A=0.0, sumy2A=0.0, sumxyA=0.0;
      float sumxB=0.0, sumyB=0.0, sumx2B=0.0, sumy2B=0.0, sumxyB=0.0;
      float varx=0.0, vary=0.0, varxy=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      float num=0.0, numA=0.0, numB=0.0;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.getx();
      smoothy = smoothsize / vtest.gety();
      smoothz = smoothsize / vtest.getz();

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    weight=1.0;
	    if (o1<smoothx)  weight*=o1/smoothx;
	    else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
	    if (o2<smoothy)  weight*=o2/smoothy;
	    else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
	    if (o3<smoothz)  weight*=o3/smoothz;
	    else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
	    if (weight<0.0)  weight=0.0;

	    valx = vref(x,y,z);
	    valy = val;
	    num += weight;
	    sumx += weight*valx;
	    sumx2 += weight*valx*valx;
	    sumy += weight*valy;
	    sumy2 += weight*valy*valy;
	    sumxy += weight*valx*valy;

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  numA+=num; sumxA+=sumx; sumyA+=sumy; 
	  sumx2A+=sumx2; sumy2A+=sumy2; sumxyA+=sumxy; 
	  sumx=0.0; sumy=0.0; sumx2=0.0; sumy2=0.0; sumxy=0.0;
	}
	numB+=numA; sumxB+=sumxA; sumyB+=sumyA; 
	sumx2B+=sumx2A; sumy2B+=sumy2A; sumxyB+=sumxyA; 
	sumxA=0.0; sumyA=0.0; sumx2A=0.0; sumy2A=0.0; sumxyA=0.0;
      }
      assert(fabs(sumxA+sumx)<1e-9);
      num = numB;
      sumx=sumxB; sumy=sumyB; sumx2=sumx2B; sumy2=sumy2B; sumxy=sumxyB;
  
      corr = 0.0;  // uncorrelated (worst) case
      if (num>2.0) {
	varxy = sumxy/(num-1.0) - (sumx*sumy)/(num*num);
	varx = sumx2/(num-1.0) - (sumx*sumx)/(num*num);
	vary = sumy2/(num-1.0) - (sumy*sumy)/(num*num);
	if ((varx>0.0) && (vary>0.0)) {
	  corr = varxy/sqrt(varx)/sqrt(vary);
	} 
      }
      return corr;
    }


  ///////////////////////////////////////////////////////////////////////
 
  float normcorr_smoothed_sinc(const volume& vref, const volume& vtest,
			       const Matrix& aff, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float corr=0.0;
      float sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
      float sumxA=0.0, sumyA=0.0, sumx2A=0.0, sumy2A=0.0, sumxyA=0.0;
      float sumxB=0.0, sumyB=0.0, sumx2B=0.0, sumy2B=0.0, sumxyB=0.0;
      float varx=0.0, vary=0.0, varxy=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      float num=0.0, numA=0.0, numB=0.0;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.getx();
      smoothy = smoothsize / vtest.gety();
      smoothz = smoothsize / vtest.getz();

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = sinc_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    weight=1.0;
	    if (o1<smoothx)  weight*=o1/smoothx;
	    else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
	    if (o2<smoothy)  weight*=o2/smoothy;
	    else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
	    if (o3<smoothz)  weight*=o3/smoothz;
	    else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
	    if (weight<0.0)  weight=0.0;

	    valx = vref(x,y,z);
	    valy = val;
	    num += weight;
	    sumx += weight*valx;
	    sumx2 += weight*valx*valx;
	    sumy += weight*valy;
	    sumy2 += weight*valy*valy;
	    sumxy += weight*valx*valy;

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  numA+=num; sumxA+=sumx; sumyA+=sumy; 
	  sumx2A+=sumx2; sumy2A+=sumy2; sumxyA+=sumxy; 
	  sumx=0.0; sumy=0.0; sumx2=0.0; sumy2=0.0; sumxy=0.0;
	}
	numB+=numA; sumxB+=sumxA; sumyB+=sumyA; 
	sumx2B+=sumx2A; sumy2B+=sumy2A; sumxyB+=sumxyA; 
	sumxA=0.0; sumyA=0.0; sumx2A=0.0; sumy2A=0.0; sumxyA=0.0;
      }
      assert(fabs(sumxA+sumx)<1e-9);
      num = numB;
      sumx=sumxB; sumy=sumyB; sumx2=sumx2B; sumy2=sumy2B; sumxy=sumxyB;
  
      corr = 0.0;  // uncorrelated (worst) case
      if (num>2.0) {
	varxy = sumxy/(num-1.0) - (sumx*sumy)/(num*num);
	varx = sumx2/(num-1.0) - (sumx*sumx)/(num*num);
	vary = sumy2/(num-1.0) - (sumy*sumy)/(num*num);
	if ((varx>0.0) && (vary>0.0)) {
	  corr = varxy/sqrt(varx)/sqrt(vary);
	} 
      }
      return corr;
    }

  
  ///////////////////////////////////////////////////////////////////////


  float leastsquares(const volume& vref, const volume& vtest,
		     const Matrix& aff)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float lsq=0.0;
      float sum=0.0, sumA=0.0, sumB=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      long int num=0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    num++;
	    valx = vref(x,y,z);
	    valy = val;
	    sum += (valx-valy)*(valx-valy);

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumA+=sum; sum=0.0;
	}
	sumB+=sumA; sumA=0.0;
      }
      assert(fabs(sumA+sum)<1e-9);
      sum = sumB;

      if (num>1) {
	lsq = sum/((float) num);
      } else {
	  // return the worst cost = (max-min)^2
	float vmin1, vmin2, vmax1, vmax2;
	get_min_max(vref,vmin1,vmax1);
	get_min_max(vtest,vmin2,vmax2);
	lsq = (Max(vmax2,vmax2)-Min(vmin1,vmin2));
	lsq = lsq*lsq;
      }
      
      return lsq;
    }


  ///////////////////////////////////////////////////////////////////////


  float leastsquares_smoothed(const volume& vref, const volume& vtest,
			      const Matrix& aff, const float smoothsize)
    {
      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1
      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      float smoothx, smoothy, smoothz, weight;
      smoothx = smoothsize / vtest.getx();
      smoothy = smoothsize / vtest.gety();
      smoothz = smoothsize / vtest.getz();

      float lsq=0.0;
      float sum=0.0, sumA=0.0, sumB=0.0;
      float valx=0.0, valy=0.0, val=0.0;
      float num=0.0, numA=0.0, numB=0.0;

      unsigned int xmin, xmax;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    weight=1.0;
	    if (o1<smoothx)  weight*=o1/smoothx;
	    else if ((xb2-o1)<smoothx) weight*=(xb2-o1)/smoothx;
	    if (o2<smoothy)  weight*=o2/smoothy;
	    else if ((yb2-o2)<smoothy) weight*=(yb2-o2)/smoothy;
	    if (o3<smoothz)  weight*=o3/smoothz;
	    else if ((zb2-o3)<smoothz) weight*=(zb2-o3)/smoothz;
	    if (weight<0.0)  weight=0.0;

	    valx = vref(x,y,z);
	    valy = val;
	    num+=weight;
	    sum += weight*(valx-valy)*(valx-valy);

	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	  sumA+=sum; sum=0.0;
	  numA+=num; num=0.0;
	}
	sumB+=sumA; sumA=0.0;
	numB+=numA; numA=0.0;
      }
      assert(fabs(sumA+sum)<1e-9);
      sum = sumB;  num = numB;
  
      
      if (num>1.0) {
	lsq = sum/num;
      } else {
	  // return the worst cost = (max-min)^2
	float vmin1, vmin2, vmax1, vmax2;
	get_min_max(vref,vmin1,vmax1);
	get_min_max(vtest,vmin2,vmax2);
	lsq = (Max(vmax2,vmax2)-Min(vmin1,vmin2));
	lsq = lsq*lsq;
      }
      
      return lsq;
    }


  ///////////////////////////////////////////////////////////////////////

  void calc_entropy(const volume& vref, const volume& vtest,
		    int *bindex,  const Matrix& aff,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2,
		    float& jointentropy, float& margentropy1,
		    float& margentropy2)
    {
      // the joint and marginal entropies between the two images are
      //  calculated here and returned
      // the last parameter, plnp, is a vector containing values of -p*log(p)
      //  which makes the calculation significantly more efficient

      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	  jointhist[i]=0;
      }
      for (int i=0; i<=no_bins; i++) {
	  marghist1[i]=0;
	  marghist2[i]=0;
      }

      long int a,b;
      float b1=no_bins/(maxtest-mintest), b0=-mintest*no_bins/(maxtest-mintest);
      float val=0.0;
 
      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    a=*bptr;
	    b=(long int) (val*b1 + b0);
	    if (b>=no_bins) b=no_bins-1;
	    if (b<0) b=0;
	    (jointhist[a*(no_bins+1) + b])++;
	    (marghist1[a])++;
	    (marghist2[b])++;

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      // note that the plnp values indexed by integers such that:
      //   plnp(n) = n/N * log(n/N)
      float p=0.0;
      int n=0, psize=plnp.Nrows();
      long int nvoxels = (long int) (vref.rows() * vref.columns() * 
				     vref.slices());
      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	n = jointhist[i];
	if (n>0) {
	  if (n<=psize)
	    jointentropy+=plnp(n);
	  else {
	    p = ((float) n) / ((float) nvoxels);
	    jointentropy+= - p*log(p);
	  }
	}
      }
      for (int i=0; i<=no_bins; i++) {
	n = marghist1[i];
	if (n>0) {
	  if (n<=psize)
	    margentropy1+=plnp(n);
	  else {
	    p = ((float) n) / ((float) nvoxels);
	    margentropy1+= - p*log(p);
	  }
	}
      }
      long int noverlap=0;
      for (int i=0; i<=no_bins; i++) {
	n = marghist2[i];
	if (n>0) {
	  noverlap += n;
	  if (n<=psize)
	    margentropy2+=plnp(n);
	  else {
	    //cerr << ":";
	    p = ((float) n) / ((float) nvoxels);
	    margentropy2+= - p*log(p);
	  }
	}
      }

      // correct for difference in total histogram size
      //  that is: noverlap vs nvoxels
      // H_1 = N_0/N_1 * H_0 + log(N_1/N_0)
      //     = N_0/N_1 * H_0 - log(N_0/N_1)
      if (noverlap > 0) {
	float nratio = ((float) nvoxels) / ((float) noverlap);
	jointentropy = nratio * jointentropy - log(nratio);
	margentropy1 = nratio * margentropy1 - log(nratio);
	margentropy2 = nratio * margentropy2 - log(nratio);
      } else {
	// Put in maximum entropy values as base cases = BAD registration
	jointentropy = 2.0*log(no_bins);
	margentropy1 = log(no_bins);
	margentropy2 = log(no_bins);
      }
      return;
    }



  float mutual_info(const volume& vref, const volume& vtest,
		    int *bindex, const Matrix& aff,
		    const float mintest, const float maxtest,
		    const int no_bins, const ColumnVector& plnp, 
		    int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      float mutualinformation = margentropy1 + margentropy2 - jointentropy;
      return mutualinformation;
    }



  float normalised_mutual_info(const volume& vref, const volume& vtest,
			       int *bindex, const Matrix& aff,
			       const float mintest, const float maxtest,
			       const int no_bins, const ColumnVector& plnp, 
			       int *jointhist, int *marghist1, int *marghist2)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
		   plnp,jointhist,marghist1,marghist2,
		   jointentropy,margentropy1,margentropy2);
      float normmi;
      if (fabs(jointentropy)<1e-9) {
	normmi = 0.0;  // BAD registration result
      } else {
	normmi = (margentropy1 + margentropy2)/jointentropy;
      }
      return normmi;
    }


  ///////////////////////////////////////////////////////////////////////

   void calc_smoothed_entropy(const volume& vref, const volume& vtest,
			      int *bindex,  const Matrix& aff,
			      const float mintest, const float maxtest,
			      const int no_bins,
			      float *jointhist, float *marghist1, 
			      float *marghist2,
			      float& jointentropy, float& margentropy1,
			      float& margentropy2,
			      const float smoothsize, const float fuzzyfrac)
   {
      // the joint and marginal entropies between the two images are
      //  calculated here and returned

      // Do everything in practice via the inverse transformation
      // That is, for every point in vref, calculate the pre-image in
      //  vtest to which it corresponds, and interpolate vtest to get the
      //  value there.
      // Also, the sampling transformations must be accounted for:
      //     T_vox1->vox2 = (T_samp2)^-1 * T_world * T_samp1

      Matrix iaffbig = vtest.sampling_matrix().i() * aff.i() *
	                     vref.sampling_matrix();  
      Matrix iaff=iaffbig.SubMatrix(1,3,1,3);
      unsigned int xb1=vref.columns()-1, yb1=vref.rows()-1, zb1=vref.slices()-1;
      float  xb2 = ((float) vtest.columns())-1.0001,
	yb2=((float) vtest.rows())-1.0001, zb2=((float) vtest.slices())-1.0001;

      float a11=iaff(1,1), a12=iaff(1,2), a13=iaff(1,3), a14=iaffbig(1,4),
	a21=iaff(2,1), a22=iaff(2,2), a23=iaff(2,3), a24=iaffbig(2,4),
	a31=iaff(3,1), a32=iaff(3,2), a33=iaff(3,3), a34=iaffbig(3,4), o1,o2,o3;

      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	  jointhist[i]=0;
      }
      for (int i=0; i<=no_bins; i++) {
	  marghist1[i]=0;
	  marghist2[i]=0;
      }

      long int a;
      float b1=no_bins/(maxtest-mintest), b0=-mintest*no_bins/(maxtest-mintest);
      float val=0.0;
 
      float smoothx, smoothy, smoothz, geomweight, wcentre, wplus, wminus, bidx;
      long int bcentre, bplus, bminus;
      smoothx = smoothsize / vtest.getx();
      smoothy = smoothsize / vtest.gety();
      smoothz = smoothsize / vtest.getz();

      unsigned int xmin, xmax;
      int *bptr;

      // The matrix algebra below has been hand-optimized from
      //  [o1 o2 o3] = a * [x y z]  at each iteration

      for (unsigned int z=0; z<=zb1; z++) { 
	for (unsigned int y=0; y<=yb1; y++) { 

	  o1= y*a12 + z*a13 + a14;  // x=0
	  o2= y*a22 + z*a23 + a24;  // x=0
	  o3= y*a32 + z*a33 + a34;  // x=0
	
	  // determine range
	  findrangex(xmin,xmax,o1,o2,o3,a11,a21,a31,xb1,yb1,zb1,xb2,yb2,zb2);

	  o1 += xmin * a11;
	  o2 += xmin * a21;
	  o3 += xmin * a31;

	  bptr = get_bindexptr(xmin,y,z,vref,bindex);

	  for (unsigned int x=xmin; x<=xmax; x++) {

  	    val = q_tri_interpolation(vtest,o1,o2,o3);
	    
	    // do the cost function record keeping...
	    geomweight=1.0;
	    if (o1<smoothx)  geomweight*=o1/smoothx;
	    else if ((xb2-o1)<smoothx) geomweight*=(xb2-o1)/smoothx;
	    if (o2<smoothy)  geomweight*=o2/smoothy;
	    else if ((yb2-o2)<smoothy) geomweight*=(yb2-o2)/smoothy;
	    if (o3<smoothz)  geomweight*=o3/smoothz;
	    else if ((zb2-o3)<smoothz) geomweight*=(zb2-o3)/smoothz;
	    if (geomweight<0.0)  geomweight=0.0;

	    // do the cost function record keeping...
	    a=*bptr;
	    bidx=val*b1 + b0;
	    bcentre=(long int) (bidx);
	    bplus = bcentre + 1;
	    bminus = bcentre - 1;
	    if (bcentre>=no_bins) {
	      bcentre=no_bins-1;
	      bplus = bcentre;
	    }
	    if (bcentre<0) {
	      bcentre=0;
	      bminus = bcentre;
	    }
	    if (bplus>=no_bins) { bplus = no_bins-1; }
	    if (bminus<0) { bminus = 0; }
	    // Fuzzy binning weights
	    bidx = fabs(bidx - (int) bidx);  // get fractional component : [0,1]
	    if (bidx<fuzzyfrac) {
	      wcentre = 0.5 + 0.5*(bidx/fuzzyfrac);
	      wminus = 1 - wcentre;
	      wplus = 0;
	    } else if (bidx>(1.0-fuzzyfrac)) {
	      wcentre = 0.5 + 0.5*((1.0-bidx)/fuzzyfrac);
	      wplus = 1 - wcentre;
	      wminus=0;
	    } else {
	      wcentre = 1;
	      wplus =0;
	      wminus=0;
	    }
	    (jointhist[a*(no_bins+1) + bcentre])+=geomweight*wcentre;
	    (marghist2[bcentre])+=geomweight*wcentre;
	    (jointhist[a*(no_bins+1) + bplus])+=geomweight*wplus;
	    (marghist2[bplus])+=geomweight*wplus;
	    (jointhist[a*(no_bins+1) + bminus])+=geomweight*wminus;
	    (marghist2[bminus])+=geomweight*wminus;
	    (marghist1[a])+=geomweight;

	    bptr++;
	    o1 += a11;
	    o2 += a21;
	    o3 += a31;
	  }
	}
      }

      float p=0.0, n=0.0;
      long int nvoxels = (long int) (vref.rows() * vref.columns() * 
				     vref.slices());
      for (long int i=0; i<((no_bins+1)*(no_bins+1)); i++) {
	n = jointhist[i];
	if (n>0) {
	  p = n / ((float) nvoxels);
	  jointentropy+= - p*log(p);
	}
      }
      for (int i=0; i<=no_bins; i++) {
	n = marghist1[i];
	if (n>0) {
	  p = n / ((float) nvoxels);
	  margentropy1+= - p*log(p);
	}
      }
      float noverlap=0;
      for (int i=0; i<=no_bins; i++) {
	n = marghist2[i];
	if (n>0) {
	  noverlap += n;
	  p = n / ((float) nvoxels);
	  margentropy2+= - p*log(p);
	}
      }

      // correct for difference in total histogram size
      //  that is: noverlap vs nvoxels
      // H_1 = N_0/N_1 * H_0 + log(N_1/N_0)
      //     = N_0/N_1 * H_0 - log(N_0/N_1)
      if (noverlap > 0) {
	float nratio = ((float) nvoxels) / ((float) noverlap);
	jointentropy = nratio * jointentropy - log(nratio);
	margentropy1 = nratio * margentropy1 - log(nratio);
	margentropy2 = nratio * margentropy2 - log(nratio);
      } else {
	// Put in maximum entropy values as base cases = BAD registration
	jointentropy = 2.0*log(no_bins);
	margentropy1 = log(no_bins);
	margentropy2 = log(no_bins);
      }
      return;
    }



  float mutual_info_smoothed(const volume& vref, const volume& vtest,
			     int *bindex, const Matrix& aff,
			     const float mintest, const float maxtest,
			     const int no_bins, 
			     float *jointhist, float *marghist1, 
			     float *marghist2,
			     const float smoothsize, const float fuzzyfrac)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_smoothed_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
			    jointhist,marghist1,marghist2,
			    jointentropy,margentropy1,margentropy2,
			    smoothsize,fuzzyfrac);
      float mutualinformation = margentropy1 + margentropy2 - jointentropy;
      return mutualinformation;
    }



  float normalised_mutual_info_smoothed(const volume& vref, const volume& vtest,
			       int *bindex, const Matrix& aff,
			       const float mintest, const float maxtest,
			       const int no_bins, float *jointhist, 
			       float *marghist1, float *marghist2,
			       const float smoothsize, const float fuzzyfrac)
    {
      float jointentropy=0.0, margentropy1=0.0, margentropy2=0.0;
      calc_smoothed_entropy(vref,vtest,bindex,aff,mintest,maxtest,no_bins,
			    jointhist,marghist1,marghist2,
			    jointentropy,margentropy1,margentropy2,
			    smoothsize,fuzzyfrac);
      float normmi;
      if (fabs(jointentropy)<1e-9) {
	normmi = 0.0;  // BAD registration result
      } else {
	normmi = (margentropy1 + margentropy2)/jointentropy;
      }
      return normmi;
    }


  ///////////////////////////////////////////////////////////////////////

  // Supporting interfaces


  float normcorr(const imagepair* ims, const Matrix& aff) 
    {
      return normcorr(ims->refvol,ims->testvol,aff);
    }

  float normcorr_smoothed(const imagepair* ims, const Matrix& aff) 
    {
      return normcorr_smoothed(ims->refvol,ims->testvol,aff, ims->smoothsize);
    }

  float normcorr_smoothed_sinc(const imagepair* ims, const Matrix& aff) 
    {
      return normcorr_smoothed_sinc(ims->refvol,ims->testvol,aff, ims->smoothsize);
    }


  float leastsquares(const imagepair* ims, const Matrix& aff) 
    {
      return leastsquares(ims->refvol,ims->testvol,aff);
    }

  float leastsquares_smoothed(const imagepair* ims, const Matrix& aff) 
    {
      return leastsquares_smoothed(ims->refvol,ims->testvol,aff, 
				   ims->smoothsize);
    }


  float woods_fn(const imagepair* ims, const Matrix& aff) 
    {
      return woods_fn(ims->refvol,ims->testvol,ims->bindex,aff,
		      ims->no_bins);
    }

  float woods_fn_smoothed(const imagepair* ims, const Matrix& aff) 
    {
      return woods_fn_smoothed(ims->refvol,ims->testvol,ims->bindex,aff,
		      ims->no_bins, ims->smoothsize);
    }


  float corr_ratio(const imagepair* ims, const Matrix& aff) 
    {
      return corr_ratio(ims->refvol,ims->testvol,ims->bindex,aff,
			ims->no_bins);
    }
  
  float corr_ratio_smoothed(const imagepair* ims, const Matrix& aff) 
    {
      return corr_ratio_smoothed(ims->refvol,ims->testvol,ims->bindex,aff,
			ims->no_bins, ims->smoothsize);
    }
  

  float mutual_info(imagepair* ims, const Matrix& aff)
    {
      return mutual_info(ims->refvol,ims->testvol,ims->bindex,aff,
			 ims->testmin,ims->testmax,
			 ims->no_bins,ims->plnp,ims->jointhist,
			 ims->marghist1,ims->marghist2);
    }


  float mutual_info_smoothed(imagepair* ims, const Matrix& aff)
    {
      return mutual_info_smoothed(ims->refvol,ims->testvol,
				  ims->bindex,aff,
				  ims->testmin,ims->testmax,
				  ims->no_bins,ims->fjointhist,
				  ims->fmarghist1,ims->fmarghist2,
				  ims->smoothsize, ims->fuzzyfrac);
    }


  float normalised_mutual_info(imagepair* ims, const Matrix& aff)
    {
      return normalised_mutual_info(ims->refvol,ims->testvol,ims->bindex,aff,
			 ims->testmin,ims->testmax,
			 ims->no_bins,ims->plnp,ims->jointhist,
			 ims->marghist1,ims->marghist2);
    }

  float normalised_mutual_info_smoothed(imagepair* ims, const Matrix& aff)
    {
      return normalised_mutual_info_smoothed(ims->refvol,ims->testvol,
					     ims->bindex,aff,
					     ims->testmin,ims->testmax,
					     ims->no_bins,ims->fjointhist,
					     ims->fmarghist1,ims->fmarghist2,
					     ims->smoothsize, ims->fuzzyfrac);
    }


  ///////////////////////////////////////////////////////////////////////////

#ifndef NO_NAMESPACE
 }
#endif



