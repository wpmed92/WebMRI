/*  invwarptest.cc

    Mark Jenkinson, FMRIB Image Analysis Group

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

#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "newimage/warpfns.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

////////////////////////////////////////////////////////////////////////////

// COMMAND LINE OPTIONS

string title="invwarptest (Version 1.2)\nCopyright(c) 2001, University of Oxford (Mark Jenkinson)";
string examples="invwarptest -w warpvol [-o invwarpvol]\ninvwarptest -w shiftmap -i inputimage -f forwardwarp";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<string> outname(string("-o,--out"), string(""),
		       string("filename for output (inverse warped) image"),
		       false, requires_argument);
Option<string> warpname(string("-w,--warp"), string(""),
			string("filename for warp/shiftmap transform (volume)"),
			true, requires_argument);
Option<string> fname(string("-f"), string(""),
			string("filename for forward warp output (volume)"),
			false, requires_argument);
Option<string> iname(string("-i"), string(""),
			string("filename for input volume (to be forward warped)"),
			false, requires_argument);


////////////////////////////////////////////////////////////////////////////


// Inverse warp stuff

bool in_tetrahedron(const ColumnVector& x, const ColumnVector& x1, 
		    const ColumnVector& x2, const ColumnVector& x3, 
		    const ColumnVector& x4)
{  // checks whether first point is inside the tetrahedron defined by
   // the other 4 vertices
  Matrix vbasis(3,3);
  vbasis.Column(1) = x2-x1;
  vbasis.Column(2) = x3-x1;
  vbasis.Column(3) = x4-x1;
  ColumnVector coeffs(3);
  if (vbasis.Determinant() < 1e-3) return false;
  coeffs = vbasis.i() * ( x - x1 );
  if ( (coeffs(1)>=0) && (coeffs(2)>=0) && (coeffs(3)>=0) &&
       (coeffs.Sum()<=1) ) {
    return true;
  } else {
    return false;
  }
}

ColumnVector midvox(const volume4D<float>& warp, int wx, int wy, int wz, 
		    const ColumnVector& dirvec)
{
  ColumnVector wmid(3);
  wmid = 0.0;
  for (int comp=1; comp<=3; comp++) {
    wmid(comp) += warp(wx,wy,wz,comp-1);
    wmid(comp) += warp(wx,wy,wz+MISCMATHS::round(dirvec(3)),comp-1);
    wmid(comp) += warp(wx,wy+MISCMATHS::round(dirvec(2)),wz,comp-1);
    wmid(comp) += warp(wx,wy+MISCMATHS::round(dirvec(2)),wz+MISCMATHS::round(dirvec(3)),comp-1);

    wmid(comp) += warp(wx+MISCMATHS::round(dirvec(1)),wy,wz,comp-1);
    wmid(comp) += warp(wx+MISCMATHS::round(dirvec(1)),wy,wz+MISCMATHS::round(dirvec(3)),comp-1);
    wmid(comp) += warp(wx+MISCMATHS::round(dirvec(1)),wy+MISCMATHS::round(dirvec(2)),wz,comp-1);
    wmid(comp) += warp(wx+MISCMATHS::round(dirvec(1)),wy+MISCMATHS::round(dirvec(2)),wz+MISCMATHS::round(dirvec(3)),comp-1);
  }
  return wmid;
}

///////////////////////////////////////////////////////////////////////////

// OLD VERSION ...

/*

	      w(1)=warp(wx,wy,wz,0); 
	      w(2)=warp(wx,wy,wz,1); 
	      w(3)=warp(wx,wy,wz,2); 
	      // find and store all adjacent voxel centres
	      Matrix wnearest(3,8);
	      Matrix dirvec(8,3);
	      dirvec << 1.0 << 1.0 << 1.0 
		     << 1.0 << 1.0 << -1.0
		     << 1.0 << -1.0 << 1.0 
		     << 1.0 << -1.0 << -1.0 
		     << -1.0 << 1.0 << 1.0 
		     << -1.0 << 1.0 << -1.0 
		     << -1.0 << -1.0 << 1.0 
		     << -1.0 << -1.0 << -1.0 ;
	      dirvec = dirvec.t();
	      ColumnVector wmid(3);
	      for (int n0=1; n0<=8; n0++) {
		wmid = midvox(warp,wx,wy,wz,dirvec.Column(n0));
		wnearest.SubMatrix(1,3,n0,n0)=wmid;
	      }
	      // loop over adjacent voxel centres (all tetrahedrons)
	      for (int n1=1; n1<=8; n1++) {
		w1 = wnearest.Column(n1);
		for (int n2=n1+1; n2<=8; n2++) {
		  w2 = wnearest.Column(n2);
		  for (int n3=n2+1; n3<=8; n3++) {
		    w3 = wnearest.Column(n3);
		    
		    if (in_tetrahedron(v,w,w1,w2,w3)) {
		      nn(wx,wy,wz)+=1.0f;
		    }



*/


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// 3D VERSIONS //
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


/*


int first_pass(const volume4D<float>& warp, volume4D<float>& invwarp,
	       volume<float>& multiplicity) 
{
  volume<float> mindist;
  mindist = invwarp[0]*0.0f;
  mindist += Sqr(warp.max());  // big dist
  multiplicity = mindist*0.0f;
  invwarp *= 0.0f;
  
  // loop over points in warp image (orig points, w)
  //  and find nearest voxel in invwarp
  for (int wz=warp.minz(); wz<=warp.maxz(); wz++) {
    cerr << "wz = " << wz << endl;
    for (int wy=warp.miny(); wy<=warp.maxy(); wy++) {
      for (int wx=warp.minx(); wx<=warp.maxx(); wx++) {
	float fx, fy, fz;
	fx = warp(wx,wy,wz,0)/warp.xdim();
	fy = warp(wx,wy,wz,1)/warp.ydim();
	fz = warp(wx,wy,wz,2)/warp.zdim();
	int x,y,z;
	x = MISCMATHS::round(fx);
	y = MISCMATHS::round(fy);
	z = MISCMATHS::round(fz);
	if (invwarp.in_bounds(x,y,z)) {
	  float dist=Sqr(fx-x) + Sqr(fy-y) + Sqr(fz-z);
	  if (dist<mindist(x,y,z)) {
	    mindist(x,y,z)=dist;
	    multiplicity(x,y,z)+=1.0f;
	    // store the original point (from warp volume)
	    invwarp.value(x,y,z,0)=wx * warp.xdim();
	    invwarp.value(x,y,z,1)=wy * warp.ydim();
	    invwarp.value(x,y,z,2)=wz * warp.zdim();
	  }
	}
      }
    }
  }

  return 0;
}


  
int second_pass(const volume4D<float>& warp, volume4D<float>& invwarp,
		volume<float>& multiplicity) 
{
  volume<float> mindist;
  mindist = invwarp[0]*0.0f;
  mindist += Sqr(warp.max());  // big dist

  // for all points which were multiply mapped onto resolve using
  //  the one most like the neighbours
  for (int wz=warp.minz(); wz<=warp.maxz(); wz++) {
    cerr << "wz = " << wz << endl;
    for (int wy=warp.miny(); wy<=warp.maxy(); wy++) {
      for (int wx=warp.minx(); wx<=warp.maxx(); wx++) {
	float fx, fy, fz;
	fx = warp(wx,wy,wz,0)/warp.xdim();
	fy = warp(wx,wy,wz,1)/warp.ydim();
	fz = warp(wx,wy,wz,2)/warp.zdim();
	int x,y,z;
	x = MISCMATHS::round(fx);
	y = MISCMATHS::round(fy);
	z = MISCMATHS::round(fz);

	if (invwarp.in_bounds(x,y,z) && (multiplicity(x,y,z)>1)) {
	  // loop over neighbours to get the average (median?)
	  float avx=0, avy=0, avz=0, n=0.0;
	  for (int z0=z-1; z0<=z+1; z0++) {
	    for (int y0=y-1; y0<=y+1; y0++) {
	      for (int x0=x-1; x0<=x+1; x0++) {
		if (multiplicity(x0,y0,z0)==1) {
		  n+=1.0;
		  avx+=invwarp(x0,y0,z0,0); 
		  avy+=invwarp(x0,y0,z0,1); 
		  avz+=invwarp(x0,y0,z0,2);
		}
	      }
	    }
	  }
	  if (n>0.5) { 
	    avx/=n;  avy/=n;  avz/=n;
	    float dist=Sqr(wx*warp.xdim()-avx)+Sqr(wy*warp.ydim()-avy)
	      +Sqr(wz*warp.zdim()-avz);
	    if (dist<mindist(x,y,z)) {
	      mindist(x,y,z)=dist;
	      // store the original point (from warp volume)
	      invwarp(x,y,z,0)=wx * warp.xdim();
	      invwarp(x,y,z,1)=wy * warp.ydim();
	      invwarp(x,y,z,2)=wz * warp.zdim();
	    }
	  }

	  
	}
      }
    }
  }

  return 0;
}



int third_pass(const volume4D<float>& warp, volume4D<float>& invwarp,
		volume<float>& multiplicity) 
{
  return 0;
}


int refine_points(const volume4D<float>& warp, volume4D<float>& invwarp,
		  const volume<float>& multiplicity) 
{

  if (warp.novoxels()>0) return 0;   // temporary abort until this is finished
  volume<float> mindist;
  mindist = invwarp[0]*0.0f;
  mindist += Sqr(warp.max());  // big dist

  float wx, wy, wz, fx, fy, fz, dist, mindist, px, py, pz;

  // 
  for (int vz=invwarp.minz(); vz<=invwarp.maxz(); vz++) {
    cerr << "vz = " << vz << endl;
    for (int vy=invwarp.miny(); vy<=invwarp.maxy(); vy++) {
      for (int vx=invwarp.minx(); vx<=invwarp.maxx(); vx++) {
	wx = invwarp(vx,vy,vz,0);
	wy = invwarp(vx,vy,vz,1);
	wz = invwarp(vx,vy,vz,2);
	fx = warp[0].interpolate(wx,wy,wz);
	fy = warp[1].interpolate(wx,wy,wz);
	fz = warp[2].interpolate(wx,wy,wz);


	// UNFINISHED STUFF ...

	if (invwarp.in_bounds(x,y,z) && (multiplicity(x,y,z)==0)) {
	  mindist=Sqr(wx*warp.xdim()-invwarp(x,y,z))+Sqr(wy*warp.ydim()-invwarp(x,y,z))+Sqr(wz*warp.zdim()-invwarp(x,y,z));
	  // search local neighbourhood for better points

	  for (int n=1; n<=8; n++) {
	    px = wx + ptoffset(n,1);
	    py = wy + ptoffset(n,2);
	    pz = wz + ptoffset(n,3);

	    fx = warp[0].interpolate(wx,wy,wz);
	    fy = warp[1].interpolate(wx,wy,wz);
	    fz = warp[2].interpolate(wx,wy,wz);

	    dist=Sqr(wx*warp.xdim()-invwarp.interpolate(px,py,pz))+Sqr(wy*warp.ydim()-invwarp.interpolate(px,py,pz))+Sqr(wz*warp.zdim()-invwarp.interpolate(px,py,pz));
	    if (dist<mindist) {
	      bestptx = px;  bestpty = py; bestptz = pz;
	    }
	  }
	  
	  invwarp(x,y,z,0)=;
	  
	}
	  
      }
    }
  }

  return 0;
}


*/

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// 1D VERSIONS //
////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



int first_pass(const volume<float>& shiftmap, volume<float>& invshift,
	       volume<float>& multiplicity, const volume<float>& mask) 
{
  // ALL SHIFTMAPS IN VOXELS, NOT MM

  // returns number of solutions to y = p + invshift(p)
  // where y + shiftmap(y) = p
  // multiplicity returns number of solutions (0,1,>1) to be fixed
  //  for (0,>1) cases in later passes
  invshift = shiftmap * 0.0f;
  multiplicity = mask - 1.0;
  volume<float> wvol = shiftmap * 0.0f;
  // invshift.extrapolationmethod(assert);  // DEBUG
  
  // loop over points in line of warp image (orig points, w)
  //  and find nearest voxel in invshift
  for (int wz=shiftmap.minz(); wz<=shiftmap.maxz(); wz++) {
    cerr << "wz = " << wz << endl;
    for (int wx=shiftmap.minx(); wx<=shiftmap.maxx(); wx++) {
      for (int wy=shiftmap.miny(); wy<shiftmap.maxy(); wy++) {  // NB less than
	if ( (mask(wx,wy,wz)>0.5) && (mask(wx,wy+1,wz)>0.5) ) {
	  float p0 = wy + shiftmap(wx,wy,wz);
	  float p1 = wy + 1.0 + shiftmap(wx,wy+1,wz);
	  
	  float pminf = Min(p0,p1);
	  float pmaxf = Max(p0,p1);
	  int pmin=(int) Min(ceil(p0),ceil(p1));
	  int pmax=(int) Max(floor(p0),floor(p1));
	  
	  for (int p=pmin; p<=pmax; p++) {
	    multiplicity(wx,p,wz) += 1.0f;
	    if (fabs(p1-p0)>1e-8) {
	      // weighting of intensity modulation at this point
	      float weight=fabs(p - pminf) / Sqr(pmaxf - pminf); 
	      if (wvol(wx,p,wz) < weight) {
		wvol(wx,p,wz) = weight;
		float alpha = (p-p0)/(p1-p0);
		invshift(wx,p,wz) = wy + alpha - p;
	      }
	    } else {
	      invshift(wx,p,wz) = 0.0;
	    }
	  }

	}

      }
    }
  }
  
  return 0;
}

  

float invpt(const volume<float>& shiftmap, int wx, int wy, int wz, 
	    float avy, const volume<float>& mask)
{
  // returns inverse point (to be stored in invshift(wx,wy,wz)) where
  //  y + invshift(y) is closest to avy
  float soln = 0.0, cand, p=avy;
  for (int y=shiftmap.miny(); y<shiftmap.maxy(); y++) {  // NB less than
    if ( (mask(wx,y,wz)>0.5) && (mask(wx,y+1,wz)>0.5) ) {
      float p0 = y + shiftmap(wx,y,wz);
      float p1 = y + 1.0 + shiftmap(wx,y+1,wz);
      
      int pmin=(int) Min(ceil(p0),ceil(p1));
      int pmax=(int) Max(floor(p0),floor(p1));
      
      if ( (wy>=pmin) && (wy<=pmax) ) {
	p = wy;
	if (fabs(p1-p0)>1e-8) {
	  float alpha = (p-p0)/(p1-p0);
	  cand = wy + alpha;
	} else {
	  cand = p;
	}
	if ( fabs(cand - avy) < fabs(soln - avy) ) {
	  soln = cand;
	}
      }
 
    }
  }
  
  return (soln - p);  // use relative shift (not absolute coord)
}


int second_pass(const volume<float>& shiftmap, volume<float>& invshift,
		volume<float>& multiplicity, const volume<float>& mask) 
{
  // Sorts out multiplicity > 1 - needs to be iterated until return=0
  // after fixing mult > 1, -1*mult is put in the corresponding spot
  // Therefore mult=1 or mult<0 can be used as a "good" point
  int nonfixed=0;
  for (int wz=invshift.minz(); wz<=invshift.maxz(); wz++) {
    cerr << "wz = " << wz << endl;
    for (int wx=invshift.minx(); wx<=invshift.maxx(); wx++) {
      for (int wy=invshift.miny(); wy<=invshift.maxy(); wy++) {
	if (multiplicity(wx,wy,wz)>1.1) {
	  float avy=0.0;  int nn=0;
	  for (int nz=Max(wz-1,0); nz<=Min(wz+1,invshift.maxz()); nz++) {
	    for (int ny=Max(wy-1,0); ny<=Min(wy+1,invshift.maxy()); ny++) {
	      for (int nx=Max(wx-1,0); nx<=Min(wx+1,invshift.maxx()); nx++) {
		if ( (fabs(multiplicity(nx,ny,nz)-1.0)<1e-8) || 
		     (multiplicity(nx,ny,nz)<-1e-8) ) {
		  // mult=1 or mult<0 (not mult=0)
		  avy+=ny+invshift(nx,ny,nz);
		  nn++;
		}
	      }
	    }
	  }
	  if (nn>=1) {
	    // found some good neighbours, so fix it!
	    multiplicity(wx,wy,wz)=-multiplicity(wx,wy,wz);
	    avy /= (float) nn;
	    // find solution that is closest to avy
	    invshift(wx,wy,wz) = invpt(shiftmap,wx,wy,wz,avy,mask);
	  } else {
	    nonfixed++;
	  }
	}
      }
    }
  }
  return nonfixed;
}



int third_pass(const volume<float>& shiftmap, volume<float>& invshift,
		volume<float>& multiplicity, const volume<float>& mask) 
{
  // Sorts out multiplicity = 0
  // Extrapolate using the overall stretch of the shiftmap
  for (int wz=invshift.minz(); wz<=invshift.maxz(); wz++) {
    cerr << "wz = " << wz << endl;
    for (int wx=invshift.minx(); wx<=invshift.maxx(); wx++) {
      float scale, cy, cp, y0=invshift.maxy()+1, y1=invshift.miny()-1;
      float p0=y0, p1=y1;
      for (int wy=invshift.miny(); wy<=invshift.maxy(); wy++) {
	// calculate global stretch in this line (average over lines?)
	// use this to guess the extrapolated values...
	if (fabs(multiplicity(wx,wy,wz)-1.0)<1e-8) {
	  if (wy<y0) {
	    y0 = wy;
	    p0 = wy + invshift(wx,wy,wz);
	  }
	  if (wy>y1) {
	    y1 = wy;
	    p1 = wy + invshift(wx,wy,wz);
	  }
	}
      }
      if (fabs(y1-y0)<1e-8) { // avoid crash
	y0=p0=invshift.miny();
	y1=p1=invshift.maxy();
      }
      scale = (p1-p0)/(y1-y0);
      cp = (p0+p1)/2.0;
      cy = (y0+y1)/2.0;
      for (int wy=invshift.miny(); wy<=invshift.maxy(); wy++) {
	if (fabs(multiplicity(wx,wy,wz))<1e-8) {  // mult = 0
	  invshift(wx,wy,wz) = scale*(wy-cy) + cp - wy;
	}
      }
    }
  }
  return 0;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



volume<float> forward_pixshift(const volume<float>& invol,
			       const volume<float>& pixshiftmap)
{
  volume<float> outvol;
  if (!samesize(invol,pixshiftmap)) {
    cerr << "Forward_pixshift:: ERROR - invol and pixshiftmap are different sizes" << endl;
    exit(1);
  }
  // find background value to use throughout
  float backgnd=invol.backgroundval();
  outvol = invol;
  outvol = backgnd;

  pixshiftmap.setextrapolationmethod(mirror);

  for (int z=0; z<invol.zsize(); z++) {
    cerr << "z = " << z << endl;
    for (int x=0; x<invol.xsize(); x++) {

      for (int y=0; y<invol.ysize()-1; y++) {
	int yA=y, yB=y+1;
	float involA = invol(x,yA,z);
	float involB = invol(x,yB,z);
	float pA = yA + pixshiftmap(x,yA,z);
	float pB = yB + pixshiftmap(x,yB,z);
	float pmin = Min(pA,pB);
	float pmax = Max(pA,pB);
	int ipmin = (int) ceil(pmin);
	int ipmax = (int) floor(pmax);
	if ((pmax - pmin)>1e-8) {
	  for (int p=ipmin; p<=ipmax; p++) {
	    float weight = (p - pA) / (pB - pA);
	    // weight /= fabs(pmax - pmin);  // intensity correction
	    // sum both forward and backward kernels
	    outvol(x,p,z) += involA + weight * (involB - involA) - backgnd;
	  }
	}
      }
    }
  }
  return outvol;
}


///////////////////////////////////////////////////////////////////////////

int invert_warp1D(const volume<float>& warp, volume<float>& invwarp)
{
  volume<float> multiplicity;
  volume<float> mask;
  mask = abs(warp);
  mask.binarise(0.001);

  first_pass(warp,invwarp,multiplicity,mask);
  // can trust all points in invwarp where multiplicity==1


  // removed the next pass as it was worse than using maximum 
  //  projection weight (current version of the first_pass)
  /*
  int nbad=0, oldnbad=1;
  while (nbad < oldnbad) {
    oldnbad = nbad;
    nbad = second_pass(warp,invwarp,multiplicity,mask);
  }
  if (nbad>0) {
    cout << "Could not resolve " << nbad << " points" << endl;
  }
  // can now trust all (except nbad) points in invwarp where multiplicity>=1
  */

  third_pass(warp,invwarp,multiplicity,mask);
  // can now trust points where multiplicity==0

  // TEMPORARY DEBUG
  save_volume(multiplicity,"multiplicity");


  return 0;

}


///////////////////////////////////////////////////////////////////////////


int invwarptest()
{

  // read in images
  volume<float> invwarp;
  volumeinfo vinfo;
  volume<float> warpvol;

  read_volume(warpvol,warpname.value(),vinfo);
  //  if (warpvol.tsize()<3) {
  //  cerr << "Warp volume does not have the correct dimensions (min. 3 volumes)"
  //	 << endl;
  //    exit(1);
  //  }

  if (iname.set()) {
    volume<float> invol;
    read_volume(invol,iname.value());
    invwarp = forward_pixshift(invol,warpvol);
  } else {
    // set up same dimensions...
    invwarp = warpvol * 0.0f;
    invert_warp1D(warpvol,invwarp);
  }

  if (outname.set()) {
    // save the results
    save_volume(invwarp,outname.value(),vinfo);
  }
  return 0;
}




int main(int argc, char *argv[])
{

  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(warpname);
    options.add(outname);
    options.add(fname);
    options.add(iname);
    options.add(verbose);
    options.add(help);
    
    options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  return invwarptest();
}

