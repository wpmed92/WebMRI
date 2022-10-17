/*  distancemap.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2003 University of Oxford  */

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
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

string title="distancemap (Version 1.0)\nCopyright(c) 2003, University of Oxford (Mark Jenkinson)";
string examples="distancemap [options] -i <inputimage> -o <outputimage>";


Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> invertinput(string("--invert"), false,
		  string("invert input image"),
		  false, no_argument);
Option<int> version(string("--version"), 2,
		  string("<int>\tselect algorithm version (1 or 2=default)"),
		  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input mask filename"),
		  true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output image filename"),
		  true, requires_argument);
Option<string> maskname(string("-m,--mask"), string(""),
		  string("mask image filename"),
		  false, requires_argument);
Option<string> labelname(string("-l,--localmax"), string(""),
		  string("local maxima output image filename"),
		  false, requires_argument);
Option<string> segmentname(string("-s,--segment"), string(""),
		  string("segmented output image filename"),
		  false, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Global variables (not options)

class rowentry { public: int x; int y; int z; float d; } ;

bool rowentry_lessthan(const rowentry& r1, const rowentry& r2) 
{
  return r1.d < r2.d ;
}

vector<rowentry> schedule;
Matrix octantsign(8,3);

////////////////////////////////////////////////////////////////////////////

// Local functions

int setup_globals(const volume<float>& vin)
{
  // octantsign gives the 8 different octant sign combinations for coord
  //  offsets
  int row=1;
  for (int p=-1; p<=1; p+=2) {
    for (int q=-1; q<=1; q+=2) {
      for (int r=-1; r<=1; r+=2) {
	octantsign(row,1)=p;
	octantsign(row,2)=q;
	octantsign(row,3)=r;
	row++;
      }
    }
  }

  // construct list of displacements (in one octant) in ascending
  // order of distance
  for (int z=vin.minz(); z<=vin.maxz(); z++) {
    for (int y=vin.miny(); y<=vin.maxy(); y++) {
      for (int x=vin.minx(); x<=vin.maxx(); x++) {
	if ( (x==0) && (y==0) && (z==0) ) 
         { // do nothing 
         } else {
	   rowentry newrow;
	   newrow.x=x;
	   newrow.y=y;
	   newrow.z=z;
	   float d2 = norm2(x*vin.xdim(),y*vin.ydim(),z*vin.zdim());
	   newrow.d=d2;
	   schedule.push_back(newrow);
         }
      }
    }
  }

  // sort schedule to get ascending d2
  sort(schedule.begin(),schedule.end(),rowentry_lessthan);
  return 0;
}

int find_nearest(const volume<float>& vin, int x, int y, int z,
	         int& x1, int& y1, int&z1)
{
  for (int r=0; r<(int) schedule.size(); r++) 
  {
    for (int p=1; p<=8; p++) {
      int dx=MISCMATHS::round(schedule[r].x*octantsign(p,1));
      int dy=MISCMATHS::round(schedule[r].y*octantsign(p,2));
      int dz=MISCMATHS::round(schedule[r].z*octantsign(p,3));
      if (vin.in_bounds(x+dx,y+dy,z+dz) && (vin.value(x+dx,y+dy,z+dz)>0.5)) 
       { 
	x1=x+dx; y1=y+dy; z1=z+dz; 
	return 0; 
       }
    }
  }
  // return furtherest point (in error)
  x1= x + MISCMATHS::round(schedule.back().x);
  y1= y + MISCMATHS::round(schedule.back().y);
  z1= z + MISCMATHS::round(schedule.back().z);
  return 1;
}


// find local maxima and returns volume with zero except at local max
volume<float> label_local_maxima(const volume<float>& vin, 
				  const volume<float>& mask) 
{
  int lmaxcount=0;
  volume<float> vout;
  vout = 0.0f*vin;
  for (int z=vin.minz(); z<=vin.maxz(); z++) {
    for (int y=vin.miny(); y<=vin.maxy(); y++) {
      for (int x=vin.minx(); x<=vin.maxx(); x++) {
	if ( (mask(x,y,z)>0.5) &&
	     vin(x,y,z)>vin(x-1,y-1,z-1) &&
	     vin(x,y,z)>vin(x,  y-1,z-1) &&
	     vin(x,y,z)>vin(x+1,y-1,z-1) &&
	     vin(x,y,z)>vin(x-1,y,  z-1) &&
	     vin(x,y,z)>vin(x,  y,  z-1) &&
	     vin(x,y,z)>vin(x+1,y,  z-1) &&
	     vin(x,y,z)>vin(x-1,y+1,z-1) &&
	     vin(x,y,z)>vin(x,  y+1,z-1) &&
	     vin(x,y,z)>vin(x+1,y+1,z-1) &&
	     vin(x,y,z)>vin(x-1,y-1,z) &&
	     vin(x,y,z)>vin(x,  y-1,z) &&
	     vin(x,y,z)>vin(x+1,y-1,z) &&
	     vin(x,y,z)>vin(x-1,y,  z) &&
	     vin(x,y,z)>=vin(x+1,y,  z) &&
	     vin(x,y,z)>=vin(x-1,y+1,z) &&
	     vin(x,y,z)>=vin(x,  y+1,z) &&
	     vin(x,y,z)>=vin(x+1,y+1,z) &&
	     vin(x,y,z)>=vin(x-1,y-1,z+1) &&
	     vin(x,y,z)>=vin(x,  y-1,z+1) &&
	     vin(x,y,z)>=vin(x+1,y-1,z+1) &&
	     vin(x,y,z)>=vin(x-1,y,  z+1) &&
	     vin(x,y,z)>=vin(x,  y,  z+1) &&
	     vin(x,y,z)>=vin(x+1,y,  z+1) &&
	     vin(x,y,z)>=vin(x-1,y+1,z+1) &&
	     vin(x,y,z)>=vin(x,  y+1,z+1) &&
	     vin(x,y,z)>=vin(x+1,y+1,z+1) ) {
	  lmaxcount++;
	  vout(x,y,z)=lmaxcount;
	}
      }
    }
  }
  return vout;
}


// create the distance map as vout
// if use_distance is false then put the value of the input volume
//  at the output location rather than the distance to it
void create_distancemap(volume<float>& vout, const volume<float>& vin,
			const volume<float>& mask, bool use_distance=true)
{
  int x1, y1, z1;
  for (int z=vout.minz(); z<=vout.maxz(); z++) {
    if (verbose.value()) { cerr << ":"; }
    for (int y=vout.miny(); y<=vout.maxy(); y++) {
      if (verbose.value()) { cerr << "."; }
      for (int x=vout.minx(); x<=vout.maxx(); x++) {
	if (mask(x,y,z)>0.5) {
	  find_nearest(vin,x,y,z,x1,y1,z1);
	  if (use_distance) {
	    vout(x,y,z)=sqrt(norm2((x1-x)*vin.xdim(),
				   (y1-y)*vin.ydim(),(z1-z)*vin.zdim()));
	  } else {
	    vout(x,y,z)=vin(x1,y1,z1);
	  }
	}
      }
    }
  }
  if (verbose.value()) { cerr << endl; }
}


// makes the minimum distance map for each voxel to the non-zero voxels
//  in the input image
int do_work(int argc, char* argv[]) 
{
  volume<float> vin, dmap, mask;
  read_volume(vin,inname.value());
  if (invertinput.value()) {
    vin = 1.0f - binarise(vin,0.5f);
  }
  setup_globals(vin);
  dmap = 0.0f*vin;
  if (maskname.set()) {
    read_volume(mask,maskname.value());
  } else {
    // set the mask as all the zero values (others have distance=0)
    mask = 1.0f - binarise(vin,0.5f);
  }

  if (verbose.value()) { cout << "Creating distance map" << endl; }
  create_distancemap(dmap,vin,mask);
  save_volume(dmap,outname.value());

  if (labelname.set() || segmentname.set()) {
    if (verbose.value()) { cout << "Finding local max" << endl; }
    volume<float> label;
    label = label_local_maxima(dmap,mask);
    if (labelname.set()) { save_volume(label,labelname.value()); }
    if (segmentname.set()) {
      if (verbose.value()) { cout << "Segmenting wrt distance" << endl; }
      create_distancemap(dmap,label,mask,false);
      save_volume(dmap,segmentname.value());
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(inname);
    options.add(outname);
    options.add(maskname);
    options.add(labelname);
    options.add(segmentname);
    options.add(invertinput);
    options.add(verbose);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
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

  // Call the local functions

  return do_work(argc,argv);
}

