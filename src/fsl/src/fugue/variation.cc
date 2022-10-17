/*  variation.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000/2001 University of Oxford  */

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

#include "unwarpfns.h"
#include "utils/options.h"
#include <map>

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;

string title="variation (Version 2.0 in c# minor)\nPhase Region Expanding Labeller for Unwrapping Discrete Estimates\nCopyright(c) 2001, University of Oxford (Mark Jenkinson)";
string examples="variation -c <rawphase> -o <unwrappedphase> [options]\nvariation -p <phasevol> -a <absvol> -o <unwrappedphase> [options]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> labelslices(string("--labelslices"), false,
		    string("does label processing in 2D (slice at a time)"),
		    false, no_argument);
Option<bool> allslices(string("-s,--slices"), false,
		       string("does all processing in 2D (slice at a time)"),
		       false, no_argument);
Option<bool> force3D(string("-f,--force3D"), false,
		     string("forces all processing to be full 3D"),
		     false, no_argument);
Option<int> connectivity(string("--connectivity"), 6,
			 string("edge connectivity of voxels"),
			 false, requires_argument);
Option<int> num_phase_split(string("-n,--numphasesplit"), 6,
			      string("number of phase partitions to use"),
			      false, requires_argument);
Option<int> startimno(string("--start"), 0,
		      string("first image number to process (default 0)"),
		      false, requires_argument);
Option<int> endimno(string("--end"), 0,
		    string("final image number to process (default Inf)"),
		    false, requires_argument);
Option<float> thresh(string("-t,--thresh"), 0,
		    string("intensity threshold for masking"),
		    false, requires_argument);
Option<string> complexvol(string("-c,--complex"), string(""),
			  string("filename of complex phase input volume"),
			  false, requires_argument);
Option<string> absvol(string("-a,--abs"), string(""),
		      string("filename of absolute input volume"),
		      false, requires_argument);
Option<string> phasevol(string("-p,--phase"), string(""),
			string("filename of raw phase input volume"),
			false, requires_argument);
Option<string> uphasevol(string("-o,--out,-u,--unwrap"), string(""),
			 string("filename for saving the unwrapped phase output"),
			 true, requires_argument);
Option<string> rawphasevol(string("-r,--rawphase"), string(""),
			   string("filename for saving the raw phase output"),
			   false, requires_argument);
Option<string> labelvol(string("-l,--labels"), string(""),
			   string("filename for saving the area labels output"),
			   false, requires_argument);
Option<string> maskvol(string("--savemask"), string(""),
			   string("filename for saving the mask volume"),
			   false, requires_argument);
Option<string> inmask(string("-m,--mask"), string(""),
		      string("filename of mask input volume"),
		      false, requires_argument);

////////////////////////////////////////////////////////////////////////////

float wrap_phase(float theta)
{
  if ( (theta<M_PI) && (theta>-M_PI) ) return theta;
  float n=MISCMATHS::round(theta/(2.0*M_PI));
  return (theta - n*2.0*M_PI);
}


void wrap_phase(volume<float>& phasevol)
{
  for (int z=phasevol.minz(); z<=phasevol.maxz(); z++) {
    for (int y=phasevol.miny(); y<=phasevol.maxy(); y++) {
      for (int x=phasevol.minx(); x<=phasevol.maxx(); x++) {
	phasevol(x,y,z) = wrap_phase(phasevol(x,y,z));
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////

typedef ColumnVector coord;
typedef ColumnVector coordpair;

bool costless(const float& c1, const float& c2) {
  return (c1<c2);
}

typedef multimap<const float,coordpair> multimapc;
typedef multimap<const float,coordpair>::iterator multimapit;


void setv(volume<float>& v, const coord& x, float val) {
  v(MISCMATHS::round(x(1)),MISCMATHS::round(x(2)),
    MISCMATHS::round(x(3))) = val;
}

float getv(const volume<float>& v, const coord& x) {
  return v(MISCMATHS::round(x(1)),MISCMATHS::round(x(2)),
	   MISCMATHS::round(x(3)));
}

bool valid(const coord& x1, const volume<float>& wrapval) {
  return fabs(getv(wrapval,x1)-0.3)>0.2;
}

void breakpair(const coordpair& cp, coord& x1, coord& x2) {
  x1 << cp(1) << cp(2) << cp(3);
  x2 << cp(4) << cp(5) << cp(6);
}

void addpt(const coord& x, const coord& x0, multimapc& costmap, 
	   const volume<float>& wrapval, const volume<float>& ph) {
  if (!valid(x0,wrapval)) { cerr << "WARNING::Invalid x0 in addpt()" << endl; return; }
  if (valid(x,wrapval)) { return; } // do nothing
  float cost = fabs(getv(ph,x) - getv(ph,x0));
  cost  = fabs(cost - 2.0*M_PI*MISCMATHS::round(cost/(2.0*M_PI)));
  coordpair cp(6);
  cp << x(1) << x(2) << x(3) << x0(1) << x0(2) << x0(3);
  costmap.insert(pair<const float, coordpair>(cost,cp));
  }

void addneighbours(const coord& x, multimapc& costmap, 
		   const volume<float>& wrapval, const volume<float>& ph) {
  int x0=MISCMATHS::round(x(1)), y0=MISCMATHS::round(x(2)), 
    z0=MISCMATHS::round(x(3));
  if (connectivity.value()==6) {
    coord newx(3);
    newx << Max(x0-1,ph.minx()) << y0 << z0;
    addpt(newx,x,costmap,wrapval,ph);
    newx << Min(x0+1,ph.maxx()) << y0 << z0;
    addpt(newx,x,costmap,wrapval,ph);
    newx << x0 << Max(y0-1,ph.miny()) << z0;
    addpt(newx,x,costmap,wrapval,ph);
    newx << x0 << Min(y0+1,ph.maxy()) << z0;
    addpt(newx,x,costmap,wrapval,ph);
    newx << x0 << y0 << Max(z0-1,ph.minz());
    addpt(newx,x,costmap,wrapval,ph);
    newx << x0 << y0 << Min(z0+1,ph.maxz());
    addpt(newx,x,costmap,wrapval,ph);
  } else if (connectivity.value()==26) { // 26 neighbour connectivity
    for (int z1=Max(z0-1,ph.minz()); z1<=Min(z0+1,ph.maxz()); z1++) {
      for (int y1=Max(y0-1,ph.miny()); y1<=Min(y0+1,ph.maxy()); y1++) {
	for (int x1=Max(x0-1,ph.minx()); x1<=Min(x0+1,ph.maxx()); x1++) {
	  coord newx(3);
	  newx << x1 << y1 << z1;
	  addpt(newx,x,costmap,wrapval,ph);
	}
      }
    }
  } else {
    cerr << "Unsupported connectivity value: " << connectivity.value() << endl;
    imthrow("Incorrect Option Setting",2);
  }
}

///////////////////////////////////////////////////////////////////////////

volume<float> ALTunwrap(const volume<float>& phin, const volume<float>& mask,
			bool verbose) 
{
  volume<float> ph(phin);
  ph *= mask;
  ColumnVector ramps;
  cerr << "Removing linear ramps" << endl;
  ramps = estimate_linear_ramps(ph,mask);
  cerr << "Estimated linear ramps" << endl;
  remove_linear_ramps(ramps,ph,mask);
  cerr << "Removed linear ramps" << endl;

  multimapc costmap;
  multimapit cit;
  volume<float> wrapval(ph);
  wrapval = 0.3;  // initially set to this value = not done
  coord x1(3), x2(3);

  // select random seed point
  cerr << "Selecting seed point" << endl;
  coord seedpoint(3);
  seedpoint << MISCMATHS::round((ph.minx() + ph.maxx())/2.0) 
	    << MISCMATHS::round((ph.miny() + ph.maxy())/2.0)
	    << MISCMATHS::round((ph.minz() + ph.maxz())/2.0);
  if (getv(mask,seedpoint)<0.5) {
    cerr << "WARNING: COV seed not in mask" << endl; // choose another point
    for (int z=mask.minz(); z<mask.maxz(); z++) {
      for (int y=mask.miny(); y<mask.maxy(); y++) {
	for (int x=mask.minx(); x<mask.maxx(); x++) {
	  if (mask(x,y,z)>0.5) { seedpoint << x << y << z; }
	}
      }
    }
  }
  setv(wrapval,seedpoint,0.0);
  cerr << "Adding neighbours" << endl;
  addneighbours(seedpoint,costmap,wrapval,ph);

  long int numit=1;
  while (!costmap.empty()) {
    if (verbose) { 
      if ((numit % 100) == 1) { 
	 cout << "Iteration: " << numit << endl;
	 cout << "Size = " << costmap.size() << endl;
	 cout << "Least cost = " << (*costmap.begin()).first << endl;
	 cout << "Greatest cost = " << (*(--(costmap.end()))).first << endl;
	 }
       numit++;
    }
    cit = costmap.begin();  // find least cost
    breakpair((*cit).second,x1,x2);  // get coordinate pair
    if ((!valid(x1,wrapval)) && (valid(x2,wrapval))) {
      setv(wrapval,x1,MISCMATHS::round( ( getv(ph,x2) - getv(ph,x1) ) / (2.0*M_PI) ) + getv(wrapval,x2));
      assert(fabs(getv(ph,x2)-getv(ph,x1)+2.0*M_PI*(getv(wrapval,x2)-getv(wrapval,x1)))<M_PI);
      addneighbours(x1,costmap,wrapval,ph); 
    }
    else if ((!valid(x2,wrapval)) && (valid(x1,wrapval))) {
      setv(wrapval,x2,MISCMATHS::round(( getv(ph,x1) - getv(ph,x2) ) / (2.0*M_PI) ) + getv(wrapval,x1));
      assert(fabs(getv(ph,x2)-getv(ph,x1)+2.0*M_PI*(getv(wrapval,x2)-getv(wrapval,x1)))<M_PI);
      addneighbours(x2,costmap,wrapval,ph); 
    }
    costmap.erase(cit);
  }

  volume<float> uph;
  uph = (ph + ((float) (2.0*M_PI))*wrapval);
  restore_linear_ramps(ramps,uph,mask);
  return uph;
}

///////////////////////////////////////////////////////////////////////////

volume<float> ALTunwrap2D(const volume<float>& phasemap, 
			  const volume<float>& mask, bool verbose)
{
  // ONLY USES 2D CONNECTIVITY

  int zmin = phasemap.minz(), zmax = phasemap.maxz();
  phasemap.activateROI();
  mask.activateROI();

  volume<float> uphase;
  uphase = phasemap;
  uphase = 0;
  uphase.activateROI();

  volume<float> pslice, uslice, olduslice, udiff, mslice, oldmslice;
  volume<float> lslice;

  for (int z=zmin; z<=zmax; z++) {
    if (verbose) cout << "SLICE NUMBER " << z << endl;

    if (z>zmin) {
      olduslice = uslice;
      oldmslice = mslice;
    }
    // get slices from phasemap and mask
    phasemap.setROIlimits(phasemap.minx(),phasemap.miny(),z,
			  phasemap.maxx(),phasemap.maxy(),z);
    pslice = phasemap.ROI();
    mask.setROIlimits(mask.minx(),mask.miny(),z,
		       mask.maxx(),mask.maxy(),z);
    lslice = mask.ROI();
    copyconvert(lslice,mslice);
    mslice.binarise(0.5);
    // unwrap 
    uslice = ALTunwrap(pslice,lslice,false); // set verbose off for unwrap
    // fix any potential phase difference
    if (z>zmin) { 
      udiff = uslice - olduslice;
      volume<float> combmslice;
      combmslice = mslice * oldmslice;
      udiff *= combmslice;
      float avpdiff;
      // Using the mean phase difference
//        avpdiff = round(udiff.sum() / combmslice.sum() / (2.0*M_PI));
      // Using the median phase difference
      float percent = 1.0 - combmslice.sum() / combmslice.nvoxels() / 2.0;
      avpdiff = MISCMATHS::round(udiff.percentile(percent) / (2.0*M_PI));
      if (avpdiff != 0.0) {
	uslice -= ((float) (avpdiff*2.0*M_PI))*mslice;
      }
    }
    // set the corresponding slice of uphase to uslice
    uphase.setROIlimits(uphase.minx(),uphase.miny(),z,
			uphase.maxx(),uphase.maxy(),z);
    uphase.copyROIonly(uslice);
  }
  
  // restore previous ROI status...
  phasemap.deactivateROI();
  mask.deactivateROI();
  
  return uphase;
}

///////////////////////////////////////////////////////////////////////////

int do_unwrapping() 
{

  volume4D<float> phasemaps, absmaps, masks;
  volumeinfo vinfo;
  
  if (complexvol.set()) {
    volume4D<float> rvol, ivol;
    read_complexvolume4D(rvol,ivol,complexvol.value(),vinfo);
    for (int n=0; n<rvol.tsize(); n++) {
      phasemaps.addvolume(phase(rvol[n],ivol[n]));
      absmaps.addvolume(abs(rvol[n],ivol[n]));
    }
  } else {
    if (verbose.value()) cerr << "Loading volumes" << endl;
    read_volume4D(phasemaps,phasevol.value(),vinfo);
    if (verbose.value()) cerr << "Phase loaded" << endl;
    read_volume4D(absmaps,absvol.value());
    if (verbose.value()) cerr << "Magnitude loaded" << endl;
  }

  bool threshset = thresh.set();
  float threshval = thresh.value();
  int imgstart=0, imgend, maskend;
  imgend = phasemaps.maxt(); 
  if (startimno.set())  imgstart = startimno.value();
  if (endimno.set())    imgend = endimno.value();
  if (inmask.set()) {
    read_volume4D(masks,inmask.value());
    if (verbose.value()) cerr << "Mask loaded" << endl;
    maskend = masks.maxt();
    if (!threshset) {
      threshset = true;
      threshval = 0.5;
    }
  }

  // sanity check on phase values
  for (int n=imgstart; n<=imgend; n++) {
    if ((phasemaps[n].max() - phasemaps[n].min())>3.0*M_PI) {
      cerr << "ERROR: input phase image exceeds allowable phase range."
	   << endl << "Allowable range is 6.283 radians.  Image range is: " 
	   << phasemaps[n].max() - phasemaps[n].min() << " radians."
	   << endl << "Aborting." << endl;
      return -1;
    }
    if ( (phasemaps[n].max()>1.001*M_PI) || (phasemaps[n].min()<-1.001*M_PI) )
      {
	if (verbose.value()) 
	  { cout << "Rewrapping phase range to [-pi,pi]" << endl; }
	wrap_phase(phasemaps[n]);
      }
  }
  
  bool label_2D = labelslices.value();
  bool unwrap_2D = allslices.value();
  // The automatic behaviour is to use 2D labels for large images
  if ((phasemaps[imgstart].xsize() > 64) || (phasemaps[imgstart].ysize() > 64))
    {
      label_2D = true;
    }
  if (force3D.value()) { label_2D = false;  unwrap_2D = false; }

  volume4D<float> uphase;

  for (int n=imgstart; n<=imgend; n++) {

    // Make mask
    volume<float> mask;
    if (inmask.set()) {
      mask = masks[Min(n,maskend)];
    } else {
      mask = absmaps[n];
    }

    if (!threshset) { threshval = basic_mask_threshold(mask); }

    if (unwrap_2D) {
      mask = make_head_mask2D(mask,threshval);
    } else {
      mask = make_head_mask(mask,threshval);
    }
    if (maskvol.set()) { save_volume(mask,maskvol.value()); }

    /*
    // Make labels
    volume<int> label;
    if (verbose.value()) 
      {cout << "Number of phase splits = " << num_phase_split.value() << endl;}
    if (unwrap_2D) {
      label = find_phase_labels2D(phasemaps[n],mask,
				  num_phase_split.value(),false);
    } else if (label_2D) {
      label = find_phase_labels2D(phasemaps[n],mask,
				  num_phase_split.value(),true);
    } else {
      label = find_phase_labels(phasemaps[n],mask,num_phase_split.value());
    }
    if (labelvol.set()) { save_volume(label,labelvol.value()); }
    */

    // Unwrap phase
    volume<float> uph;
    if (unwrap_2D) {
      uph = ALTunwrap2D(phasemaps[n],mask,verbose.value());
    } else {
      uph = ALTunwrap(phasemaps[n],mask,verbose.value());
    }

    uphase.addvolume(uph);
  }

  // Save outputs
  save_volume4D(uphase,uphasevol.value(),vinfo);

  if (rawphasevol.set()) {
    save_volume4D(phasemaps,rawphasevol.value());
  }

  return 0;
}




int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(complexvol);
    options.add(absvol);
    options.add(phasevol);
    options.add(uphasevol);
    options.add(num_phase_split);
    //options.add(labelslices);
    options.add(allslices);
    options.add(force3D);
    options.add(thresh);
    options.add(inmask);
    options.add(startimno);
    options.add(endimno);
    options.add(maskvol);
    options.add(rawphasevol);
    //options.add(labelvol);
    options.add(connectivity);
    options.add(verbose);
    options.add(help);
    
    options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
    if ( ( (complexvol.set()) && (absvol.set()) ) || 
	 ( (complexvol.set()) && (phasevol.set()) ) ) 
      {
	options.usage();
	cerr << endl 
	     << "Cannot specify both --complex AND --phase or --abs."
	     << endl;
	exit(EXIT_FAILURE);
      }
 
    if ( ( (phasevol.set()) && (absvol.unset()) ) || 
	 ( (phasevol.unset()) && (absvol.set()) ) ) 
      {
	options.usage();
	cerr << endl 
	     << "Both --phase AND --abs must be used together." 
	     << endl;
	exit(EXIT_FAILURE);
      }
    
    if (num_phase_split.value() < 2)
      {
	options.usage();
	cerr << endl 
	     << "Always set --numphasesplit greater than 1." 
	     << endl << "NOT " << num_phase_split.value() << endl;
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  return do_unwrapping();
}











