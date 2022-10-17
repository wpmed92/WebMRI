/*  avwstats++.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2003-2005 University of Oxford  */

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


// Like avwstats but better!  :)

#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "newimage/costfns.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

bool masks_used=false;
bool lthr_used=false;
bool uthr_used=false;

void print_usage(const string& progname) {
  cout << "Usage: avwstats++ <input> [options]" << endl << endl;
  cout << endl;
  cout << "-l <lthresh> : set lower threshold" << endl;
  cout << "-u <uthresh> : set upper threshold" << endl;
  cout << "-r           : output <robust min intensity> <robust max intensity>" << endl;
  cout << "-R           : output <min intensity> <max intensity>" << endl;
  cout << "-e           : output mean entropy ; mean(-i*ln(i))" << endl;
  cout << "-E           : output mean entropy (of nonzero voxels)" << endl;
  cout << "-v           : output <voxels> <volume>" << endl;
  cout << "-V           : output <voxels> <volume> (for nonzero voxels)" << endl;
  cout << "-m           : output mean" << endl;
  cout << "-M           : output mean (for nonzero voxels)" << endl;
  cout << "-s           : output standard deviation" << endl;
  cout << "-S           : output standard deviation (for nonzero voxels)" << endl;
  cout << "-w           : output smallest ROI <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> containing nonzero voxels" << endl;
  cout << "-x           : output co-ordinates of maximum voxel" << endl;
  cout << "-X           : output co-ordinates of minimum voxel" << endl;
  cout << "-c           : output centre-of-gravity (cog) in mm coordinates" << endl;
  cout << "-C           : output centre-of-gravity (cog) in voxel coordinates" << endl;
  cout << "-p <n>       : output nth percentile (n between 0 and 100)" << endl;
  cout << "-P <n>       : output nth percentile (for nonzero voxels)" << endl;
  cout << "-a           : use absolute values of all image intensities"<< endl;
  cout << "-k <mask>    : use the specified image (filename) for masking - overrides lower and upper thresholds" << endl;
  cout << "-h <nbins>   : output a histogram (for the thresholded/masked voxels only) with nbins" << endl; 
  cout << "-H <nbins> <min> <max>   : output a histogram (for the thresholded/masked voxels only) with nbins and histogram limits of min and max" << endl; 
  cout << endl;
  cout << "Note - thresholds are not inclusive ie lthresh<allowed<uthresh" << endl;
}


// Some specialised nonzero functions just for speedup
//  (it avoids generating masks when not absolutely necessary)

long int nonzerocount(const volume4D<float>& vol)
{
  long int totn=0;
  for (int t=vol.mint(); t<=vol.maxt(); t++) {
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (vol(x,y,z,t)!=0.0) {
	    totn++;
	  }
	}
      }
    }
  }
  return totn;
}

double nonzeromean(const volume4D<float>& vol)
{
  double totv=0.0;
  long int totn=0;
  for (int t=vol.mint(); t<=vol.maxt(); t++) {
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (vol(x,y,z,t)!=0.0) {
	    totv+=(double) vol(x,y,z,t);
	    totn++;
	  }
	}
      }
    }
  }
  double meanval=0.0;
  if (totn>0) {
    meanval=totv/totn;
  }
  return meanval;
}

double nonzerostddev(const volume4D<float>& vol)
{
  double totv=0.0, totvv=0.0;
  long int totn=0;
  for (int t=vol.mint(); t<=vol.maxt(); t++) {
    for (int z=vol.minz(); z<=vol.maxz(); z++) {
      for (int y=vol.miny(); y<=vol.maxy(); y++) {
	for (int x=vol.minx(); x<=vol.maxx(); x++) {
	  if (vol(x,y,z,t)!=0.0) {
	    float v=vol(x,y,z,t);
	    totvv+=(double) v*v;
	    totv+=(double) v;
	    totn++;
	  }
	}
      }
    }
  }
  double var=0.0;
  if (totn>1) {
    double meanval = totv/totn;
    var = (totvv - totn*meanval*meanval)/(totn-1);
  }
  return std::sqrt(var);
}

int generate_masks(volume4D<float> &mask, volume4D<float> &masknz, const volume4D<float> &vin,
		   float& lthr, float& uthr) 
{
  if (!lthr_used) { lthr=vin.min()-1; }
  if (!uthr_used) { uthr=vin.max()+1; }
  mask = binarise(vin,lthr,uthr,exclusive);
  masknz = mask * (1.0f - binarise(vin,0.0f, 0.0f));
  return 0;
}

int generate_masks(const volume4D<float> &mask, volume4D<float> &masknz, const volume4D<float> &vin)
{
  masknz = mask * (1.0f - binarise(vin,0.0f, 0.0f));
  return 0;
}


int fmrib_main_float(int argc, char* argv[]) 
{

  cout.setf(ios::dec); 
  cout.setf(ios::fixed, ios::floatfield); 
  cout.setf(ios::left, ios::adjustfield); 
  cout.precision(6);  

  volume4D<float> vin, vol, mask, masknz;
  read_volume4D(vol,argv[1]);

  float lthr=0, uthr=0;  // these initial values are not used
  if (masks_used) {
    vin = vol;
    generate_masks(mask,masknz,vin,lthr,uthr);
    vol = vin * mask;
  }

  int narg=2;
  string sarg;
 
  while (narg<argc) {
    sarg=argv[narg];

    if (sarg=="-m") {
      if (masks_used) cout <<  vol.mean(mask) << " ";
      else cout << vol.mean() << " ";
    } else if (sarg=="-M") {
      if (masks_used) cout << vol.mean(masknz) << " ";
      else {
	double nzmean=0;
	nzmean = nonzeromean(vol);
	cout << nzmean << " ";
      }
    } else if (sarg=="-X") {
      if (masks_used) {
	cout << vol.mincoordx(mask) << " " << vol.mincoordy(mask) << " " 
	   << vol.mincoordz(mask) << " ";
      } else {
	cout << vol.mincoordx() << " " << vol.mincoordy() << " " 
	   << vol.mincoordz() << " ";
      }
    } else if (sarg=="-x") { 
      if (masks_used) {
	cout << vol.maxcoordx(mask) << " " << vol.maxcoordy(mask) << " " 
	     << vol.maxcoordz(mask) << " ";
      } else {
	cout << vol.maxcoordx() << " " << vol.maxcoordy() << " " 
	     << vol.maxcoordz() << " ";
      }
    } else if (sarg=="-w") {
      if (!masks_used) { 
	if (vin.nvoxels()<1) { vin = vol; }
	masks_used=true;
	generate_masks(mask,masknz,vin,lthr,uthr); 
	vol = vin * mask; 
      }
      int xmin=masknz.maxx(),xmax=masknz.minx(),ymin=masknz.maxy(),ymax=masknz.miny(),zmin=masknz.maxz(),zmax=masknz.minz(),tmin=masknz.maxt(),tmax=masknz.mint();
      
      for(int t=masknz.mint();t<=masknz.maxt();t++) {
	for(int z=masknz.minz();z<=masknz.maxz();z++) {
	  for(int y=masknz.miny();y<=masknz.maxy();y++) {
	    for(int x=masknz.minx();x<=masknz.maxx();x++) {
	      if (masknz(x,y,z,t)>0.5) {
		// if (masknz(x,y,z)>0.5) {
		if (x<xmin) xmin=x;
		if (x>xmax) xmax=x;
		if (y<ymin) ymin=y;
		if (y>ymax) ymax=y;
		if (z<zmin) zmin=z;
		if (z>zmax) zmax=z;
		if (t<tmin) tmin=t;
		if (t>tmax) tmax=t;
	      }
	    }
	  }
	}
      }
      cout << xmin << " " << 1+xmax-xmin << " " << ymin << " " << 1+ymax-ymin << " " << zmin << " " << 1+zmax-zmin << " " << tmin << " " << 1+tmax-tmin << " ";
    } else if (sarg=="-e") {
      if (!masks_used) { 
	if (vin.nvoxels()<1) { vin = vol; }
	masks_used=true;
	generate_masks(mask,masknz,vin,lthr,uthr); 
	vol = vin * mask; 
      }
      ColumnVector hist;
      int nbins=1000;
      double entropy=0;
      hist = vol.histogram(nbins,mask);
      double ntot = hist.Sum();
      for (int j=1; j<=nbins; j++) {
	if (hist(j)>0) {
	  entropy -= (hist(j)/ntot) * log(hist(j)/ntot);	
	}
      }
      entropy /= log((double) nbins);
      cout << entropy << " ";
    } else if (sarg=="-E") { 
      ColumnVector hist;
      int nbins=1000;
      double entropy=0;
      hist = vol.histogram(nbins,masknz);
      double ntot = hist.Sum();
      for (int j=1; j<=nbins; j++) {
	if (hist(j)>0) {
	  entropy -= (hist(j)/ntot) * log(hist(j)/ntot);	
	}
      }
      entropy /= log((double) nbins);
      cout << entropy << " ";
    } else if (sarg=="-k") {
      narg++;
      if (narg>=argc) {
	cerr << "Must specify an argument to -k" << endl;
	exit(2);
      }
      read_volume4D(mask,argv[narg]);
      if (!samesize(mask[0],vol[0])) {
	cerr << "Mask and image must be the same size" << endl;
	exit(3);
      }
      if ( mask.tsize() > vol.tsize() ) {
	cerr << "Mask and image must be the same size" << endl;
	exit(3);
      }
      if ( mask.tsize() != vol.tsize() ) {
	// copy the last max volume until the correct size is reached
	while (mask.tsize() < vol.tsize() ) {
   	  mask.addvolume(mask[mask.maxt()]);
        }
      }
      if (!masks_used) {
	masks_used=true;
	vin = vol;
      }
      float th= 0.5;
      if (th!=0) {
        mask.binarise(th);
      } else {
        mask.binarise(1);
      }
      generate_masks(mask,masknz,vin);
      vol = vin * mask;
    } else if (sarg=="-l") {
      narg++;
      if (narg<argc) {
        lthr = atof(argv[narg]);
      } else {
	cerr << "Must specify an argument to -l" << endl;
	exit(2);
      }
      lthr_used=true;
      if (!masks_used) {
	masks_used=true;
	vin = vol;
      }
      generate_masks(mask,masknz,vin,lthr,uthr);
      vol = vin * mask;
    } else if (sarg=="-u") {
      narg++;
      if (narg<argc) {
        uthr = atof(argv[narg]);
      } else {
	cerr << "Must specify an argument to -u" << endl;
	exit(2);
      }
      uthr_used=true;
      if (!masks_used) {
	masks_used=true;
	vin = vol;
      }
      generate_masks(mask,masknz,vin,lthr,uthr);
      vol = vin * mask;
    } else if (sarg=="-a") {
      vol = abs(vin);
    } else if (sarg=="-v") {
      if (masks_used) {
	cout << (long int) mask.sum() << " " 
	     << mask.sum() * vol.xdim() * vol.ydim() * vol.zdim() << " ";
      } else {
	cout << (long int) vol.nvoxels() << " "
	     << vol.nvoxels() * vol.xdim() * vol.ydim() * vol.zdim() << " ";
      }
    } else if (sarg=="-V") {
      if (masks_used) {
	cout << (long int) masknz.sum() << " " 
	     << masknz.sum() * vol.xdim() * vol.ydim() * vol.zdim() << " ";
      } else {
	long int nzvox;
	nzvox = nonzerocount(vol);
	cout << nzvox << " " 
	     << nzvox * vol.xdim() * vol.ydim() * vol.zdim() << " ";
      }
    } else if (sarg=="-d") {
	// hidden debug option!
      cout << vol.sum() << " ";
    } else if (sarg=="-s") {
	if (masks_used) cout << vol.stddev(mask) << " ";
	else cout << vol.stddev() << " ";
    } else if (sarg=="-S") {
      if (masks_used) {
	cout << vol.stddev(masknz) << " ";
      } else {
	cout << nonzerostddev(vol) << " ";
      }
    } else if (sarg=="-r") {
      if (masks_used) {
	float rmin=vol.robustmin(mask);
	float rmax=vol.robustmax(mask);
	if (rmin>rmax) { 
	  float tmp = rmax;
	  rmax = rmin;
	  rmin = tmp;
	}
	cout << rmin << " " << rmax << " ";
      } else {
	cout << vol.robustmin() << " " << vol.robustmax() << " ";
      }
    } else if (sarg=="-R") {
	if (masks_used) cout << vol.min(mask) << " " << vol.max(mask) << " ";
	else cout << vol.min() << " " << vol.max() << " ";
    } else if (sarg=="-c") {
	ColumnVector cog(4);
	// convert from fsl mm to voxel to sform coord
	cog.SubMatrix(1,3,1,1) = vol[0].cog();
	cog(4) = 1.0;
	if (vol[0].sform_code()!=NIFTI_XFORM_UNKNOWN) {
	  cog = vol[0].sform_mat() * (vol[0].sampling_mat()).i() * cog; 
	}
	cout << cog(1) << " " << cog(2) << " " << cog(3) << " " ;
    } else if (sarg=="-C") {
    ColumnVector cog(4);
	// convert from fsl mm to voxel coord
	cog.SubMatrix(1,3,1,1) = vol[0].cog();
	cog(4) = 1.0;
	cog = (vol[0].sampling_mat()).i() * cog;
	cout << cog(1) << " " << cog(2) << " " << cog(3) << " " ;
    } else if (sarg=="-p") {
      float n;
      narg++;
      if (narg<argc) {
        n = atof(argv[narg]);
      } else {
	cerr << "Must specify an argument to -p" << endl;
	exit(2);
      }
      if ( (n<0) || (n>100) ) {
    	cerr << "Percentile must be between 0 and 100" << endl;
    	exit(1);
      }
      if (masks_used) cout << vol.percentile((float) n/100.0, mask) << " ";
      else cout << vol.percentile((float) n/100.0) << " ";
    } else if (sarg=="-P") { 
      float n;
      narg++;
      if (narg<argc) {
        n = atof(argv[narg]);
      } else {
	cerr << "Must specify an argument to -P" << endl;
	exit(2);
      }
      if ( (n<0) || (n>100) ) {
    	cerr << "Percentile must be between 0 and 100" << endl;
    	exit(1);
      }
      if (!masks_used) {
	if (vin.nvoxels()<1) { vin = vol; }
	masks_used=true;
	generate_masks(mask,masknz,vin,lthr,uthr); 
	vol = vin * mask; 
      }
      cout << vol.percentile((float) n/100.0,masknz) << " ";
    } else if (sarg=="-h") {
      float n;
      narg++;
      if (narg<argc) {
        n = atof(argv[narg]);
      } else {
	cerr << "Must specify the number of bins" << endl;
	exit(2);
      }
      int nbins = (int) n;
      if (nbins<1) {
    	cerr << "Must specify at least 1 bin" << endl;
    	exit(1);
      }
      if (masks_used) {
	cout << vol.histogram(nbins,vol.min(),vol.max(),mask) << " ";
      } else {
	cout << vol.histogram(nbins,vol.min(),vol.max()) << " ";
      }
   } else if (sarg=="-H") {
      float n;
      narg++;
      if (narg<argc) {
        n = atof(argv[narg]);
      } else {
	cerr << "Must specify the number of bins" << endl;
	exit(2);
      }
      int nbins = (int) n;
      if (nbins<1) {
    	cerr << "Must specify at least 1 bin" << endl;
    	exit(1);
      }
      float min=0;
      narg++;
      if (narg<argc) {
        min = atof(argv[narg]);
      } else {
	cerr << "Must specify the histogram minimum intensity" << endl;
	exit(2);
      }
      float max=0;
      narg++;
      if (narg<argc) {
        max = atof(argv[narg]);
      } else {
	cerr << "Must specify the histogram maximum intensity" << endl;
	exit(2);
      }
      if (masks_used) {
	cout << vol.histogram(nbins,min,max,mask) << " ";
      } else {
	cout << vol.histogram(nbins,min,max) << " ";
      }
    } else {
	cerr << "Unrecognised option: " << sarg << endl;
	exit(3);
    }
  
    narg++;
  }

  cout << endl;
  return 0;
}



int main(int argc,char *argv[])
{

  Tracer tr("main");

  string progname=argv[0];
  if (argc<3) { 
    print_usage(progname);
    return 1; 
  }

  return fmrib_main_float(argc,argv);

}

