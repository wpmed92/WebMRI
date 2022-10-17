//     avwmaths++.cc Image processing routines, some basic, some not so basic...
//     Steve Smith, David Flitney, Stuart Clare and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2006 University of Oxford  
//     CCOPYRIGHT  

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"


#if defined ( __CYGWIN__ ) ||  defined (__sun)
extern "C" { 
#include <ieeefp.h> 
}
#endif

using namespace MISCMATHS;
using namespace NEWIMAGE;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

void print_usage(const string& progname) 
{
  cout << "Usage: avwmaths [-dt <datatype>] <first_input> [operations and inputs] <output> [-odt <datatype>]" << endl;

  cout << "\nDatatype information:" << endl;
  cout << " -dt sets the native (math operations) datatype (default float for all except double images)" << endl;
  cout << " -odt sets the output datatype (default as original image)" << endl;
  cout << " Possible datatypes are: char short int float double" << endl;
  cout << " Additionally \"-dt input\" will set the native datatype to that of the original image" << endl;

  cout << "\nBinary operations:" << endl;
  cout << "  (some inputs can be either an image or a number)" << endl;
  cout << " -add   : add following input to current image" << endl;
  cout << " -sub   : subtract following input from current image" << endl;
  cout << " -mul   : multiply current image by following input" << endl;
  cout << " -div   : divide current image by following input" << endl;
  cout << " -mas   : use (following image>0) to mask current image" << endl;
  cout << " -thr   : use following number to threshold current image (zero anything below the number)" << endl;
  cout << " -thrp  : use following percentage (0-100) of ROBUST RANGE to threshold current image (zero anything below the number)" << endl;
  cout << " -thrP  : use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold below" << endl;
  cout << " -uthr  : use following number to upper-threshold current image (zero anything above the number)" << endl;
  cout << " -uthrp : use following percentage (0-100) of ROBUST RANGE to upper-threshold current image (zero anything above the number)" << endl;
  cout << " -uthrP : use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold above" << endl;
  cout << " -max   : take maximum of following input and current image" << endl;
  cout << " -min   : take minimum of following input and current image" << endl;

  cout << "\nBasic unary operations:" << endl;
  cout << " -exp   : exponential" << endl;
  cout << " -log   : natural logarithm" << endl;
  cout << " -sqr   : square" << endl;
  cout << " -sqrt  : square root" << endl;
  cout << " -abs   : absolute value" << endl;
  cout << " -bin   : use (current image>0) to binarise" << endl;
  cout << " -index : replace each nonzero voxel with a unique (subject to wrapping) index number" << endl;
  cout << " -grid <value> <spacing> : add a 3D grid of intensity <value> with grid spacing <spacing>" << endl;
  cout << " -edge  : edge strength" << endl;
  cout << " -cluster <height_power> <size_power> <connectivity>: cluster-enhance (connectivity = 6 for 3D, 26 for skeletons)" << endl;
  //cout << " -clusterhist <connectivity>: cluster-enhance-histogram (connectivity = 6 for 3D, 26 for skeletons)" << endl;
  cout << " -nan   : replace NaNs (improper numbers) with 0" << endl;
  cout << " -nanm  : make NaN (improper number) mask with 1 for NaN voxels, 0 otherwise" << endl;
  cout << " -inm <mean> :  (-i i ip.c) intensity normalisation (per 3D volume mean)" << endl;
  cout << " -ing <mean> :  (-I i ip.c) intensity normalisation, global 4D mean)" << endl;

  cout << "\nKernel operations:" << endl;
  cout << " -kernel 3D : 3x3x3 box centered on target voxel (set as default kernel)" << endl;
  cout << " -kernel 2D : 3x3x1 box centered on target voxel" << endl;
  cout << " -kernel box    <size>     : all voxels in a box of width <size> centered on target voxel" << endl;
  cout << " -kernel boxv   <size>     : <size>x<size>x<size> box centered on target voxel, CAUTION: size should be an odd number" << endl;
  cout << " -kernel gauss  <sigma>    : gaussian kernel (sigma in mm, not voxels)" << endl;
  cout << " -kernel sphere <size>     : all voxels in a sphere of radius <size> mm centered on target voxel" << endl;
  cout << " -kernel file   <filename> : use external file as kernel" << endl;

  cout << "\nSpatial Filtering operations: N.B. all these options use the kernel specified by -kernel" << endl;
  cout << " -dilM    : Mean Dilation of zero voxels  (using non-zero voxels in kernel)" << endl;
  cout << " -dilD    : Modal Dilation of zero voxels (using non-zero voxels in kernel)" << endl;
  cout << " -dilF    : Maximum filtering of all voxels" << endl;
  cout << " -ero     : Erode by zeroing non-zero voxels when zero voxels found in kernel" << endl;
  cout << " -eroF    : Minimum filtering of all voxels" << endl;
  cout << " -fmedian : Median Filtering " << endl;
  cout << " -fmean   : Mean filtering, kernel weighted (conventionally used with gauss kernel)" << endl;
  cout << " -fmeanu  : Mean filtering, kernel weighted, un-normalised (gives edge effects)" << endl;

  cout << "\nDimensionality reduction operations:" << endl;
  cout << "  (the \"T\" can be replaced by X, Y or Z to collapse across a different dimension)" << endl;
  cout << " -Tmean   : mean across time" << endl;
  cout << " -Tstd    : standard deviation across time" << endl;
  cout << " -Tmax    : max across time" << endl;
  cout << " -Tmaxn   : time index of max across time" << endl;
  cout << " -Tmin    : min across time" << endl;
  cout << " -Tmedian : median across time" << endl;
  cout << " -Tperc <percentage> : nth percentile (0-100) of FULL RANGE across time" << endl;
  cout << " -Tar1    : temporal AR(1) coefficient (use -odt float and probably demean first)" << endl;

  cout << "\nMulti-argument operations:" << endl;
  cout << " -roi <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> : zero outside roi" << endl;
  cout << " -bptf  <hp_sigma> <lp_sigma> : (-t in ip.c) Bandpass temporal filtering; nonlinear highpass and Gaussian linear lowpass (with sigmas in volumes, not seconds); set either sigma<0 to skip that filter" << endl;

  cout << "\ne.g. avwmaths input_volume -add input_volume2 output_volume" << endl;
  cout << "     avwmaths input_volume -add 2.5 output_volume" << endl;
  cout << "     avwmaths input_volume -add 2.5 -mul input_volume2 output_volume\n" << endl;
}

bool isNumber( const string& x )
{
bool flag=true;
 if (!(isdigit(x.at(0)) || x.at( 0 ) == '-' || x.at( 0 ) == '.')) flag=false;
 for (unsigned int i = 1; i < x.size(); i++ )
   if (!(isdigit(x.at(i)) || x.at( i ) == '.')) flag=false;
 return flag; 
} 

template <class T>
int fmrib_main(int argc, char *argv[], short output_dt)
{
  volume4D<T> input_volume;
  volumeinfo vinfo;
  //for (int i = 2; i < argc-1; i++) cout << argv[i] << endl;
  volume<float> kernel;
  kernel=box_kernel(3,3,3);
  bool seperable=false;
  read_volume4D(input_volume,string(argv[1]),vinfo);

  for (int i = 2; i < argc-1; i++)  //main loop
  {    
    volume4D<T> temp_volume;
    /********************Dimensionality Reduction*******************/
    /***************************************************************/
    if (isupper((int)argv[i][1]) && argv[i][0] == '-')  //if first letters are -capital - dimensionality reduction...
    { 
      int xoff=1,yoff=1,zoff=1,toff=1,nsize;
      if (argv[i][1] == 'T') toff=input_volume.tsize(); 
      if (argv[i][1] == 'Z') zoff=input_volume.zsize();  
      if (argv[i][1] == 'Y') yoff=input_volume.ysize();  
      if (argv[i][1] == 'X') xoff=input_volume.xsize();  
      temp_volume=input_volume;
      input_volume.reinitialize(input_volume.xsize()/xoff,input_volume.ysize()/yoff,input_volume.zsize()/zoff,input_volume.tsize()/toff); 
      input_volume.copyproperties(temp_volume);
      nsize=xoff*yoff*zoff*toff;
      volume<T> column_volume(nsize,1,1); //will be size of appropriate dimension, as only 1 arg is non-unitary
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	    {
              for (int j=0;j<nsize;j++) column_volume.value(j,0,0)=temp_volume(x+j*(xoff!=1),y+j*(yoff!=1),z+j*(zoff!=1),t+j*(toff!=1));
              //This goes along the appropriate axis (non unitary offset variable) and fills a "column" volume with data
	      if (string(argv[i]+2) == "max")    input_volume.value(x,y,z,t)=column_volume.max(); 
	      if (string(argv[i]+2) == "min")    input_volume.value(x,y,z,t)=column_volume.min();
              if (string(argv[i]+2) == "mean")   input_volume.value(x,y,z,t)=(T)column_volume.mean();
              if (string(argv[i]+2) == "std")    input_volume.value(x,y,z,t)=(T)column_volume.stddev();
	      if (string(argv[i]+2) == "maxn")   input_volume.value(x,y,z,t)=column_volume.maxcoordx();
              if (string(argv[i]+2) == "median") input_volume.value(x,y,z,t)=column_volume.percentile(0.5);
	      if (string(argv[i]+2) == "perc")   input_volume.value(x,y,z,t)=column_volume.percentile(atof(argv[i+1])/100.0);
              if (string(argv[i]+2) == "ar1") 
              {
                column_volume-=(T)column_volume.mean();
                double sumsq=column_volume.sumsquares();
                input_volume(x,y,z,t)=0;
		if(sumsq!=0) for (int k=1;k<nsize;k++) input_volume(x,y,z,t)+=(T)(column_volume(k,0,0)*column_volume(k-1,0,0)/sumsq);
	      }
	    }
       if (string(argv[i]+2) == "perc") i++;
       
    }
    /********************Binary Operations**************************/
    /***************************************************************/
    else if (string(argv[i])=="-mas")
    {  
        read_volume4D(temp_volume,string(argv[++i]));
        temp_volume.binarise(0,temp_volume.max()+1,exclusive); // needed to binarise max() value + 1 due to
        for (int t=0;t<input_volume.tsize();t++)                     //potential issue with exclusive binarise
          input_volume[t]*=temp_volume[t%temp_volume.tsize()]; //this gives compatibility with 3 and 4D masks
    }                                                            //without needing to have a specific volume3D variable
    /***************************************************************/
    else if (string(argv[i])=="-thr") input_volume.threshold((T)atof(argv[++i]),input_volume.max()+1,inclusive);
    /***************************************************************/
    else if (string(argv[i])=="-thrp") 
    {
      T lowerlimit =(T)(input_volume.robustmin()+(atof(argv[++i])/100.0)*(input_volume.robustmax()-input_volume.robustmin())); 
      input_volume.threshold(lowerlimit,input_volume.max()+1,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-thrP") 
    {
      volume4D<T> mask(input_volume);
      mask.binarise(0,input_volume.max()+1,exclusive);
      T lowerlimit =(T)(input_volume.robustmin(mask)+(atof(argv[++i])/100.0)*(input_volume.robustmax(mask)-input_volume.robustmin(mask))); 
      input_volume.threshold(lowerlimit,input_volume.max()+1,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-uthr") input_volume.threshold(input_volume.min()-1,(T)atof(argv[++i]),inclusive);
    /***************************************************************/
    else if (string(argv[i])=="-uthrp") 
    {
      T upperlimit = (T)(input_volume.robustmin()+(atof(argv[++i])/100.0)*(input_volume.robustmax()-input_volume.robustmin())); 
       input_volume.threshold(input_volume.min()-1,upperlimit,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-uthrP") 
    {
       volume4D<T> mask(input_volume);
       mask.binarise(0,input_volume.max()+1,exclusive);
       T upperlimit = (T)(input_volume.robustmin(mask)+(atof(argv[++i])/100.0)*(input_volume.robustmax(mask)-input_volume.robustmin(mask))); 
       input_volume.threshold(input_volume.min()-1,upperlimit,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-kernel") 
    {
       kernel.destroy();
       float xdim=input_volume.xdim();
       float ydim=input_volume.ydim();
       float zdim=input_volume.zdim();
       if(string(argv[i+1])=="-2D")      kernel=box_kernel(3,3,1);
       else if(string(argv[i+1])=="-3D") kernel=box_kernel(3,3,3);
       else
       {
	 float size=atof(argv[i+2]);
	 if(string(argv[i+1])=="box")      kernel=box_kernel(size,xdim,ydim,zdim);
	 if(string(argv[i+1])=="boxv")   kernel=box_kernel((int)size,(int)size,(int)size);
         else if(string(argv[i+1])=="gauss")  kernel=gaussian_kernel3D(size,xdim,ydim,zdim);
         else if(string(argv[i+1])=="sphere") kernel=spherical_kernel(size,xdim,ydim,zdim);
	 else if(string(argv[i+1])=="file")   read_volume(kernel,string(argv[i+2]));
         if(string(argv[i+1])=="box" || string(argv[i+1])=="gauss") seperable=true;       
         else seperable=false;  
	  i++;
       }
       i++;
       save_volume(kernel,"kernel");
       
    }
    /***************************************************************/
    else if (string(argv[i])=="-add"){
      if (isNumber(string(argv[++i]))) input_volume+=(T)atof(argv[i]); 
      else if (FslFileExists(argv[i])) 
      {  
	read_volume4D(temp_volume,string(argv[i]));
        for (int t=0;t<input_volume.tsize();t++) input_volume[t]+=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-sub"){
      if (isNumber(string(argv[++i]))) input_volume-=(T)atof(argv[i]);
      else if (FslFileExists(argv[i])) 
      {  
	read_volume4D(temp_volume,string(argv[i]));
        for (int t=0;t<input_volume.tsize();t++) input_volume[t]-=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-mul"){
      if (isNumber(string(argv[++i]))) input_volume*=(T)atof(argv[i]);
      else if (FslFileExists(argv[i])) 
      {  
	read_volume4D(temp_volume,string(argv[i]));
        for (int t=0;t<input_volume.tsize();t++) input_volume[t]*=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-div"){
      if (isNumber(string(argv[++i]))) {if (atof(argv[i])!=0) input_volume/=(T)atof(argv[i]);}
      else if (FslFileExists(argv[i])) 
      {  
        read_volume4D(temp_volume,string(argv[i]));
        for(int t=0;t<input_volume.tsize();t++)      
	{
          int t2=t%temp_volume.tsize();     
          for(int z=0;z<input_volume.zsize();z++)
            for(int y=0;y<input_volume.ysize();y++)	    
	      for(int x=0;x<input_volume.xsize();x++)
		  if(temp_volume(x,y,z,t2)!=0) input_volume.value(x,y,z,t) /= temp_volume.value(x,y,z,t2);
                  else input_volume.value(x,y,z,t)=(T)0.0;
        }          
      }}
    /***************************************************************/
    else if (string(argv[i])=="-max" || string(argv[i])=="-min")
    {
      T param=0;
      bool max=false;
      bool file=false;
      if (string(argv[i])=="-max") max=true;
      if (isNumber(string(argv[++i]))) param=(T)atof(argv[i]);
      else if (FslFileExists(argv[i])) 
      {                           
	read_volume4D(temp_volume,string(argv[i]));
        file=true;
      }     
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	    {
              if (max && file) input_volume.value(x,y,z,t)=MAX(input_volume.value(x,y,z,t),temp_volume.value(x,y,z,t%temp_volume.tsize()));
              if (max && !file) input_volume.value(x,y,z,t)=MAX(input_volume.value(x,y,z,t),param);
              if (!max && file) input_volume.value(x,y,z,t)=MIN(input_volume.value(x,y,z,t),temp_volume.value(x,y,z,t%temp_volume.tsize()));
              if (!max && !file) input_volume.value(x,y,z,t)=MIN(input_volume.value(x,y,z,t),param);
            }
    }
    /*********************Unary Operations**************************/
    /***************************************************************/
    else if (string(argv[i])=="-sqrt")
    {
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	    {
              if (input_volume.value(x,y,z,t)> 0) input_volume.value(x,y,z,t)=(T)sqrt(input_volume.value(x,y,z,t));
              else  input_volume.value(x,y,z,t)= 0; 
            }
    }
    /***************************************************************/
    else if (string(argv[i])=="-sqr") input_volume*=input_volume;
    /***************************************************************/
    else if (string(argv[i])=="-exp")
    {
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
               input_volume.value(x,y,z,t)=(T)exp((double)input_volume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-log")
    {
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
              if (input_volume.value(x,y,z,t)> 0) input_volume.value(x,y,z,t)=(T)log((double)input_volume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-abs") input_volume=abs(input_volume);
    /***************************************************************/
    else if (string(argv[i])=="-bin") input_volume.binarise(0,input_volume.max()+1,exclusive); 
    /***************************************************************/
    else if (string(argv[i])=="-index")
    {
      int indexval=0;
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	      if (input_volume.value(x,y,z,t)>0) 
		{
		  input_volume.value(x,y,z,t)=(T)indexval;
		  indexval++;
		}
    }
    /***************************************************************/
    else if (string(argv[i])=="-grid")
    {
      double gridvalue = atof(argv[++i]);
      int gridspacing = atoi(argv[++i]);
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	      if ( x%gridspacing==0 || y%gridspacing==0 || z%gridspacing==0 )
		input_volume.value(x,y,z,t)=(T)gridvalue;
    }
    /*****************SPATIAL FILTERING OPTIONS*********************/
    /***********************Mean Dilation***************************/
    else if (string(argv[i])=="-dilM")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"dilateM");
    /***********************Modal Dilation**************************/
    else if (string(argv[i])=="-dilD")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"dilateD");
    /***********************MJ Dilation**************************/
    else if (string(argv[i])=="-dilF")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"dilate");
    /***********************Steves Erosion**************************/
    else if (string(argv[i])=="-ero")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"erodeS");
    /**************************MJ Erosion**************************/
    else if (string(argv[i])=="-eroF")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"erode");
    /***********************Median Filtering***********************/
    else if (string(argv[i])=="-fmedian")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"median");
    /******************Mean Filtering*************************/
    else if (string(argv[i])=="-fmean")
       input_volume=generic_convolve(input_volume,kernel,seperable,true);
    /******************Mean Filtering Unnormalised************/
    else if (string(argv[i])=="-fmeanu")
       input_volume=generic_convolve(input_volume,kernel,false,false);
    /*****************END OF FILTERING OPTIONS***************/
    else if (string(argv[i])=="-edge")
       input_volume=edge_strengthen(input_volume);
    else if (string(argv[i])=="-cluster")
      {
	float height_power = atof(argv[++i]);
	float size_power = atof(argv[++i]);
	int connectivity = atoi(argv[++i]);

	for(int t=0;t<input_volume.tsize();t++)           
	  {
	    float maxval=input_volume[t].max();
	    volume<float> clusterenhance;
	    copyconvert(input_volume[0],clusterenhance);
	    clusterenhance=0;

	    float delta=maxval/100;  // i.e. 100 subdivisions of the max input stat height
	    for (float thresh=delta; thresh<maxval; thresh+=delta)
	      {
		volume<float> clusters;
		copyconvert(input_volume[t],clusters);
		clamp(clusters,clusters.min()-1,thresh+delta);
		clusters.binarise(thresh);

		ColumnVector clustersizes;  
		volume<int>tmpvol=connected_components(clusters,clustersizes,connectivity);
		clustersizes = pow(clustersizes,size_power) * pow(thresh,height_power);
		for(int z=0;z<input_volume.zsize();z++)
		  for(int y=0;y<input_volume.ysize();y++)	    
		    for(int x=0;x<input_volume.xsize();x++)
		      if (tmpvol(x,y,z)>0)
			clusterenhance(x,y,z) += clustersizes(tmpvol(x,y,z));
	      }
	    //clusterenhance = (clusterenhance/clusterenhance.max())*maxval;
	    copyconvert(clusterenhance,input_volume[t]);
	  }
      }
    else if (string(argv[i])=="-clusterhist")
      {
	int connectivity = atoi(argv[++i]);
	int maxvol=50;
	float maxval=input_volume.max();
	cerr << "maxval=" << maxval << endl;
	float delta=maxval/100;  // i.e. 100 subdivisions of the max input stat height
	Matrix hist(maxvol , 101);

	for(int t=0;t<input_volume.tsize();t++)           
	  {
	    cerr << "t=" << t << endl;
	    for (float h=1; h<100; h++)
	      {
		float thresh=h*delta;
		volume<float> clusters;
		copyconvert(input_volume[t],clusters);
		clamp(clusters,clusters.min()-1,thresh+delta);
		clusters.binarise(thresh);

		ColumnVector clustersizes;  
		connected_components(clusters,clustersizes,connectivity);
		for(int r=1; r<=clustersizes.Nrows(); r++)
		  hist(Min(clustersizes(r),maxvol),h)++;
	      }
	  }

	cout << hist.Rows(1,maxvol) << endl;
      }
    /******************************************************/
    else if (string(argv[i])=="-nanm")
     {   
       for(int t=0;t<input_volume.tsize();t++)           
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               if ( finite((double)input_volume.value(x,y,z,t))) input_volume.value(x,y,z,t)=0;
	       else input_volume.value(x,y,z,t)=1;
     }
     /******************************************************/
    else if (string(argv[i])=="-nan")
     {   
       for(int t=0;t<input_volume.tsize();t++)           
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               if (!finite((double)input_volume.value(x,y,z,t))) input_volume.value(x,y,z,t)=0;     
     }
     /******************************************************/
    else if (string(argv[i])=="-roi")
     {   
       for(int t=0;t<input_volume.tsize();t++)           
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               if((x<atoi(argv[i+1])) || (x>=atoi(argv[i+1])+atoi(argv[i+2])) || (y<atoi(argv[i+3])) || (y>=atoi(argv[i+3])+atoi(argv[i+4])) || (z<atoi(argv[i+5])) || (z>=atoi(argv[i+5])+atoi(argv[i+6])) || (t<atoi(argv[i+7])) || (t>=atoi(argv[i+7])+atoi(argv[i+8])) )
                 input_volume.value(x,y,z,t)=0;
       i+=8;
     }
     /*******************IP functions***********************/
    else if (string(argv[i])=="-inm")
     { 
       double target,tmean;
       target = atof(argv[++i]);
       volume4D<T> mask(input_volume);   
       mask.binarise(0,mask.max()+1,exclusive); 
       for(int t=0;t<input_volume.tsize();t++)     
       {
         tmean=target/input_volume[t].mean(mask[t]);
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               input_volume.value(x,y,z,t)=(T)(input_volume.value(x,y,z,t)*tmean);
       }
     }
    else if (string(argv[i])=="-ing")
     { 
       double tmean,target;
       target = atof(argv[++i]);
       volume4D<T> mask(input_volume);   
       mask.binarise(0,mask.max()+1,exclusive); 
       tmean=target/input_volume.mean(mask);
       for(int t=0;t<input_volume.tsize();t++)     
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               input_volume.value(x,y,z,t)=(T)(input_volume.value(x,y,z,t)*tmean);
     }
    else if(string(argv[i])=="-bptf")
     {
       input_volume=bandpass_temporal_filter(input_volume,atof(argv[i+1]),atof(argv[i+2]));
       i+=2;
     }
    else { cout << "\n Error in command line: unknown option \"" << argv[i] << "\"\n" << endl; print_usage("blah"); return 1; }
     /******************************************************/
  } 

  if (dtype(input_volume)>=DT_FLOAT && output_dt < DT_FLOAT)
  {
    for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
              input_volume.value(x,y,z,t)=(T) MISCMATHS::round(input_volume.value(x,y,z,t));
  }

  FslSetCalMinMax(&vinfo,input_volume.min(),input_volume.max());
  save_volume4D_dtype(input_volume,string(argv[argc-1]),output_dt,vinfo,true);
  return 0;
}


int main(int argc,char *argv[])
{
  if (argc < 2) 
  { 
    print_usage(string(argv[0]));
    return 1; 
  }
  short original_dt;
  if(string(argv[1]) =="-datatype" || string(argv[1])== "-dt")  original_dt = dtype(string(argv[3]));
  else original_dt = dtype(string(argv[1]));
  short output_dt=original_dt;
  if(string(argv[argc-2])=="-output_datatype" || string(argv[argc-2])== "-odt") //output datatype
  {
    if(string(argv[argc-1])=="char")        output_dt =  DT_UNSIGNED_CHAR;  
    else if(string(argv[argc-1])=="short")  output_dt =  DT_SIGNED_SHORT;
    else if(string(argv[argc-1])=="int")    output_dt =  DT_SIGNED_INT;
    else if(string(argv[argc-1])=="float")  output_dt =  DT_FLOAT;
    else if(string(argv[argc-1])=="double") output_dt =  DT_DOUBLE;
    else {cout << "Error: Unknown datatype \"" << argv[argc-1] << "\" - Possible datatypes are: char short int float double" << endl; return 1;}
    argc-=2;
  }
  if(string(argv[1])=="-datatype" || string(argv[1])== "-dt") //input datatype
  {
    short input_dt=-1;
     if(string(argv[2])=="input") input_dt=original_dt;
     else if(string(argv[2])=="char" || input_dt == DT_UNSIGNED_CHAR)     return fmrib_main<char>(argc-2, argv+2,output_dt);
     else if(string(argv[2])=="short" || input_dt == DT_SIGNED_SHORT)return fmrib_main<short>(argc-2, argv+2,output_dt);
     else if(string(argv[2])=="int" || input_dt == DT_SIGNED_INT)    return fmrib_main<int>(argc-2, argv+2,output_dt);
     else if(string(argv[2])=="float" || input_dt == DT_FLOAT)   return fmrib_main<float>(argc-2, argv+2,output_dt); 
     else if(string(argv[2])=="double" || input_dt == DT_DOUBLE) return fmrib_main<double>(argc-2, argv+2,output_dt); 
     else {cout << "Error: Unknown datatype \"" << argv[2] <<  "\" - Possible datatypes are: char short int float double input" << endl; return 1;}
  }
  else if (dtype(string(argv[1]))==DT_DOUBLE) return fmrib_main<double>(argc,argv,output_dt);
  else return fmrib_main<float>(argc,argv,output_dt);
}



