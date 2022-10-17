//     avwmerge++.cc concatenate AVW files into a single output
//     Steve Smith, David Flitney, Stuart Clare and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"

using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: avwmerge <-x/y/z/t/a> <output> <file1 file2 .......>" << endl;
  cout << "     -t : concatenate images in time" << endl;
  cout << "     -x : concatenate images in the x direction"  << endl;
  cout << "     -y : concatenate images in the y direction"  << endl;
  cout << "     -z : concatenate images in the z direction" << endl;
  cout << "     -a : auto-choose: single slices -> volume, volumes -> 4D (time series)"  << endl;
}

template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input_volume;
  volumeinfo vinfo;
  int i,j,k,t,vol,direction,dimerror=0,xdimtot=0,ydimtot=0,zdimtot=0,tdimtot=0,xoffset=0,yoffset=0,zoffset=0,toffset=0;
 
  if (!strcmp(argv[1], "-t"))       direction=0;
  else if (!strcmp(argv[1], "-x"))  direction=1;
  else if (!strcmp(argv[1], "-y"))  direction=2;
  else if (!strcmp(argv[1], "-z"))  direction=3; 
  else if (!strcmp(argv[1], "-a"))  direction=4;
  else 
  {
    print_usage(string(argv[0]));
    return(1);
  }
  read_volume4D_hdr_only(input_volume,string(argv[3]),vinfo);
  xdimtot=input_volume.xsize(); 
  ydimtot=input_volume.ysize(); 
  zdimtot=input_volume.zsize();
  tdimtot=input_volume.tsize();
  if(direction==4)
  {
    if( (zdimtot<2) && (tdimtot<2) ) direction=3;
     else direction=0;
  }
  input_volume.destroy();    //Remove when new newimage comes out

  for(vol = 4; vol < argc; vol++)
  {        
    read_volume4D_hdr_only(input_volume,string(argv[vol]));
    if (direction==0) tdimtot+=input_volume.tsize(); 
  //    cerr << tdimtot << endl; //This will give errors in tdimtot if input_volume.destroy not present
    if (direction==1) xdimtot+=input_volume.xsize(); 
    if (direction==2) ydimtot+=input_volume.ysize();
    if (direction==3) zdimtot+=input_volume.zsize();
    input_volume.destroy();     //Remove when new newimage comes out
  }
  volume4D<T> output_volume(xdimtot,ydimtot,zdimtot,tdimtot);
  read_volume4D(input_volume,string(argv[3]));  
  output_volume.copyproperties(input_volume);

  for(vol = 3; vol < argc; vol++)
  {   
    if (vol>3) read_volume4D(input_volume,string(argv[vol]));  
    if (direction == 0 && (input_volume.xsize() != xdimtot || input_volume.ysize() != ydimtot || input_volume.zsize() != zdimtot)) dimerror=1;
    if (direction == 1 && (input_volume.ysize() != ydimtot || input_volume.zsize() != zdimtot || input_volume.tsize() != tdimtot)) dimerror=1;
    if (direction == 2 && (input_volume.xsize() != xdimtot || input_volume.zsize() != zdimtot || input_volume.tsize() != tdimtot)) dimerror=1; 
    if (direction == 3 && (input_volume.xsize() != xdimtot || input_volume.ysize() != ydimtot || input_volume.tsize() != tdimtot)) dimerror=1;
    if (dimerror)
    {
      cerr << "Error in size-match along non-concatenated dimension" << endl; 
      return 1;
    }

             
    for(t=0;t<input_volume.tsize();t++)           
      for(k=0;k<input_volume.zsize();k++)
        for(j=0;j<input_volume.ysize();j++)	    
          for(i=0;i<input_volume.xsize();i++)
            output_volume.value(i+xoffset,j+yoffset,k+zoffset,t+toffset)=input_volume.value(i,j,k,t);
    if (direction==0)  toffset+=input_volume.tsize();  
    if (direction==1)  xoffset+=input_volume.xsize();  
    if (direction==2)  yoffset+=input_volume.ysize(); 
    if (direction==3)  zoffset+=input_volume.zsize(); 
    input_volume.destroy();   //Remove when new newimage comes out      
  }

  save_volume4D(output_volume,string(argv[2]),vinfo);
  return 0;
}


int main(int argc,char *argv[])
{
  if (argc < 4) 
  { 
    print_usage(string(argv[0]));
    return 1; 
  }
  return call_fmrib_main(dtype(string(argv[3])),argc,argv); 
}


