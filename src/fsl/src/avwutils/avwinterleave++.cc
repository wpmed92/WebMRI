//     avwinterleave++.cc - combine two interleaved frames
//     Steve Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2001-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"

using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: avwinterleave <in1> <in2> <out> [-i]" << endl;
  cout << "-i : reverse slice order" << endl;
}

template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume<T> input_volume1;
  volume<T> input_volume2;
  volumeinfo vinfo;
  int i,j,k,invert=0;

  if (argc == 5 &&!strcmp(argv[4], "-i")) invert=1;

  read_volume(input_volume1,string(argv[1]),vinfo);
  read_volume(input_volume2,string(argv[2]));  
  
  if (input_volume2.xsize() != input_volume1.xsize()  || input_volume2.ysize() != input_volume1.ysize()  || input_volume2.zsize() != input_volume1.zsize() ) 
  {
    cerr << "Error in size-match along non-concatenated dimension" << endl; 
    return 1;
  }
  volume<T> output_volume(input_volume1.xsize(),input_volume1.ysize(),2*input_volume1.zsize());
  output_volume.setdims(input_volume1.xdim(),input_volume1.ydim(),input_volume1.zdim());
  
  for(k=0;k<input_volume1.zsize();k++)
    for(j=0;j<input_volume1.ysize();j++)	    
      for(i=0;i<input_volume1.xsize();i++)
      {
        if (invert)
	{
	  output_volume.value(i,j,2*(input_volume1.zsize()-k)-1)=input_volume2.value(i,j,k);
	  output_volume.value(i,j,2*(input_volume1.zsize()-1-k))=input_volume1.value(i,j,k);
	}
        else
	{
	  output_volume.value(i,j,2*k)=input_volume1.value(i,j,k);
	  output_volume.value(i,j,2*k+1)=input_volume2.value(i,j,k);
	}
      }
  FslSetCalMinMax(&vinfo,output_volume.min(),output_volume.max());
  save_volume(output_volume,string(argv[3]),vinfo);
  return 0;
}


int main(int argc,char *argv[])
{
  if (argc < 3) 
  {
    print_usage(string(argv[0]));
    return 1; 
  }
  return call_fmrib_main(dtype(string(argv[1])),argc,argv); 
}


