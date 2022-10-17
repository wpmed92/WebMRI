//     avwsplit++.cc - split 4D files into 3D files for SPM
//     David Flitney and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"
#include "miscmaths/miscmaths.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: avwsplit <inputfile.hdr>" << endl;
  cout << "       avwsplit <inputfile.hdr> [basename]" << endl;
}


template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input_volume;
  volumeinfo vinfo;
  volume<T> output_volume;
  string input_name=string(argv[1]);
  string output_name="vol";
  read_volume4D(input_volume,input_name,vinfo);
  int t;
  if (argc==3) output_name=string(argv[2]);
  for(t=0;t<=input_volume.maxt();t++)
  {
    output_volume=input_volume[t];
    save_volume(output_volume,(output_name+num2str(t,4)),vinfo);
  }
  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  string progname=argv[0];
  if (argc != 3 && argc != 2) 
  { 
    print_usage(progname);
    return 1; 
  }
   
  string iname=string(argv[1]);
  return call_fmrib_main(dtype(iname),argc,argv); 
}


