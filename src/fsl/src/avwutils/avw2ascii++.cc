//     avw2ascii++.c - convert AVW to raw ASCII text
//     Stephen Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2001-2005 University of Oxford  
//     COPYRIGHT  
#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"
#include "miscmaths/miscmaths.h"
#include <fstream>

using namespace NEWIMAGE;
using namespace MISCMATHS;

void print_usage(const string& progname) {
  cout << endl;
  cout << "Usage: avw2ascii <input> <output>" << endl;
}


template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input_volume;
  string input_name=string(argv[1]);
  read_volume4D(input_volume,input_name);
  string output_name=string(argv[2]);
  ofstream output_file;
  int i,j,k,t;
  for(t=0;t<=input_volume.maxt();t++)
  {
    output_file.open((output_name+num2str(t,5)).c_str(),ofstream::out);
    for(k=0;k<=input_volume.maxz();k++)
    {
      for(j=0;j<=input_volume.maxy();j++)
      {
        for(i=0;i<=input_volume.maxx();i++)
	{
          output_file << input_volume(i,j,k,t) << " ";
        }
        output_file << endl;
      }
      output_file << endl;
    }
    output_file.close();
  } 
  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  string progname=argv[0];
  if (argc != 3) 
  { 
    print_usage(progname);
    return 1; 
  }
   
  string iname=string(argv[1]);
  return call_fmrib_main(dtype(iname),argc,argv); 
}


