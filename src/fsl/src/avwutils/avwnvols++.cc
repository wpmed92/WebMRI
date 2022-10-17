//     avwnvols++.cc   REALLY advanced program for counting the volumes in a 4D file
//     Stephen Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 1999-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"

using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: avwnvols <input>" << endl;
}


template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> inputvol;
  string input_name=string(argv[1]);
  read_volume4D_hdr_only(inputvol,input_name);
  cout << inputvol.tsize() << endl;
  return 0;
}



int main(int argc,char *argv[])
{
  string progname=argv[0];
  if (argc != 2)
  {
    print_usage(progname);
    return 1;
  }
  if (!FslFileExists(argv[1])) 
  { 
    cout << "0" << endl;
    return 1; 
  }
  else 
  { 
    string iname=string(argv[1]);
    return call_fmrib_main(dtype(iname),argc,argv); 
  }
}

