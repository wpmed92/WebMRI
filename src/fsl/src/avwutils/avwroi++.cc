//     avwroi++.cc  extract cuboid ROI and/or timeseries from image
//     Stephen Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 1999-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"

using namespace NEWIMAGE;

void print_usage(const string& progname) {
  cout << endl;
  cout << "Usage: avwroi <input> <output> <xmin> <xsize> <ymin> <ysize> <zmin> <zsize>" << endl;
  cout << "       avwroi <input> <output> <tmin> <tsize>\n" << endl;
  cout << "       avwroi <input> <output> <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize>" << endl;
  cout << "Note: counting (in both time and space) starts with 0 not 1!" << endl;
}


template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input_vol;
  volumeinfo vinfo;
  volume4D<T> output_vol;
  string input_name=string(argv[1]);
  string output_name=string(argv[2]);
  read_volume4D(input_vol,input_name,vinfo);

  int minx = -1,miny = -1,minz = -1,mint = -1,maxx = -1,maxy = -1,maxz= -1,maxt = -1,argindex = 3;
  //Note there is currently no explicit bounds testing e.g. if maxx>inputvol.maxx() etc
  if (argc==5) //4D Timeseries only
  {
    minx=0;
    maxx=input_vol.maxx();
    miny=0;
    maxy=input_vol.maxy();
    minz=0;
    maxz=input_vol.maxz();
    mint=atoi(argv[argindex++]);
    maxt=atoi(argv[argindex++]);
    maxt+=mint-1;
  }
  else if (argc==9)  //3D Region of interest 
  {
    minx=atoi(argv[argindex++]);   //N.B. could compact some of these lines with the ones below...
    maxx=atoi(argv[argindex++]);
    miny=atoi(argv[argindex++]);
    maxy=atoi(argv[argindex++]);
    minz=atoi(argv[argindex++]);
    maxz=atoi(argv[argindex++]);
    mint=0;
    maxt=input_vol.maxt();
    maxx+=minx-1;
    maxy+=miny-1;
    maxz+=minz-1;
  }
  else if (argc==11)   //4D Timeseries and Region of interest
  {
    minx=atoi(argv[argindex++]);
    maxx=atoi(argv[argindex++]);
    miny=atoi(argv[argindex++]);
    maxy=atoi(argv[argindex++]);
    minz=atoi(argv[argindex++]);
    maxz=atoi(argv[argindex++]);
    mint=atoi(argv[argindex++]);
    maxt=atoi(argv[argindex++]); 
    maxx+=minx-1;
    maxy+=miny-1;
    maxz+=minz-1;
    maxt+=mint-1;
  }
  input_vol.setROIlimits(minx,miny,minz,mint,maxx,maxy,maxz,maxt);
  input_vol.activateROI();  //is this needed?
  output_vol=input_vol.ROI();
  save_volume4D(output_vol,output_name,vinfo);
  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  string progname=argv[0];
  if (argc !=5 && argc !=9 && argc !=11) 
  { 
    print_usage(progname);
    return 1; 
  }
   
  string iname=string(argv[1]);
  return call_fmrib_main(dtype(iname),argc,argv); 
}


