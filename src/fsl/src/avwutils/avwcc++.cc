//     avwmcc++.cc cross-correlate two time series
//     Steve Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2001-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"
#include <iomanip>

using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "\nUsage: avwcc <first_input> <second_input> [cc_thresh]" << endl;
}

template <class T>
int fmrib_main(int argc, char *argv[])
{
  volume4D<T> input_volume1;
  volume4D<T> input_volume2;
  string input_name1=string(argv[1]);
  string input_name2=string(argv[2]);
  read_volume4D(input_volume1,input_name1);
  read_volume4D(input_volume2,input_name2);
  int t1,t2,i,j,k;
  double ss1,ss2;
  double thresh=0.1;
  double score;

  if (input_volume1.maxx() != input_volume2.maxx() ||  input_volume1.maxy() != input_volume2.maxy()  ||  input_volume1.maxz() != input_volume2.maxz())
  {
    cerr << "Error: Mismatch in image dimensions" << endl; 
    return 1;
  }
  if (argc==4)  thresh=atof(argv[3]);

  for(t1=0;t1<=input_volume1.maxt();t1++)
  {
    ss1=sqrt(input_volume1[t1].sumsquares());  
    for(t2=0;t2<=input_volume2.maxt();t2++)
    {
       ss2=sqrt(input_volume2[t2].sumsquares());  
       score=0;
       for(k=0;k<=input_volume1.maxz();k++)
         for(j=0;j<=input_volume1.maxy();j++)
           for(i=0;i<=input_volume1.maxx();i++)
	     score+=input_volume1(i,j,k,t1)*input_volume2(i,j,k,t2); 
       score=fabs(score)/(ss1*ss2); 
       if (score>thresh)
         cout << setw(3) << t1+1 << " " << setw(3) << t2+1 << " " <<  setiosflags (ios::fixed) << setprecision(2) << score << endl;
    }
  }

  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  string progname=argv[0];
  if (argc < 3 || argc > 4) 
  { 
    print_usage(progname);
    return 1; 
  }
   
  string iname=string(argv[1]);
  return call_fmrib_main(dtype(iname),argc,argv); 
}
