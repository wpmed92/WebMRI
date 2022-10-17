#include <iostream>
#include <fstream>
#include <iomanip>
#include "newimage/newimageall.h"
#include "utils/log.h"



using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
//using namespace NEWMAT;
//////////////////////////
/////////////////////////

void read_masks(vector<string>& masks, const string& filename){
  ifstream fs(filename.c_str());
  string tmp;
  if(fs){
    fs>>tmp;
    while(!fs.eof()){
      masks.push_back(tmp);
      fs>>tmp;
    }
  }
  else{
    cerr<<filename<<" does not exist"<<endl;
    exit(0);
  }
}


ColumnVector vox_to_vox(const ColumnVector& xyz1,const ColumnVector& dims1,const ColumnVector& dims2,const Matrix& xfm){
  ColumnVector xyz1_mm(4),xyz2_mm,xyz2(3);
  xyz1_mm<<xyz1(1)*dims1(1)<<xyz1(2)*dims1(2)<<xyz1(3)*dims1(3)<<1;
  xyz2_mm=xfm*xyz1_mm;
  xyz2_mm=xyz2_mm/xyz2_mm(4);
  xyz2<<xyz2_mm(1)/dims2(1)<<xyz2_mm(2)/dims2(2)<<xyz2_mm(3)/dims2(3);
  return xyz2;
}


int main ( int argc, char **argv ){
  if(argc<6){
    cerr<<"Usage: indexer <mask_for_going_ahead> <file_of_volume_names> <x> <y> <z>"<<endl;
    cerr<<endl;
    cerr<< "<x> <y> <z> in mni mm coords"<<endl;
    return 0;
  }
 
  vector<string> masknames;
  volume<int> inside_mask;
  read_volume(inside_mask,argv[1]);
  read_masks(masknames,argv[2]);
  vector<volume<float> > cortex_masks;
  volume<float> tmpcort;
  for( unsigned int m = 0; m < masknames.size(); m++ ){
    read_volume(tmpcort,masknames[m]);
    cortex_masks.push_back(tmpcort);
  }
  
  float x_mm=atof(argv[3]);
  float y_mm=atof(argv[4]);
  float z_mm=atof(argv[5]);
  float x_orig_mm=90;
  float y_orig_mm=124;
  float z_orig_mm=72;
  int x_roi_start_vox=45;
  int y_roi_start_vox=70;
  int z_roi_start_vox=50;
  int numsubjects=11;
  float xvox=(x_mm+x_orig_mm)/cortex_masks[0].xdim()-x_roi_start_vox;
  float yvox=(y_mm+y_orig_mm)/cortex_masks[0].ydim()-y_roi_start_vox;
  float zvox=(z_mm+z_orig_mm)/cortex_masks[0].zdim()-z_roi_start_vox;
  
  if(inside_mask((int)round(xvox),(int)round(yvox),(int)round(zvox))==0){
    cout<<"Sorry - Your input is not in our defined thalamus"<<endl;
    for(unsigned int i=0;i<masknames.size();i++){
      cout<<0<<endl;
    }
  }  
  else{
    cout<<"Input inside thalamus"<<endl;
    for(unsigned int i=0;i<masknames.size();i++){
      if(cortex_masks[i].interpolate(xvox,yvox,zvox)==0){
	cout<<std::setprecision(0)<<std::fixed;
      }
      else{
	cout<<std::setprecision(2)<<std::fixed;
      }
      cout<<cortex_masks[i].interpolate(xvox,yvox,zvox)/numsubjects<<endl;
    }
  }
  return 0;
}


 














