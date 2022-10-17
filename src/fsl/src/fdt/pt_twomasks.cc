/*  Copyright (C) 2004 University of Oxford  */

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

#include "pt_twomasks.h"
#include "pt_seeds_to_targets.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace TRACTVOLS;
using namespace mesh;

void twomasks_symm(){
  
  probtrackOptions& opts = probtrackOptions::getInstance();
  Log& logger = LogSingleton::getInstance();
  
  float xst,yst,zst,x,y,z;
  int nparticles=opts.nparticles.value();
  int nsteps=opts.nsteps.value();
  ///////////////////////////////////
  volume<int> mask;
  read_volume(mask,opts.maskfile.value());  
  volume<int> skipmask;
  if(opts.skipmask.value()!="") read_volume(skipmask,opts.skipmask.value());
  float lcrat=5;
  volume4D<float> loopcheck;
  if(opts.loopcheck.value()){
    loopcheck.reinitialize(int(ceil(mask.xsize()/lcrat)+1),int(ceil(mask.ysize()/lcrat)+1),int(ceil(mask.zsize()/lcrat)+1),3);
    loopcheck=0;
  }
  
  volume<int> Seeds,mask2;
  read_volume(Seeds,opts.seedfile.value());

  if(opts.mask2.value()!=""){
    read_volume(mask2,opts.mask2.value());
    Seeds=Seeds+(2*mask2);
  }

  volume<int> RUBBISH; //just stops the streamline but counts where it has already been
  volume<int> UBERRUBBISH; //kills streamline and excludes where it has already been
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value());
  }
  if(opts.uberrubbishfile.value()!=""){
    read_volume(UBERRUBBISH,opts.uberrubbishfile.value());
  }
  
 
  volume<int> prob;
  volume<int> beenhere;
  beenhere=Seeds;
  beenhere=0;
  prob=beenhere;
//   int numseeds=0;
//   for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
//     for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
//       for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
// 	if(Seeds.value(Wx,Wy,Wz)>0){
// 	  numseeds++;
// 	}
//       }
//     }
//   }


  cout<<"Loading MCMC volumes"<<endl;
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
 
  Matrix Seeds_to_DTI;
  
  if(opts.seeds_to_dti.value()!=""){
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
  }
  else{
    Seeds_to_DTI=Identity(4);
  }
  
  

  Matrix path(nsteps,3);
  path=1;
  
  float tmp2;
  ColumnVector th_ph_f;

    for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
      cout<<Sz<<endl;
      for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
	for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	  if(Seeds(Sx,Sy,Sz)>0){
	    int	 maskno=Seeds(Sx,Sy,Sz);
	    ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	    xyz_seeds << Sx << Sy << Sz;
	    dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();
	    
	    xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	    Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
	    
	    for( int p = 0; p < nparticles ; p++ ){
	      for(int direc=1;direc<=2;direc++){
		x=xst;y=yst;z=zst;
		part.change_xyz(x,y,z);
		if(direc==2){
		  part.restart_reverse();  //go again in the opposite direction
		}
		int partlength=0;
		bool keepflag=false;  /// only keep it if this direction went through the masks
		for( int it = 1 ; it <= nsteps/2; it++){
		  if( (mask( round(part.x()), round(part.y()), round(part.z())) > 0) ){
		    
		    ///////////////////////////////////
		    //loopchecking
		    ///////////////////////////////////
		    if(opts.loopcheck.value()){
		      float oldrx=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),0);
		      float oldry=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),1);
		      float oldrz=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),2);
		      if(part.rx()*oldrx+part.ry()*oldry+part.rz()*oldrz<0)
			{
			  // p--;
			  break;
			}
		  
		      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),0)=part.rx();
		      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),1)=part.ry();
		      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),2)=part.rz();  
		  
		    }
		    
		    
		    x=part.x();y=part.y();z=part.z();
		    xyz_dti <<x<<y<<z;
		    xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		    int x_s =(int)round((float)xyz_seeds(1));
		    int y_s =(int)round((float)xyz_seeds(2));
		    int z_s =(int)round((float)xyz_seeds(3));
		    
		    if(opts.uberrubbishfile.value()!=""){
		      if(UBERRUBBISH(x_s,y_s,z_s)>0) {
			keepflag=false;
			break;
		      }
		    }
		    if(opts.rubbishfile.value()!=""){
		      if(RUBBISH(x_s,y_s,z_s)>0) break;
		    }
		    path(it,1)=x_s; 
		    path(it,2)=y_s;
		    path(it,3)=z_s;
		    partlength+=1;
		    if( Seeds(x_s,y_s,z_s)==(maskno*-1 + 3) ){ //mult by -1 and add 3 makes 1 into 2 and 2 into 1
		      keepflag=true;
		    }
		    if(opts.skipmask.value() == ""){
		      th_ph_f=vols.sample(part.x(),part.y(),part.z());
		    }
		    else{
		      if(skipmask(x_s,y_s,z_s)==0)
			th_ph_f=vols.sample(part.x(),part.y(),part.z());
		    }
		    
		    
		    tmp2=rand(); tmp2/=RAND_MAX;
		    if(th_ph_f(3)>tmp2){
		      if(!part.check_dir(th_ph_f(1),th_ph_f(2),opts.c_thr.value())){
			break;
		      }

		      if((th_ph_f(1)!=0&&th_ph_f(2)!=0)){
			if( (mask( round(part.x()), round(part.y()), round(part.z())) != 0) ){
			  if(!opts.modeuler.value())
			    part.jump(th_ph_f(1),th_ph_f(2));
			  else
			    {
			      ColumnVector test_th_ph_f;
			      part.testjump(th_ph_f(1),th_ph_f(2));
			      test_th_ph_f=vols.sample(part.testx(),part.testy(),part.testz());
			      part.jump(test_th_ph_f(1),test_th_ph_f(2));
			    }
			}
		    
			

		      
		    
		    
		      }
		    }
		
		
		  }
	      
		} // Close Step Number Loop
		if(opts.loopcheck.value()){
		  loopcheck=0;
		}
	  
	      
	      // Do the counting after a particle has finished
	      // And only if keepflag is set
		if(keepflag){
       		  for(int s=1;s<=partlength;s++){
		    int x_s=int(path(s,1));
		    int y_s=int(path(s,2));
		    int z_s=int(path(s,3));
		    if(beenhere(x_s,y_s,z_s)==0){
		      prob(x_s,y_s,z_s)+=1;
		      beenhere(x_s,y_s,z_s)=1;
		    }
		  } 
		  
		  indexset(beenhere,path,0);
		}
	      } //close direction loop

	      part.reset(); 

	   	    
	    
	    } // Close Particle Number Loop


	  }

	}// close x loop
      }//close y loop
   

    }//close z loop
  save_volume(prob,logger.appendDir(opts.outfile.value()));
}


void waypoints(){
  probtrackOptions& opts = probtrackOptions::getInstance();
  Log& logger = LogSingleton::getInstance();
  
  float xst,yst,zst,x,y,z;
  int nparticles=opts.nparticles.value();
  int nsteps=opts.nsteps.value();
  ///////////////////////////////////
  volume<int> mask;
  read_volume(mask,opts.maskfile.value());  
  volume<int> skipmask;
  if(opts.skipmask.value()!="") read_volume(skipmask,opts.skipmask.value());
  float lcrat=5;
  volume4D<float> loopcheck;
  if(opts.loopcheck.value()){
    loopcheck.reinitialize(int(ceil(mask.xsize()/lcrat)+1),int(ceil(mask.ysize()/lcrat)+1),int(ceil(mask.zsize()/lcrat)+1),3);
    loopcheck=0;
  }
  
  volume<int> Seeds,mask2;
  read_volume(Seeds,opts.seedfile.value());
  
  
  vector<string> masknames;
  if(fsl_imageexists(opts.mask2.value())){
    masknames.push_back( opts.mask2.value() );
  }
  else{
    read_masks(masknames,opts.mask2.value());
  }
  vector<volume<int> > waymasks;
  volume<int> tmpway;
  vector<bool> passed_flags; //store whether the current sample has passed through each of the waypoints
  for( unsigned int m = 0; m < masknames.size(); m++ ){
    if(opts.verbose.value()>0)
      cout<<masknames[m]<<endl;
    read_volume(tmpway,masknames[m]);
    waymasks.push_back(tmpway);
    passed_flags.push_back(false);
  }
  

  volume<int> RUBBISH; //just stops the streamline but counts where it has already been
  volume<int> UBERRUBBISH; //kills streamline and excludes where it has already been
  bool uberkeepflag=true;
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value());
  }
  if(opts.uberrubbishfile.value()!=""){
    read_volume(UBERRUBBISH,opts.uberrubbishfile.value());
  }
  
  
  volume<int> prob;
  volume<int> beenhere;
  beenhere=Seeds;
  beenhere=0;
  prob=beenhere;
//   int numseeds=0;
//   for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
//     for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
//       for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
// 	if(Seeds.value(Wx,Wy,Wz)>0){
// 	  numseeds++;
// 	}
//       }
//     }
//   }


  cout<<"Loading MCMC volumes"<<endl;
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
 
  Matrix Seeds_to_DTI;
  
  if(opts.seeds_to_dti.value()!=""){
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
  }
  else{
    Seeds_to_DTI=Identity(4);
  }
  
  

  Matrix path(nsteps,3);
  path=1;
  
  float tmp2;
  ColumnVector th_ph_f;
  int keeptotal=0;
    for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
      cout<<Sz<<endl;
      for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
	for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	  if(Seeds(Sx,Sy,Sz)>0){
	    ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	    xyz_seeds << Sx << Sy << Sz;
	    dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();
	  
	    xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	    Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
	    
	    for( int p = 0; p < nparticles ; p++ ){

	      for(int direc=1;direc<=2;direc++){
		x=xst;y=yst;z=zst;
		part.change_xyz(x,y,z);
		if(direc==2){
		  part.restart_reverse();  //go again in the opposite direction
		}
		int partlength=0;
		for(unsigned int pf=0;pf<passed_flags.size();pf++) {
		  passed_flags[pf]=false;  /// only keep it if this direction went through all the masks
		}
		uberkeepflag=true;
		for( int it = 1 ; it <= nsteps/2; it++){
		  if( (mask( round(part.x()), round(part.y()), round(part.z())) > 0) ){
		    
		    ///////////////////////////////////
		    //loopchecking
		    ///////////////////////////////////
		    if(opts.loopcheck.value()){
		      float oldrx=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),0);
		      float oldry=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),1);
		      float oldrz=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),2);
		      if(part.rx()*oldrx+part.ry()*oldry+part.rz()*oldrz<0)
			{
			  // p--;
			  break;
			}
		  
		      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),0)=part.rx();
		      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),1)=part.ry();
		      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),2)=part.rz();  
		  
		    }
		    
		    
		    x=part.x();y=part.y();z=part.z();
		    xyz_dti <<x<<y<<z;
		    xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		    int x_s =(int)round((float)xyz_seeds(1));
		    int y_s =(int)round((float)xyz_seeds(2));
		    int z_s =(int)round((float)xyz_seeds(3));
		    if(opts.uberrubbishfile.value()!=""){
		      if(UBERRUBBISH(x_s,y_s,z_s)>0) {
			uberkeepflag=false;
			break;
		      }
		    }
		    if(opts.rubbishfile.value()!=""){
		      if(RUBBISH(x_s,y_s,z_s)>0) break;
		    }
		    path(it,1)=x_s; 
		    path(it,2)=y_s;
		    path(it,3)=z_s;
		    partlength+=1;

		    //update every passed_flag
		    for( unsigned int wm=0;wm<waymasks.size();wm++ ){
		      if( waymasks[wm](x_s,y_s,z_s)>0 ) 
			passed_flags[wm]=true;
		    }


		    if(opts.skipmask.value() == ""){
		      th_ph_f=vols.sample(part.x(),part.y(),part.z());
		    }
		    else{
		      if(skipmask(x_s,y_s,z_s)==0)
			th_ph_f=vols.sample(part.x(),part.y(),part.z());
		    }
		    
		    
		    tmp2=rand(); tmp2/=RAND_MAX;
		    if(th_ph_f(3)>tmp2){
		      if(!part.check_dir(th_ph_f(1),th_ph_f(2),opts.c_thr.value())){
			break;
		      }

		      if((th_ph_f(1)!=0&&th_ph_f(2)!=0)){
			if( (mask( round(part.x()), round(part.y()), round(part.z())) != 0) ){
			  if(!opts.modeuler.value())
			    part.jump(th_ph_f(1),th_ph_f(2));
			  else
			    {
			      ColumnVector test_th_ph_f;
			      part.testjump(th_ph_f(1),th_ph_f(2));
			      test_th_ph_f=vols.sample(part.testx(),part.testy(),part.testz());
			      part.jump(test_th_ph_f(1),test_th_ph_f(2));
			    }
			}
		    
			

		      
		    
		    
		      }
		    }
		
		
		  }
	      
		} // Close Step Number Loop
		if(opts.loopcheck.value()){
		  loopcheck=0;
		}
	  
	      
	      // Do the counting after a particle has finished
	      // And only if all passed_flags are set
		bool keepflag=uberkeepflag;
		for(unsigned int pf=0;pf<passed_flags.size();pf++){
		  keepflag=(keepflag && passed_flags[pf]);
		}
		
		if(keepflag){
		  keeptotal++;
       		  for(int s=1;s<=partlength;s++){
		    int x_s=int(path(s,1));
		    int y_s=int(path(s,2));
		    int z_s=int(path(s,3));
		    if(beenhere(x_s,y_s,z_s)==0){
		      prob(x_s,y_s,z_s)+=1;
		      beenhere(x_s,y_s,z_s)=1;
		    }
		  } 
		  indexset(beenhere,path,0);
		}
	      } //close direction loop

	      part.reset(); 

	   	    
	    
	    } // Close Particle Number Loop


	  }

	}// close x loop
      }//close y loop
   

    }//close z loop
    ColumnVector keeptotvec(1);
    keeptotvec(1)=keeptotal;
    
    save_volume(prob,logger.appendDir(opts.outfile.value()));
    write_ascii_matrix(keeptotvec,logger.appendDir("waytotal"));
}
