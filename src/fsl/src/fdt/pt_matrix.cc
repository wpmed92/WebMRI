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

#include "pt_matrix.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace TRACTVOLS;
using namespace mesh;

void matrix2(){
  
  probtrackOptions& opts =probtrackOptions::getInstance();
  Log& logger = LogSingleton::getInstance();

  ////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  //  logger.makeDir(opts.logdir.value(),"logfile",true,false);
  ////////////////////////////////////
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
  //////////////////////////////////////////////
  // Segmented Volumes                        //
  //////////////////////////////////////////////
 //   vector<string> masknames;
//    read_masks(masknames,opts.cortexfile.value());
//    vector<volume<int> > cortex_masks;
//    vector<volume<int> > thal_segs;
//    volume<int> tmpcort;
//    cerr<<"Number of masks "<<masknames.size()<<endl;
//    for( unsigned int m = 0; m < masknames.size(); m++ ){
//      cerr<<"Reading "<<masknames[m]<<endl;
//      read_volume(tmpcort,masknames[m]);
//      cortex_masks.push_back(tmpcort);
//      tmpcort=0;
//      thal_segs.push_back(tmpcort);
//  }


  volume<int> Seeds;
  read_volume(Seeds,opts.seedfile.value());

  volume<int> RUBBISH;
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value());
  }

  volume<int> prob;
  volume<int> beenhere;
  prob=Seeds;prob=0;

  int numseeds=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds.value(Wx,Wy,Wz)>0){
	  numseeds++;
	}
      }
    }
  }

  volume<int> ConMat;
  volume<int> CoordMat(numseeds,3,1);
  volume<int> CoordMat_tracts_om; //for storing tractspace coords for othermatrix
  volume<int> lrmask;

  
  Matrix tempy;
  volume4D<int> lookup;
  read_volume(lrmask,opts.lrmask.value());
  beenhere=lrmask;beenhere=0;
  int numnz=0;
  for(int Wz=lrmask.minz();Wz<=lrmask.maxz();Wz++){
    for(int Wy=lrmask.miny();Wy<=lrmask.maxy();Wy++){
      for(int Wx=lrmask.minx();Wx<=lrmask.maxx();Wx++){
	if(lrmask.value(Wx,Wy,Wz)>0){
	  numnz++;
	}
      }
    }
  }
  if(numnz> pow(2,(float)sizeof(short)*8-1)-1){
    cerr<<"Output matrix too big for AVW - stopping."<<endl;
    cerr<<" Remember - you can store your tracts in "<<endl;
    cerr<<" low res even if you want your seeds in high res"<<endl;
    cerr<<" Just subsample the structural space mask"<<endl;
    cerr<<" Although, it must stay in line with the seeds"<<endl;
    exit(-1);
  }
  
  //    cerr<<"WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  //      cerr<<"You will need significantly more than "<<numseeds*numnz*2<<" bytes of free memory!!"<<endl;
  
  ConMat.reinitialize(numseeds,numnz,1);
  ConMat=0;
  tempy.ReSize(numnz,1);
  for(int i=1;i<=numnz;i++){tempy(i,1)=i-1;}
  lookup.addvolume(lrmask);
  lookup.setmatrix(tempy.t(),lrmask);
  
  CoordMat_tracts_om.reinitialize(numnz,3,1);//store the tract space coordmat.
  int mytrow=0;
  for(int Wz=lrmask.minz();Wz<=lrmask.maxz();Wz++){
    for(int Wy=lrmask.miny();Wy<=lrmask.maxy();Wy++){
      for(int Wx=lrmask.minx();Wx<=lrmask.maxx();Wx++){
	if(lrmask(Wx,Wy,Wz)>0){
	  CoordMat_tracts_om(mytrow,0,0)=Wx;
	  CoordMat_tracts_om(mytrow,1,0)=Wy;
	  CoordMat_tracts_om(mytrow,2,0)=Wz;
	  mytrow++;
	}
      }
    }
  }
  
  
  
  int myrow=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds(Wx,Wy,Wz)>0){
	  CoordMat(myrow,0,0)=Wx;
	  CoordMat(myrow,1,0)=Wy;
	  CoordMat(myrow,2,0)=Wz;
	  myrow++;
	}
      }
    }
  }
  
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
  int other_conrow=0;
  int other_concol=0;

  for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
    cout<<Sz<<endl;
    for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
      for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	if(Seeds(Sx,Sy,Sz)>0){
	  
	  ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti,xyz_othermatrix_tracts(3),dim_othermatrix_tracts(3);
	  xyz_seeds << Sx << Sy << Sz;
	  dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();
	  dim_othermatrix_tracts <<lrmask.xdim()<<lrmask.ydim()<<lrmask.zdim();
	  
	  
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
		
		  ////////////////////////////////////////

		  x=part.x();y=part.y();z=part.z();
		  xyz_dti <<x<<y<<z;
		  xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		  int x_s =(int)round((float)xyz_seeds(1));
		  int y_s =(int)round((float)xyz_seeds(2));
		  int z_s =(int)round((float)xyz_seeds(3));

		  if(opts.rubbishfile.value()!=""){
		    if(RUBBISH(x_s,y_s,z_s)>0) break;
		  }
		  
		  //find out where we are in the space of "lrmask" (which is in alignment with Seeds, but not 
		  // necessarily the same resolution
		  xyz_othermatrix_tracts=vox_to_vox(xyz_dti,vols.dimensions(),dim_othermatrix_tracts,Seeds_to_DTI.i()); 
		  int x_omt=(int)round((float)xyz_othermatrix_tracts(1));
		  int y_omt=(int)round((float)xyz_othermatrix_tracts(2));
		  int z_omt=(int)round((float)xyz_othermatrix_tracts(3));
		  path(it+(direc-1)*nsteps/2,1)=x_omt; 
		  path(it+(direc-1)*nsteps/2,2)=y_omt;
		  path(it+(direc-1)*nsteps/2,3)=z_omt;
		  
		  // Find out where this lrmask voxel is in the unwrapped matrix
		  other_concol=lookup(x_omt,y_omt,z_omt,0);
		  //load up the matrix
		  if(other_concol!=0){
		    if(beenhere(x_omt,y_omt,z_omt)==0){
		      ConMat(other_conrow,other_concol,0)+=1;
		      beenhere(x_omt,y_omt,z_omt)=1;
		    }
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
	    }//close direc loop
	    indexset(beenhere,path,0);
	    part.reset();

	  } // Close Particle Number Loop
	  other_conrow++;
	}
      }
    }
    


  } //Close Seed number Loop
  
  save_volume(ConMat,logger.appendDir(opts.outfile.value()));
  save_volume(CoordMat,logger.appendDir("coords_for_"+opts.outfile.value()));
  save_volume(CoordMat_tracts_om,logger.appendDir("tract_space_coords_for_"+opts.outfile.value()));
  save_volume4D(lookup,logger.appendDir("lookup_tractspace_"+opts.outfile.value()));
  
}




void matrix1(){
  
  probtrackOptions& opts =probtrackOptions::getInstance();
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


  volume<int> Seeds;
  read_volume(Seeds,opts.seedfile.value());

  volume<int> RUBBISH;
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value());
  }

  
  volume<int> prob;
  volume<int> beenhere;
  prob=Seeds;prob=0;
  beenhere=Seeds;beenhere=0;
 
//  prob=tmpcort;beenhere=tmpcort;
  int numseeds=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds.value(Wx,Wy,Wz)>0){
	  numseeds++;
	}
      }
    }
  }

  volume<int> ConMat;
  volume<int> CoordMat(numseeds,3,1);
  ConMat.reinitialize(numseeds,numseeds,1);
  ConMat=0;
  
  
  int myrow=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds(Wx,Wy,Wz)>0){
	  CoordMat(myrow,0,0)=Wx;
	  CoordMat(myrow,1,0)=Wy;
	  CoordMat(myrow,2,0)=Wz;
	  myrow++;
	}
      }
    }
  }
  
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
  int Conrow=0;
  int hitcount=0;

  for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
    cout<<Sz<<endl;
    for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
      for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	if(Seeds(Sx,Sy,Sz)>0){
	  //	  other_conrow++;
	  //	  int repeatcount=0;
	  ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	  xyz_seeds << Sx << Sy << Sz;
	  dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();	  
	  xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	  xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	  Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
	  hitcount=0;
	  for( int p = 0; p < nparticles ; p++ ){
	    
	    for(int direc=1;direc<=2;direc++){
	      x=xst;y=yst;z=zst;
	      part.change_xyz(x,y,z);	    
	      if(direc==2){
		part.restart_reverse();  //go again in the opposite direction
	      }
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
		
		  ////////////////////////////////////////


		  //		if(opts.verbose.value()>2){
		  //		  logger << part;
		  //		} 

		 
		  x=part.x();y=part.y();z=part.z();
		  xyz_dti <<x<<y<<z;
		  xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		  int x_s =(int)round((float)xyz_seeds(1));
		  int y_s =(int)round((float)xyz_seeds(2));
		  int z_s =(int)round((float)xyz_seeds(3));
		  
		  if(opts.rubbishfile.value()!=""){
		    if(RUBBISH(x_s,y_s,z_s)>0) break;
		  }
		  
		  if(beenhere(x_s,y_s,z_s)==0){
		    prob(x_s,y_s,z_s)+=1;
		    beenhere(x_s,y_s,z_s)=1;
		  }
		  
		  if(opts.skipmask.value() == ""){
		    th_ph_f=vols.sample(part.x(),part.y(),part.z());
		  }
		  else{
		    if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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
	    }
	    indexset(beenhere,path,0);
	    part.reset();
	   
	  } // Close Particle Number Loop
	
	  
	  ///Store connectivity information in massive matrix -  matrix mode.  
	  int Concol=0;
	  for(int Wz=prob.minz();Wz<=prob.maxz();Wz++){
	    for(int Wy=prob.miny();Wy<=prob.maxy();Wy++){
	      for(int Wx=prob.minx();Wx<=prob.maxx();Wx++){
		if(Seeds(Wx,Wy,Wz)>0){
		  if(prob(Wx,Wy,Wz)>0){
		    ConMat(Conrow,Concol,0)=prob(Wx,Wy,Wz);
		  }
		  Concol++;
		}
		prob(Wx,Wy,Wz)=0;
		
	      }
	    }
	  }
	  
	  Conrow++;	
	}
	
      }
    }
  } //Close Seed number Loop
  save_volume(ConMat,logger.appendDir(opts.outfile.value()));
  save_volume(CoordMat,logger.appendDir("coords_for_"+opts.outfile.value()));
}




void maskmatrix(){
  
  probtrackOptions& opts =probtrackOptions::getInstance();
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

  volume<int> Seeds;
  read_volume(Seeds,opts.seedfile.value());

  volume<int> RUBBISH;
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value());
  }

  
  volume<int> prob;
  volume<int> beenhere;
  prob=Seeds;prob=0;
  beenhere=prob;
 

  int numseeds=0;
  int  maxclusternum=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds.value(Wx,Wy,Wz)>0){
	  numseeds++;
	  if(Seeds.value(Wx,Wy,Wz) > maxclusternum){
	    maxclusternum=Seeds.value(Wx,Wy,Wz);
	  }
	}
      }
    }
  }

  Matrix meanconmat(maxclusternum,maxclusternum),maxconmat(maxclusternum,maxclusternum);
  maxconmat=0;meanconmat=0;
    
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
  int hitcount=0;

  
  
  for(int seedclust=1;seedclust<=maxclusternum;seedclust++){
    int clustsize=0;
    ColumnVector clustsums(maxclusternum),clustmaxes(maxclusternum);
    clustsums=0;clustmaxes=0;
    for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
      cout<<seedclust<<" "<<Sz<<endl;
      for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
	for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	  if(Seeds(Sx,Sy,Sz)==seedclust){

	    clustsize++;
	    ColumnVector flags(maxclusternum); 
	    ColumnVector clustcounts(maxclusternum);  clustcounts=0;


	    ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	    xyz_seeds << Sx << Sy << Sz;
	    dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();	  
	    xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	    Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
	    hitcount=0;
	    for( int p = 0; p < nparticles ; p++ ){
	     flags=0;
	      for(int direc=1;direc<=2;direc++){
		x=xst;y=yst;z=zst;
		part.change_xyz(x,y,z);	    
		if(direc==2){
		  part.restart_reverse();  //go again in the opposite direction
		}
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
		    if(opts.rubbishfile.value()!=""){
		      if(RUBBISH(x_s,y_s,z_s)>0) break;
		    }
		    
		    int cn=Seeds(x_s,y_s,z_s);
		    if(cn!=0){
		      if(flags(cn)==0){
			clustcounts(cn)=clustcounts(cn)+1;flags(cn)=1;
		      }
		    }
		    
		    if(beenhere(x_s,y_s,z_s)==0){
		      prob(x_s,y_s,z_s)+=1;
		      beenhere(x_s,y_s,z_s)=1;
		    }
		    
		    if(opts.skipmask.value() == ""){
		      th_ph_f=vols.sample(part.x(),part.y(),part.z());
		    }
		    else{
		      if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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
	      }
	      indexset(beenhere,path,0);
	      part.reset();
	      
	    } // Close Particle Number Loop
	
	  
	    clustsums+=clustcounts;
	    for(int cn=1;cn<=maxclusternum;cn++){
	      if(clustcounts(cn)>clustmaxes(cn)){
		clustmaxes(cn)=clustcounts(cn);
	      }
	    }
	    
	  } 
	  
	}//close x seed
      } //close y seed
    }  //close z seed
    for(int targclust=1;targclust<=maxclusternum;targclust++){
      meanconmat(seedclust,targclust)=clustsums(targclust)/clustsize;
      maxconmat(seedclust,targclust)=clustmaxes(targclust);
    }
    
  } //close cluster number
  write_ascii_matrix(meanconmat,logger.appendDir("mean_"+opts.outfile.value()));
  write_ascii_matrix(maxconmat,logger.appendDir("max_"+opts.outfile.value()));
}


