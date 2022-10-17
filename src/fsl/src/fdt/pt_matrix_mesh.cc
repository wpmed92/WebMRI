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

#include "pt_matrix_mesh.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace TRACTVOLS;
using namespace mesh;

void mesh_matrix2(){
  
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

  
  volume<int> beenhere;
  Mesh mseeds;

  if(opts.meshfile.value()!=""){
    mseeds.load(opts.meshfile.value());    
    mseeds.load_fs_label(opts.seedfile.value());
  }
  
  
  volume<int> RUBBISH;
  ColumnVector dim_rubbish;
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value()); //rubbishfile should be in MNI space.
    dim_rubbish.ReSize(3);
    dim_rubbish <<RUBBISH.xdim()<<RUBBISH.ydim()<<RUBBISH.zdim();
  }
  
  

  int numseeds=0;
  
  if(opts.meshfile.value()!=""){
  for (vector<Mpoint*>::iterator i = mseeds._points.begin(); i!=mseeds._points.end(); i++ )
      if((*i)->get_value() >0) numseeds++;
  }


  volume<int> ConMat;
  //  volume<int> CoordMat(numseeds,3,1); don't need Coordmat thing for mesh as info is in Seedfile
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
  
  
  
//   int myrow=0;
//   for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
//     for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
//       for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
// 	if(Seeds(Wx,Wy,Wz)>0){
// 	  CoordMat(myrow,0,0)=Wx;
// 	  CoordMat(myrow,1,0)=Wy;
// 	  CoordMat(myrow,2,0)=Wz;
// 	  myrow++;
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
  int other_conrow=0;
  int other_concol=0;

  ColumnVector mni_origin(3),fs_coord_mm(3),xyz_dti,xyz_rubbish,xyz_othermatrix_tracts(3),dim_othermatrix_tracts(3);
  dim_othermatrix_tracts <<lrmask.xdim()<<lrmask.ydim()<<lrmask.zdim();
  mni_origin << 92 << 128 << 37;
  
  for (vector<Mpoint*>::iterator i = mseeds._points.begin(); i!=mseeds._points.end(); i++ ){
    if((*i)->get_value() >0){
      
      
      fs_coord_mm<<(*i)->get_coord().X<<(*i)->get_coord().Y<<(*i)->get_coord().Z; 
      xyz_dti=mni_to_imgvox(fs_coord_mm,mni_origin, Seeds_to_DTI,vols.dimensions()); //xyz_dti in voxels, not mm
      xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3); //xyz_dti in voxels,not mm
      
      Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);

      for( int p = 0; p < nparticles ; p++ ){
	
	//Don't have a direction loop as always want to track in from cortex.
	
	x=xst;y=yst;z=zst;
	part.change_xyz(x,y,z);	    
	part.set_dir((*i)->local_normal().X,(*i)->local_normal().Y,(*i)->local_normal().Z);//Set the start dir so that we track inwards from cortex
	for( int it = 1 ; it <= nsteps; it++){
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

	    if(opts.rubbishfile.value()!=""){
	      xyz_rubbish=vox_to_vox(xyz_dti,vols.dimensions(),dim_rubbish,Seeds_to_DTI.i());
	      int x_s =(int)round((float)xyz_rubbish(1));
	      int y_s =(int)round((float)xyz_rubbish(2));
	      int z_s =(int)round((float)xyz_rubbish(3));
	      if(RUBBISH(x_s,y_s,z_s)>0) break;
	    }
	    
	    //find out where we are in the space of "lrmask" (which is in alignment with Seeds, but not 
	    // necessarily the same resolution
	    xyz_othermatrix_tracts=vox_to_vox(xyz_dti,vols.dimensions(),dim_othermatrix_tracts,Seeds_to_DTI.i()); 
	    int x_omt=(int)round((float)xyz_othermatrix_tracts(1));
	    int y_omt=(int)round((float)xyz_othermatrix_tracts(2));
	    int z_omt=(int)round((float)xyz_othermatrix_tracts(3));
	    path(it,1)=x_omt; 
	    path(it,2)=y_omt;
	    path(it,3)=z_omt;
	    
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
	      xyz_rubbish=vox_to_vox(xyz_dti,vols.dimensions(),dim_rubbish,Seeds_to_DTI.i());
	      int x_s =(int)round((float)xyz_rubbish(1));
	      int y_s =(int)round((float)xyz_rubbish(2));
	      int z_s =(int)round((float)xyz_rubbish(3));
	      
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
      
	indexset(beenhere,path,0);
	part.reset();
	
      } // Close Particle Number Loop
      other_conrow++;
      
    }
    
    
  } //Close Seed number Loop
  
  save_volume(ConMat,logger.appendDir(opts.outfile.value()));
  //  save_volume(CoordMat,logger.appendDir("coords_for_"+opts.outfile.value()));
  save_volume(CoordMat_tracts_om,logger.appendDir("tract_space_coords_for_"+opts.outfile.value()));
  save_volume4D(lookup,logger.appendDir("lookup_tractspace_"+opts.outfile.value()));
  
}





void mesh_lengths(){
  probtrackOptions& opts =probtrackOptions::getInstance();
 
  ////////////////////////////
  Log& logger = LogSingleton::getInstance();
  if(opts.verbose.value()>1){
    logger.makeDir("particles","particle0",true,false);
  }
  ////////////////////////////////////
  float xst,yst,zst,x,y,z;
  int nparticles=opts.nparticles.value();
  int nsteps=opts.nsteps.value();
  ///////////////////////////////////
  volume<int> mask;
  volume<int> RUBBISH;
  volume<int> skipmask;
  volume<int> prob,beenhere;
  read_volume(mask,opts.maskfile.value());
  if(opts.seedref.value()!=""){
    read_volume(prob,opts.seedref.value());
    prob=0;beenhere=prob;
    
  }
  else{
    prob=mask;prob=0;
    beenhere=prob;
  }
  
  Matrix Seeds_to_DTI;
  read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value()); // Here seeds_to_dti should take the standard volume to diff space 
  
  float lcrat=5;
  volume4D<float> loopcheck((int)ceil(mask.xsize()/lcrat)+1,(int)ceil(mask.ysize()/lcrat)+1,(int)ceil(mask.zsize()/lcrat)+1,3);
  loopcheck=0;
  
  if(opts.rubbishfile.value()!="") read_volume(RUBBISH,opts.rubbishfile.value());
  if(opts.skipmask.value()!="") read_volume(skipmask,opts.skipmask.value());
  
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
  
  Matrix path(nsteps,3);
  path=1;

  float tmp2;
  float randtmp1,randtmp2,randtmp3;
  ColumnVector th_ph_f;
  
  Mesh mseeds;
  int ftype=mseeds.load(opts.meshfile.value()); 
  mseeds.load_fs_label(opts.seedfile.value());
  ColumnVector mni_origin(3),fs_coord_mm(3),xyz_dti,xyz_seeds,dim_seeds(3);
  dim_seeds<<prob.xdim()<<prob.ydim()<<prob.zdim(); //In seedref space if exists. Else in dti space
  mni_origin << 92 << 128 << 37;
  
  
  
  for (vector<Mpoint*>::iterator i = mseeds._points.begin(); i!=mseeds._points.end(); i++ ){
    if((*i)->get_value() >0){
      
      fs_coord_mm<<(*i)->get_coord().X<<(*i)->get_coord().Y<<(*i)->get_coord().Z; 
      xyz_dti=mni_to_imgvox(fs_coord_mm,mni_origin, Seeds_to_DTI,vols.dimensions()); //xyz_dti in voxels, not mm
      xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3); //xyz_dti in voxels,not mm
      
      
      
      Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
      
      int length=0;
      for( int p = 0; p < nparticles ; p++ ){
	if(opts.verbose.value()>0){
	  cout<<"particle number "<<p<<endl;
	}
	
	if(opts.verbose.value()>1)
	  logger.setLogFile("particle"+num2str(p));
	
	//Don't have a direction loop as in other cases, as always want to track in from cortex.
	
	x=xst;y=yst;z=zst;
	part.change_xyz(x,y,z);	    
	part.set_dir((*i)->local_normal().X,(*i)->local_normal().Y,(*i)->local_normal().Z);//Set the start dir so that we track inwards from cortex
	for( int it = 1 ; it <= nsteps; it++){
	  if( (mask( round(part.x()), round(part.y()), round(part.z())) == 1) ){
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
	    
	    
	    if(opts.verbose.value()>1){
	      logger << part;
	    } 
	    
	    
	    int x_s,y_s,z_s;
	    
	    if(opts.seedref.value()!=""){
	      x=part.x();y=part.y();z=part.z();
	      xyz_dti <<x<<y<<z;
	      xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
	      x_s =(int)round((float)xyz_seeds(1));
	      y_s =(int)round((float)xyz_seeds(2));
	      z_s =(int)round((float)xyz_seeds(3));
	    }
	    else{
	      x_s=(int)round(part.x());
	      y_s=(int)round(part.y());
	      z_s=(int)round(part.z());
	    }
	      
	    if(opts.rubbishfile.value()!="")
	      {
		if(RUBBISH(x_s,y_s,z_s)>0) break;
	      }
	      
	    path(it,1)=x_s; 
	    path(it,2)=y_s;
	    path(it,3)=z_s;

	    if(beenhere(x_s,y_s,z_s)==0){
	      prob(x_s,y_s,z_s)+=1;
	      beenhere(x_s,y_s,z_s)=1;
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
	      if(!part.check_dir(th_ph_f(1),th_ph_f(2),opts.c_thr.value()) && it!=1){ 
		//Don't do curvature checking on the first step as we have set the old direction to the surface normal 
		break;
	      }
	      
	      if((th_ph_f(1)!=0&&th_ph_f(2)!=0)){
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
	    length++;
	      
	  }
	    
	} // Close Step Number Loop
	  
	if(opts.loopcheck.value()){
	  loopcheck=0;}
	  
	part.reset();
	indexset(beenhere,path,0);
	
      } // Close Particle Number Loop
      (*i)->set_value(length/nparticles);

      //      string thisout=opts.outfile.value()+num2str(xst)+(string)"_"+num2str(yst)+(string)"_"+num2str(zst);
      // save_volume(prob,thisout);
    }
  } //Close Seed number Loop
  mseeds.save_fs_label(logger.appendDir(opts.outfile.value()));
}

