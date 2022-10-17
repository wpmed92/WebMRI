#include "streamlines.h"



namespace TRACT{


 
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
  
  
  Streamliner::Streamliner():opts(probtrackxOptions::getInstance()),logger(LogSingleton::getInstance()),
			     vols(opts.usef.value()){
    
    read_volume(m_mask,opts.maskfile.value());
    m_part.initialise(0,0,0,0,0,0,opts.steplength.value(),m_mask.xdim(),m_mask.ydim(),m_mask.zdim(),false);
    if(opts.skipmask.value()!="") read_volume(m_skipmask,opts.skipmask.value());
    m_lcrat=5;
    if(opts.loopcheck.value()){
      m_loopcheck.reinitialize(int(ceil(m_mask.xsize()/m_lcrat)+1),int(ceil(m_mask.ysize()/m_lcrat)+1),int(ceil(m_mask.zsize()/m_lcrat)+1),3);
      m_loopcheck=0;
    }
    if(opts.rubbishfile.value()!=""){
      read_volume(m_rubbish,opts.rubbishfile.value());
    }
      
    vector<string> masknames;
    if(opts.waypoints.value()!=""){
      if(fsl_imageexists(opts.waypoints.value())){
	masknames.push_back( opts.waypoints.value() );
      }
      else{
	read_masks(masknames,opts.waypoints.value());
      }

      for( unsigned int m = 0; m < masknames.size(); m++ ){
	volume<int>* tmpptr =new volume<int>;
	if(opts.verbose.value()>0)
	  cout<<masknames[m]<<endl;
	read_volume(*tmpptr,masknames[m]);
	m_waymasks.push_back(tmpptr);
	m_passed_flags.push_back(false);
	m_own_waymasks.push_back(true);
      }
    } 
    if(opts.seeds_to_dti.value()!=""){
      read_ascii_matrix(m_Seeds_to_DTI,opts.seeds_to_dti.value());
    }
    else{
      m_Seeds_to_DTI=Identity(4);
    }
    vols.initialise(opts.basename.value());
    m_path.reserve(opts.nparticles.value());
    m_x_s_init=0;
    m_y_s_init=0;
    m_z_s_init=0;
  }
  
  
  bool Streamliner::streamline(const float& x_init,const float& y_init,const float& z_init, const ColumnVector& dim_seeds,const int& fibst){ 
    
    //fibst tells tractvolsx which fibre to start with if there are more than one..
    //x_init etc. are in seed space...
    vols.reset(fibst);
    m_x_s_init=x_init; //seed x position in voxels
    m_y_s_init=y_init; // and y
    m_z_s_init=z_init; // and z
    ColumnVector xyz_seeds(3);
    xyz_seeds<<x_init<<y_init<<z_init;
    ColumnVector xyz_dti;
    ColumnVector th_ph_f;
    float xst,yst,zst,x,y,z,tmp2;
    xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),m_Seeds_to_DTI);    
    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
    m_path.clear();
    x=xst;y=yst;z=zst;
    m_part.change_xyz(x,y,z);
    int partlength=0;
    //NB - this only goes in one direction!!
    for(unsigned int pf=0;pf<m_passed_flags.size();pf++) {
      m_passed_flags[pf]=false;  /// only keep it if this streamline went through all the masks
    }
    
    for( int it = 1 ; it <= opts.nsteps.value()/2; it++){
      if( (m_mask( round(m_part.x()), round(m_part.y()), round(m_part.z())) > 0) ){
	///////////////////////////////////
	//loopchecking
	///////////////////////////////////
	if(opts.loopcheck.value()){
	  float oldrx=m_loopcheck((int)round(m_part.x()/m_lcrat),(int)round(m_part.y()/m_lcrat),(int)round(m_part.z()/m_lcrat),0);
	  float oldry=m_loopcheck((int)round(m_part.x()/m_lcrat),(int)round(m_part.y()/m_lcrat),(int)round(m_part.z()/m_lcrat),1);
	  float oldrz=m_loopcheck((int)round(m_part.x()/m_lcrat),(int)round(m_part.y()/m_lcrat),(int)round(m_part.z()/m_lcrat),2);
	  if(m_part.rx()*oldrx+m_part.ry()*oldry+m_part.rz()*oldrz<0)
	    {
	      break;
	    }
	    
	  m_loopcheck((int)round(m_part.x()/m_lcrat),(int)round(m_part.y()/m_lcrat),(int)round(m_part.z()/m_lcrat),0)=m_part.rx();
	  m_loopcheck((int)round(m_part.x()/m_lcrat),(int)round(m_part.y()/m_lcrat),(int)round(m_part.z()/m_lcrat),1)=m_part.ry();
	  m_loopcheck((int)round(m_part.x()/m_lcrat),(int)round(m_part.y()/m_lcrat),(int)round(m_part.z()/m_lcrat),2)=m_part.rz();  
	    
	}
	
	if(opts.verbose.value()>1)
	  logger<<m_part;
	
	x=m_part.x();y=m_part.y();z=m_part.z();
	xyz_dti <<x<<y<<z;
	xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,m_Seeds_to_DTI.i());
	int x_s =(int)round((float)xyz_seeds(1));
	int y_s =(int)round((float)xyz_seeds(2));
	int z_s =(int)round((float)xyz_seeds(3));
	  
	//update every passed_flag
	for( unsigned int wm=0;wm<m_waymasks.size();wm++ ){
	  if( (*m_waymasks[wm])(x_s,y_s,z_s)>0 ) {
	    m_passed_flags[wm]=true;
	  }
	}
	m_path.push_back(xyz_seeds);
	//	  m_path(it,1)=x_s; 
	//	  m_path(it,2)=y_s;
	//	  m_path(it,3)=z_s;
	partlength++;
	

	if(opts.rubbishfile.value()!=""){
	  if(m_rubbish(x_s,y_s,z_s)>0) break;
	}
	  
	
	  
	if(opts.skipmask.value() == ""){
	  th_ph_f=vols.sample(m_part.x(),m_part.y(),m_part.z(),m_part.rx(),m_part.ry(),m_part.rz());
	}
	else{
	  if(m_skipmask(x_s,y_s,z_s)==0)
	    th_ph_f=vols.sample(m_part.x(),m_part.y(),m_part.z(),m_part.rx(),m_part.ry(),m_part.rz());
	}
	  
	  
	tmp2=rand(); tmp2/=RAND_MAX;
	if(th_ph_f(3)>tmp2){
	  if(!m_part.check_dir(th_ph_f(1),th_ph_f(2),opts.c_thr.value())){
	    break;
	  }
	    
	  if((th_ph_f(1)!=0&&th_ph_f(2)!=0)){
	    if( (m_mask( round(m_part.x()), round(m_part.y()), round(m_part.z())) != 0) ){
	      if(!opts.modeuler.value())
		m_part.jump(th_ph_f(1),th_ph_f(2));
	      else
		{
		  ColumnVector test_th_ph_f;
		  m_part.testjump(th_ph_f(1),th_ph_f(2));
		  test_th_ph_f=vols.sample(m_part.testx(),m_part.testy(),m_part.testz(),m_part.rx(),m_part.ry(),m_part.rz());
		  m_part.jump(test_th_ph_f(1),test_th_ph_f(2));
		}
	    }
	    
	    
	  }
	}
	  
	  
      }
	
    } // Close Step Number Loop
    if(opts.loopcheck.value()){
      m_loopcheck=0;
    }
    bool accept_path=true;
    if(m_passed_flags.size()!=0){
      for(unsigned int i=0; i<m_passed_flags.size();i++)
	if(!m_passed_flags[i])
	  accept_path=false;
    }   
    return accept_path;
  }
  
  
  void Counter::initialise(){
    if(opts.simpleout.value()||opts.matrix1out.value()){
      initialise_path_dist();
    }
    if(opts.s2tout.value()){
      initialise_seedcounts();
    }
    if(opts.matrix1out.value()){
      initialise_matrix1();
    }
    if(opts.matrix2out.value()){
      initialise_matrix2();
    }
    if(opts.maskmatrixout.value()){
      initialise_maskmatrix();
    }
  }
  
  void Counter::initialise_seedcounts(){
    
    volume<int> tmp;
    //vector<int> tmpvec;
    read_masks(m_targetmasknames,opts.targetfile.value());
    m_targflags.resize(m_targetmasknames.size(),0);
    //m_particle_numbers.resize(m_targetmasknames.size());
    //tmpvec.reserve(opts.nparticles.value());
    cout<<"Number of masks "<<m_targetmasknames.size()<<endl;
    //are they initialised to zero?
    for(unsigned int m=0;m<m_targetmasknames.size();m++){
      read_volume(tmp,m_targetmasknames[m]);
      m_targetmasks.push_back(tmp);
      tmp=0;
      m_seedcounts.push_back(tmp);
      //m_particle_numbers.push_back(tmpvec);
    }
  }
  

  void Counter::initialise_matrix1(){
    m_Conrow=0;
    int numseeds=0;
    for(int Wz=m_seeds.minz();Wz<=m_seeds.maxz();Wz++)
      for(int Wy=m_seeds.miny();Wy<=m_seeds.maxy();Wy++)
	for(int Wx=m_seeds.minx();Wx<=m_seeds.maxx();Wx++)
	  if(m_seeds.value(Wx,Wy,Wz)>0)
	    numseeds++;
      
    m_ConMat.reinitialize(numseeds,numseeds,1);
    m_CoordMat.reinitialize(numseeds,3,1);
    int myrow=0;
      
    for(int Wz=m_seeds.minz();Wz<=m_seeds.maxz();Wz++){
      for(int Wy=m_seeds.miny();Wy<=m_seeds.maxy();Wy++){
	for(int Wx=m_seeds.minx();Wx<=m_seeds.maxx();Wx++){
	  if(m_seeds(Wx,Wy,Wz)>0){
	    m_CoordMat(myrow,0,0)=Wx;
	    m_CoordMat(myrow,1,0)=Wy;
	    m_CoordMat(myrow,2,0)=Wz;
	    myrow++;
	  }
	}
      }
    }
    
  }
  
  void Counter::initialise_matrix2(){
    m_Conrow2=0;
    read_volume(m_lrmask,opts.lrmask.value());
    m_beenhere2.reinitialize(m_lrmask.xsize(),m_lrmask.ysize(),m_lrmask.zsize());
    m_lrdim.ReSize(3);
    m_lrdim<<m_lrmask.xdim()<<m_lrmask.ydim()<<m_lrmask.zdim();
    int numseeds=0,numnz=0;
    for(int Wz=m_seeds.minz();Wz<=m_seeds.maxz();Wz++)
      for(int Wy=m_seeds.miny();Wy<=m_seeds.maxy();Wy++)
	for(int Wx=m_seeds.minx();Wx<=m_seeds.maxx();Wx++)
	  if(m_seeds.value(Wx,Wy,Wz)>0)
	    numseeds++;
    
    for(int Wz=m_lrmask.minz();Wz<=m_lrmask.maxz();Wz++)
      for(int Wy=m_lrmask.miny();Wy<=m_lrmask.maxy();Wy++)
	for(int Wx=m_lrmask.minx();Wx<=m_lrmask.maxx();Wx++)
	  if(m_lrmask.value(Wx,Wy,Wz)>0)
	    numnz++;
    
    
    if(numnz> pow(2,(float)sizeof(short)*8-1)-1){
      cerr<<"Output matrix too big for AVW - stopping."<<endl;
      cerr<<" Remember - you can store your tracts in "<<endl;
      cerr<<" low res even if you want your seeds in high res"<<endl;
      cerr<<" Just subsample the structural space mask"<<endl;
      cerr<<" Although, it must stay in line with the seeds"<<endl;
      exit(-1);
    }
    m_ConMat2.reinitialize(numseeds,numnz,1);
    m_CoordMat2.reinitialize(numseeds,3,1);
    m_CoordMat_tract2.reinitialize(numnz,3,1);
    
    Matrix tempy(numnz,1);
    for(int i=1;i<=numnz;i++){tempy(i,1)=i-1;}
    m_lookup2.addvolume(m_lrmask);
    m_lookup2.setmatrix(tempy.t(),m_lrmask);
      
    int mytrow=0;
    for(int Wz=m_lrmask.minz();Wz<=m_lrmask.maxz();Wz++)
      for(int Wy=m_lrmask.miny();Wy<=m_lrmask.maxy();Wy++)
	for(int Wx=m_lrmask.minx();Wx<=m_lrmask.maxx();Wx++)
	  if(m_lrmask(Wx,Wy,Wz)>0){
	    m_CoordMat_tract2(mytrow,0,0)=Wx;
	    m_CoordMat_tract2(mytrow,1,0)=Wy;
	    m_CoordMat_tract2(mytrow,2,0)=Wz;
	    mytrow++;
	  }
      
    int myrow=0;
    for(int Wz=m_seeds.minz();Wz<=m_seeds.maxz();Wz++)
      for(int Wy=m_seeds.miny();Wy<=m_seeds.maxy();Wy++)
	for(int Wx=m_seeds.minx();Wx<=m_seeds.maxx();Wx++)
	  if(m_seeds(Wx,Wy,Wz)>0){
	    m_CoordMat2(myrow,0,0)=Wx;
	    m_CoordMat2(myrow,1,0)=Wy;
	    m_CoordMat2(myrow,2,0)=Wz;
	    myrow++;
	  }
      
      
  }
  
  void Counter::count_streamline(){
    if(opts.simpleout.value()||opts.matrix1out.value()){
      update_pathdist();
    }
    if(opts.s2tout.value()){
      update_seedcounts();
    }
    if(opts.matrix2out.value()){
      update_matrix2_row();
    }
    if(opts.maskmatrixout.value()){
      update_maskmatrix();
    }
  }
  
  void Counter::count_seed(){
    if(opts.matrix1out.value()){
      update_matrix1();
    }
    if(opts.matrix2out.value()){
      next_matrix2_row();
    }
  }
  
    
  void Counter::clear_streamline(const bool& forwardflag,const bool& backwardflag){
    if(opts.simpleout.value()||opts.matrix1out.value()){
      reset_beenhere(forwardflag,backwardflag);
    }
    if(opts.s2tout.value()){
      reset_targetflags();
    }
    if(opts.matrix2out.value()){
      reset_beenhere2(forwardflag,backwardflag);
    }
    if(opts.maskmatrixout.value()){
      //Do whatever it it you have to do!!
    }
  }
  
  void Counter::update_pathdist(){
    const vector<ColumnVector>& path=m_stline.get_path_ref();
    for(unsigned int i=0;i<path.size();i++){
      int x_s=int(round(float(path[i](1)))),y_s=int(round(float(path[i](2)))),z_s=int(round(float(path[i](3))));
      if(m_beenhere(x_s,y_s,z_s)==0){
	m_prob(x_s,y_s,z_s)+=1;
	m_beenhere(x_s,y_s,z_s)=1;
      }
    }
    
  }

  void Counter::reset_beenhere(const bool& forwardflag,const bool& backwardflag){
    if(forwardflag){
      for(unsigned int i=0;i<m_path.size();i++){
	int x_s=int(round(float(m_path[i](1)))),y_s=int(round(float(m_path[i](2)))),z_s=int(round(float(m_path[i](3))));
	m_beenhere(x_s,y_s,z_s)=0;
      }
    }
    if(backwardflag){
      const vector<ColumnVector>& path=m_stline.get_path_ref();
      for(unsigned int i=0;i<path.size();i++){
	int x_s=int(round(float(path[i](1)))),y_s=int(round(float(path[i](2)))),z_s=int(round(float(path[i](3))));
	m_beenhere(x_s,y_s,z_s)=0;
      }
    }
  }
  
  
  void Counter::update_seedcounts(){
    const vector<ColumnVector>& path=m_stline.get_path_ref();
    int xseedvox=int(round(m_stline.get_x_seed()));
    int yseedvox=int(round(m_stline.get_y_seed()));
    int zseedvox=int(round(m_stline.get_z_seed()));
    for(unsigned int i=0;i<path.size();i++){
      int x_s=int(round(float(path[i](1)))),y_s=int(round(float(path[i](2)))),z_s=int(round(float(path[i](3))));
      for(unsigned int m=0;m<m_targetmasknames.size();m++){
	if(m_targetmasks[m](x_s,y_s,z_s)>0 && m_targflags[m]==0){
	  m_seedcounts[m](xseedvox,yseedvox,zseedvox)=m_seedcounts[m](xseedvox,yseedvox,zseedvox)+1;
	  m_targflags[m]=1;
	  //m_particle_numbers[m].push_back(particle_number);
	}
      }
    }
    
  }
  

  
  void Counter::update_matrix1(){
    //after each particle, update_pathdist(), only run this after each voxel
    int Concol=0;
    for(int Wz=m_prob.minz();Wz<=m_prob.maxz();Wz++){
      for(int Wy=m_prob.miny();Wy<=m_prob.maxy();Wy++){
	for(int Wx=m_prob.minx();Wx<=m_prob.maxx();Wx++){
	  if(m_seeds(Wx,Wy,Wz)>0){
	    if(m_prob(Wx,Wy,Wz)>0){
	      m_ConMat(m_Conrow,Concol,0)=m_prob(Wx,Wy,Wz);
	    }
	    Concol++;
	  }
	  m_prob(Wx,Wy,Wz)=0;
	  
	}
      }
    }
    
    m_Conrow++;
  }
  
  void Counter::update_matrix2_row(){
    //run this one every streamline - not every voxel..
    const vector<ColumnVector>& path=m_stline.get_path_ref();
    for(unsigned int i=0;i<path.size();i++){
      ColumnVector xyz_seeds=path[i];
      ColumnVector xyz_lr=vox_to_vox(xyz_seeds,m_seedsdim,m_lrdim,m_I);
      int x_lr=int(round(float(xyz_lr(1)))),y_lr=int(round(float(xyz_lr(2)))),z_lr=int(round(float(xyz_lr(3))));
      int Concol2=m_lookup2(x_lr,y_lr,z_lr,0);
      if(Concol2!=0){
	if(m_beenhere2(x_lr,y_lr,z_lr)==0){
	  m_ConMat2(m_Conrow2,Concol2,0)+=1;
	  m_beenhere2(x_lr,y_lr,z_lr)=1;
	}
      }
      
    }
    
  }
  
  
  
  void Counter::reset_beenhere2(const bool& forwardflag,const bool& backwardflag){
    if(forwardflag){
      for(unsigned int i=0;i<m_path.size();i++){
	ColumnVector xyz_seeds=m_path[i];
	ColumnVector xyz_lr=vox_to_vox(xyz_seeds,m_seedsdim,m_lrdim,m_I);
	int x_lr=int(round(float(xyz_lr(1)))),y_lr=int(round(float(xyz_lr(2)))),z_lr=int(round(float(xyz_lr(3))));
	m_beenhere2(x_lr,y_lr,z_lr)=0;
      }
    }
    if(backwardflag){
      const vector<ColumnVector>& path=m_stline.get_path_ref();
      for(unsigned int i=0;i<path.size();i++){
	ColumnVector xyz_seeds=path[i];
	ColumnVector xyz_lr=vox_to_vox(xyz_seeds,m_seedsdim,m_lrdim,m_I);
	int x_lr=int(round(float(xyz_lr(1)))),y_lr=int(round(float(xyz_lr(2)))),z_lr=int(round(float(xyz_lr(3))));
	m_beenhere2(x_lr,y_lr,z_lr)=0;
      }  
    }
    
  }

  void Counter::save(){
    if(opts.simpleout.value()){
      save_pathdist();
    }
    if(opts.s2tout.value()){
      save_seedcounts();
    }
    if(opts.matrix1out.value()){
      save_matrix1();
    }
    if(opts.matrix2out.value()){
      save_matrix2();
    }
    if(opts.maskmatrixout.value()){
      save_maskmatrix();
    }
    
  }
  
  void Counter::save_pathdist(){  
    save_volume(m_prob,logger.appendDir("fdt_paths"));
  }
  
  void Counter::save_pathdist(string add){  //for simple mode
    string thisout=opts.outfile.value();
    make_basename(thisout);
    thisout+=add;
    save_volume(m_prob,thisout);
  }

  void Counter::save_seedcounts(){
    for(unsigned int m=0;m<m_targetmasknames.size();m++){
      string tmpname=m_targetmasknames[m];
      
      int pos=tmpname.find("/",0);
      int lastpos=pos;
      
      while(pos>=0){
	lastpos=pos;
	pos=tmpname.find("/",pos);
	// replace / with _
	tmpname[pos]='_';
      }
      
      //only take things after the last pos
      tmpname=tmpname.substr(lastpos+1,tmpname.length()-lastpos-1);
      
      save_volume(m_seedcounts[m],logger.appendDir("seeds_to_"+tmpname));
    }
  }
  
  void Counter::save_matrix1(){
    save_volume(m_ConMat,logger.appendDir("fdt_matrix1"));
    save_volume(m_CoordMat,logger.appendDir("coords_for_fdt_matrix1"));
  }

  void Counter::save_matrix2(){
    save_volume(m_ConMat2,logger.appendDir("fdt_matrix2"));
    save_volume(m_CoordMat2,logger.appendDir("coords_for_fdt_matrix2"));
    save_volume(m_CoordMat_tract2,logger.appendDir("tract_space_coords_for_fdt_matrix2"));
    save_volume4D(m_lookup2,logger.appendDir("lookup_tractspace_fdt_matrix2"));
  
  }
  
  void Seedmanager::run(const float& x,const float& y,const float& z,bool onewayonly, int fibst){
    //onewayonly for mesh things..
    cout <<x<<" "<<y<<" "<<z<<endl;
    if(fibst == -1){
      fibst=m_seeds(int(round(x)),int(round(y)),int(round(z)))-1;//fibre to start with is taken from seed volume..
    }
    if(opts.randfib.value()){
      float tmp=rand()/RAND_MAX;
      if(tmp>0.5)
	fibst=0;
      else
	fibst=1;// fix this for > 2 fibres
    }
    
    for(int p=0;p<opts.nparticles.value();p++){
      if(opts.verbose.value()>1)
	logger.setLogFile("particle"+num2str(p));
      
      m_stline.reset();
      bool forwardflag=false,backwardflag=false;
      if(!onewayonly){
	if(m_stline.streamline(x,y,z,m_seeddims,fibst)){ //returns whether to count the streamline or not
	  forwardflag=true;
	  m_counter.store_path();
	  m_counter.count_streamline();
	}
	m_stline.reverse();
      }
      if(m_stline.streamline(x,y,z,m_seeddims,fibst)){
	backwardflag=true;
	m_counter.count_streamline();
      }
     
      m_counter.clear_streamline(forwardflag,backwardflag); 
    }
    m_counter.count_seed();
    
    
  }



}
  
  

