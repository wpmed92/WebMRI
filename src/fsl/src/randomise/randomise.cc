/*  randomise.cc

    Tim Behrens & Steve Smith & Matthew Webster (FMRIB) & Tom Nichols (UMich)

    Copyright (C) 2004-2006 University of Oxford  */

/*  CCOPYRIGHT  */

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "libprob.h"
#include "ranopts.h"
#include <algorithm>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

using namespace RANDOMISE;




bool check_dims(const volume4D<float>& data,const volume<float>& mask,const  Matrix& dm,const  Matrix& confounds,const  Matrix& tc,const  Matrix& fc){
  bool  ret=true;
  if (!samesize(data[0],mask)) { 
    cerr << "mask dimensions do not match input data dimensions!" << endl;
    ret=false;
  }
  if (dm.Nrows()!=data.tsize()){
      cerr << "number of rows in design matrix doesn't match number of \"time points\" in input data!" << endl;
      ret=false;
  }
  if(confounds.Nrows()!=0){
    if(confounds.Nrows()!=dm.Nrows()){
      cerr << "number of rows in confound matrix doesn't match number of rows in design matrix!" << endl;
      ret=false;
    }
  }
  if (tc.Ncols()!=dm.Ncols()){
    cerr << "number of columns in t-contrast matrix doesn't match number of columns in design matrix!" << endl;
    ret=false;
  }
  if(fc.Ncols() !=0){
    if (fc.Ncols()!=tc.Nrows()){
      cerr << "number of columns in f-contrast matrix doesn't match number of rows in t-contrast matrix!" << endl;
      ret=false;
    }
  }
  return ret;
}


bool next_mult_vec(ColumnVector& mult){
  bool finish=false;
  int n=mult.Nrows();

  if (mult(1)==0)
    mult=-1; // initialise
  else {
    while(!finish){
      if(mult(n)==-1){
	mult(n)=1;
	if (n<mult.Nrows()) mult.Rows(n+1,mult.Nrows())=-1;
	finish=true;
      }
      else
	n--;
    }
  }

  return (!(mult.Minimum()==1));
}



void rand_mult_vec(ColumnVector& mult){

  for(int i=1;i<=mult.Nrows();i++){
    float tmp=rand();
    tmp/=RAND_MAX;
    if(tmp > 0.5)
      mult(i)=1;
    else
      mult(i)=-1;
  }
  
}


void dm_mult(Matrix& dm,const ColumnVector& mult){
  for(int i=1;i<=dm.Nrows();i++){
    dm.Row(i)*=mult(i);
  }
}


void dm_mult(const Matrix& dm, Matrix& dm_new, const ColumnVector& mult){
  for(int i=1;i<=dm.Nrows();i++){
    dm_new.Row(i)=dm.Row(i)*mult(i);
  }
}



void dm_permute(const Matrix& dm_in, Matrix& dm_out, const ColumnVector& r){

  if(r.Nrows() != dm_in.Nrows()){
    cerr<<"permutation vector has wrong number of elements"<<endl;
    exit(-1);
  }
  
  if ( (dm_out.Nrows()!=dm_in.Nrows()) || (dm_out.Ncols()!=dm_in.Ncols()) ){
    dm_out.ReSize(dm_in.Nrows(),dm_in.Ncols());
  }

  for(int row=1;row<=r.Nrows();row++){
    dm_out.Row(row) << dm_in.Row(int(r(row)));
  }

}




void rand_perm_vec(ColumnVector& r){
  int n=r.Nrows();
  vector<pair<float,int> > tmpvec(n);
  
  
  for(int i=0;i<n;i++){
    tmpvec[i].first=rand();
    tmpvec[i].second=i+1;
  }
  
  sort(tmpvec.begin(),tmpvec.end());
  for(int i=1;i<=n;i++){
    r(i)=tmpvec[i-1].second;
  }
  
}


void next_perm_vec(ColumnVector& r){

  if(r(1)<0){
    // initialise and return
    r = -r;
  }
  else{
     bool move=false;
    int n=r.Nrows();
    int j=1,i=0;
    ColumnVector tmp;
    while(!move){
      
      if(r(n-i)>r(n-j)){
	tmp=r.Rows(n-j,n-i-1);
	r(n-j)=r(n-i);
	r.Rows(n-j+1,n-i)=tmp;
	tmp=r.Rows(n-j+1,r.Nrows());
	SortAscending(tmp);
	r.Rows(n-j+1,r.Nrows())<<tmp;            
	move=true;
      }
      else
	i++;

      if(i==j){
	j++;i=0;
      }
    }

  }
}




void make_dm_labels(const Matrix& dm,const ColumnVector& r, ColumnVector& labels){
  if(labels.Nrows() != r.Nrows()){
    labels.ReSize(r.Nrows()); 
  }
  ColumnVector labels_orig(r.Nrows());
  vector<RowVector> label_lu;
  bool found_it=false;
  for(int i=1;i<=dm.Nrows();i++){
    found_it=false;
    for(unsigned int l=0;l<label_lu.size();l++){
      if(dm.Row(i)==label_lu[l]){
	labels_orig(i)=l+1;
	found_it=true;
      }
    }
    if(!found_it){
      label_lu.push_back(dm.Row(i));
      labels_orig(i)=label_lu.size();
    }
  }
  
  // Now permute labels 
  for(int i=1; i<=labels.Nrows(); i++) labels(i)=labels_orig(int(r(i)));
}



bool check_perm(vector<ColumnVector>& oldperms,const ColumnVector& newperm){
  bool isin=false;
  for(int i=oldperms.size()-1; i>=0; i--){
    if(newperm==oldperms[i]){
      isin=true;
//       ColumnVector tmp=oldperms[i];
//       for(int j=1;j<=tmp.Nrows();j++)
//       {
//         cout << tmp(j) << " " << newperm(j) << endl;
//       }
//       cout << endl;
      break;

    }
  } 
  if(!isin){
    oldperms.push_back(newperm);
  }
  return !isin;
}

void ols_var_sm(const Matrix& data,const Matrix& des,const Matrix& tc, Matrix& cope,Matrix& varcope,const volume<float>& mask,const volume<float>& mask_sm,float sigma_mm){
  // ols_var_sm
  // data is t x v
  // des is t x ev (design matrix)
  // tc is cons x ev (contrast matrix)
  // cope and varcope will be cons x v
  // but will be resized if they are wrong
  // hence may be passed in uninitialised
  // TB 2004

  if(data.Nrows() != des.Nrows()){
    cerr <<"RANDOMISE::ols_var_sm - data and design have different number of time points"<<endl;
    exit(-1);
  }
  if(des.Ncols() != tc.Ncols()){
    cerr <<"RANDOMISE::ols_var_sm - design and contrast matrix have different number of EVs"<<endl;
    exit(-1);
  }  
  volume4D<float> sigsqvol;
  
  Matrix pdes = pinv(des);
  Matrix prevar=diag(tc*pdes*pdes.t()*tc.t());
  Matrix R=Identity(des.Nrows())-des*pdes;
  float tR=R.Trace();
  Matrix pe=pdes*data;
  cope=tc*pe;
  Matrix res=data-des*pe;
  Matrix sigsq=sum(SP(res,res))/tR;
  sigsqvol.setmatrix(sigsq,mask);
  sigsqvol[0]=smooth(sigsqvol[0],sigma_mm);
  sigsqvol[0]/=mask_sm;
  sigsq=sigsqvol.matrix(mask);
  varcope=prevar*sigsq;

}


double compute_nperms(const Matrix& dm){
  
  ColumnVector labels;
  ColumnVector tmp(dm.Nrows());
  for(int i=1;i<=tmp.Nrows();i++){tmp(i)=i;};
  
    make_dm_labels(dm,tmp,labels);

    int num_labels=int(labels.MaximumAbsoluteValue());
    ColumnVector label_counts(num_labels);
    label_counts=0;
    for(int i=1; i<=labels.Nrows(); i++){
      label_counts(int(labels(i)))+=1;
    }
    float yo = lgam(dm.Nrows()+1);
    for(int i=1; i<=num_labels; i++)
      yo -= lgam(label_counts(i)+1);
    
    return std::floor(exp(yo)+0.5);

}

void delete_clusters(volume<int>& cvol,const ColumnVector& clustersizes, int size_thresh){
  for(int z=0;z<cvol.zsize();z++){
    for(int y=0;y<cvol.ysize();y++){
      for(int x=0;x<cvol.xsize();x++){
	if(cvol(x,y,z)!=0){
	  if(clustersizes(cvol(x,y,z)) < size_thresh ){
	    cvol(x,y,z)=0;
	  }
	}
      }
    }
  }
}


int main(int argc,char *argv[]){

  Log& logger = LogSingleton::getInstance();
  ranopts& opts = ranopts::getInstance();
  opts.parse_command_line(argc,argv,logger);

  Matrix dm, confound, tc, fc;
  double n_exhaust_perms;

  if(!opts.one_samp.value()){

    // read dm, tc and fc, and estimate max-perms if required
    dm=read_vest(opts.dm_file.value());

    n_exhaust_perms=compute_nperms(dm);

    if (opts.verbose.value() || opts.how_many_perms.value())
      cout << n_exhaust_perms << " permutations required for exhaustive test" << endl;
    if (opts.how_many_perms.value())
      return 0;

    tc=read_vest(opts.tc_file.value());

    //    if (opts.fc_file.value()!=""){
    //  fc=read_vest(opts.fc_file.value());
    //}
  }
  else {

    // group-mean only; create necessary dm and tc, estimate max-perms if required
    volume4D<float> data;
    read_volume4D_hdr_only(data, opts.in_fileroot.value());

    dm.ReSize(data.tsize(),1);
    dm=1;
    tc.ReSize(1,1);
    tc=1;

    n_exhaust_perms=std::pow(2.0,data.tsize());

    if (opts.verbose.value() || opts.how_many_perms.value())
      cout << n_exhaust_perms << " inversion permutations required for exhaustive test" << endl;
    if (opts.how_many_perms.value())
      return 0;
  }

  
  if(opts.confound_file.value()!=""){
    confound=read_vest(opts.confound_file.value());
  }

  volume<float> mask;
  volume<int> tmpvol;

  // read data and convert to matrix, optionally using supplied mask  
  Matrix datam;
  {// scope in which data exists in volume4D format
    
    volume4D<float> data;
    read_volume4D(data,opts.in_fileroot.value());
    
    if (opts.maskname.value()!=""){
      read_volume(mask,opts.maskname.value());
      mask.binarise(0.0001);
    } else {
      mask = meanvol(data);
      mask.binarise(0.0001);
    }    
    if(!check_dims(data,mask,dm,confound,tc,fc)){
      exit(-1);
    }
    
    datam=data.matrix(mask);   

    if (opts.demean_data.value())
      datam=remmean(datam);

    if (opts.verbose.value())
      cout << "data loaded" << endl;
  }

#ifdef grot  
  volume<float> mask;
  read_volume(mask,opts.maskname.value());
  volume<int> tmpvol;
  // read data and convert to matrix, optionally using supplied mask  
  Matrix datam;
  {// scope in which data exists in volume4D format
    volume<short> smask;
    read_volume(smask,opts.maskname.value());
    
    volume4D<short> data;
    read_volume4D(data,opts.in_fileroot.value());
    
    datam=data.matrix(smask);   

    if (opts.demean_data.value())
      datam=remmean(datam);

    if (opts.verbose.value())
      cout << "data loaded" << endl;
  }
#endif


  // prepare smoothed mask for use (as a convolution renormaliser) in variance smoothing if required
  volume<float> mask_sm;
  if(opts.var_sm_sig.value()>0){
    mask_sm=smooth(mask,opts.var_sm_sig.value());
  }

  Matrix cope, tstat, tstat_orig, fstat, fstat_orig, varcope;

  //How  many permutations will we do
  int n_perms;
  bool exhaustive=0;
  if(opts.n_perm.value()==0){
    n_perms=int(n_exhaust_perms);
    exhaustive=1;
    if (opts.verbose.value()) cout<<"will do all "<<n_perms<<" unique permutations"<<endl;
  }
  else if(opts.n_perm.value()>=n_exhaust_perms){
    n_perms=int(n_exhaust_perms);
    exhaustive=1;
    if (opts.verbose.value()) cout<<"will do "<<n_perms<<" permutations"<<endl;
  }
  else{
    n_perms=opts.n_perm.value();
    exhaustive=0;
    if (opts.verbose.value()) cout<<"doing "<<n_perms<<" permutations as requested"<<endl;
  }
  // containers for different inference distribution
  Matrix maxdist, maxdistC, maxdistCmass, vox_numbigger;
  volume4D<float> tstat4D, tmp_tstat4D;

  
  // resize the containers for the relevant inference distributions
  maxdist.ReSize(tc.Nrows(),n_perms);
  maxdist=0;
  vox_numbigger.ReSize(tc.Nrows(),datam.Ncols()); //number of permuted  which are bigger than original
  vox_numbigger=0;
  if ( opts.cluster_thresh.value()>0 |opts.clustermass_thresh.value()> 0 ) { //cluster thresholding
    maxdistC.ReSize(tc.Nrows(),n_perms);
    maxdistC=0;
    maxdistCmass.ReSize(tc.Nrows(),n_perms);
    maxdistCmass=0;
    tstat4D.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),tc.Nrows());
  }
  


  // remove confound space from data  
  if(opts.confound_file.value()!=""){
    datam=(Identity(confound.Nrows())-confound*pinv(confound))*datam;
  }
  //Things needed to generate the permutation or multiplication vectors etc.
  Matrix dmperm=dm;
  ColumnVector one_to_n(dm.Nrows());
  for(int i=1;i<=one_to_n.Nrows();i++) one_to_n(i)=i;
  ColumnVector permvec = -one_to_n;            // tell next_perm_vec to initialise on first pass
  ColumnVector multvec(dm.Nrows()); multvec=0; // tell next_mult_vec to initialise on first pass
  ColumnVector labels;
  vector<ColumnVector> oldperms; // store old perm/mult vecs to avoid redundancy.
  ColumnVector clustersizes,clustermasses;  
  bool accept;
  bool correct_perm_done=0;

  for(int perm=1; perm<=n_perms; perm++){

    if (opts.verbose.value()) cout << "starting permutation " << perm << endl;
    // Get next design matrix (one of the 4 techniques depending on options)
    if(opts.one_samp.value()){
      if(exhaustive){
	next_mult_vec(multvec);
      }
      else{
	accept=false;
	while(!accept){
        
	  rand_mult_vec(multvec);
	  make_dm_labels(multvec,one_to_n,labels);
	  accept=check_perm(oldperms,labels);
}
      }

      dm_mult(dm,dmperm,multvec);
    }
    else{
      if(exhaustive){
        accept=false;
        while(!accept){
          //DEBUG
          if (perm == 1) { next_perm_vec(permvec); make_dm_labels(dm,permvec,labels); }
          else
	    {
	      ColumnVector oldlabels=labels;
               next_perm_vec(labels);
                for (int k=1;k<=labels.Nrows();k++)		{
                   if (labels(k)!=oldlabels(k))
  		 {
  		      for (int l=1;l<=labels.Nrows();l++) 
  		      {
                          if (labels(l)!=oldlabels(l) && oldlabels(l)==labels(k) )
 			{
 			  double oldvec = permvec(l);
                           permvec(l)=permvec(k);
                           permvec(k)=oldvec;
                           oldvec=oldlabels(l);
                           oldlabels(l)=oldlabels(k);
                           oldlabels(k)=oldvec;  
                         }
                        }
 		 }
	       }
	    }
          //end of guaranteed acceptable perms...
          accept=check_perm(oldperms,labels);
        } 
      }
      else{
	accept=false;
	while(!accept){
	  rand_perm_vec(permvec);     
	  make_dm_labels(dm,permvec,labels);
	  accept=check_perm(oldperms,labels);
	}
      }
      dm_permute(dm,dmperm,permvec);
    }

    // make sure that the "real" permutation is used
    if (dm==dmperm)
      correct_perm_done=1;
    if ( (perm==n_perms) && !correct_perm_done)
      dmperm=dm;

    //cout << dmperm << endl;
  
    if(perm==1){ //compute stats for the original case
      if(opts.var_sm_sig.value()==0)
	ols(datam,dm,tc,cope,varcope);   
      else
	ols_var_sm(datam,dm,tc,cope,varcope,mask,mask_sm,opts.var_sm_sig.value());   
      tstat_orig=SD(cope,sqrt(varcope));
    }

    // compute stats for current permutation
    if(opts.var_sm_sig.value()==0)
      ols(datam,dmperm,tc,cope,varcope);   
    else
      ols_var_sm(datam,dmperm,tc,cope,varcope,mask,mask_sm,opts.var_sm_sig.value());   
    tstat=SD(cope,sqrt(varcope));

    
    maxdist.Column(perm)<<max(tstat.t()).t(); // max stat recording

    vox_numbigger+= geqt(tstat,tstat_orig); // voxelwise stats recording

    if ( opts.cluster_thresh.value() > 0 ) { //cluster thresholding
      tstat4D.setmatrix(tstat,mask);
      tstat4D.binarise(opts.cluster_thresh.value());

      for(int t=0;t<tstat4D.tsize();t++){

	connected_components(tstat4D[t],clustersizes,26);
	
	if ( clustersizes.Nrows() > 0 ) {
	
	 
	  // Store Max cluster size for each perm
	  
	  maxdistC(t+1,perm)=int(clustersizes.MaximumAbsoluteValue());
	}
      }
    }
    if ( opts.clustermass_thresh.value() > 0 ) { //cluster mass thresholding
      tstat4D.setmatrix(tstat,mask);
      tmp_tstat4D=tstat4D;
      tstat4D.binarise(opts.clustermass_thresh.value());
      
      for(int t=0;t<tstat4D.tsize();t++){


	tmpvol=connected_components(tstat4D[t],clustersizes,26);
	clustermasses.ReSize(clustersizes.Nrows());
	clustermasses=0;
	
	for(int z=0; z<tmpvol.zsize(); z++)
	  for(int y=0; y<tmpvol.ysize(); y++)
	    for(int x=0; x<tmpvol.xsize(); x++)
	      if(tmpvol(x,y,z)>0)
		clustermasses(tmpvol(x,y,z))=clustermasses(tmpvol(x,y,z))+tmp_tstat4D[t](x,y,z);
	
	if ( clustersizes.Nrows() > 0 ) {
	  // Store Max cluster mass for each perm
	  maxdistCmass(t+1,perm)=int(clustermasses.MaximumAbsoluteValue());
	}
      }
    }
    //End of for loop
  }

  volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
  RowVector tmpvec;
  Matrix this_t_stat(1, tstat_orig.Ncols());

  for(int row=1; row<=tc.Nrows(); row++){

    // save unthresholded tstat
    output.setmatrix(tstat_orig.Row(row),mask);
    save_volume4D(output,opts.out_fileroot.value()+"_tstat"+num2str(row));

    // max stat
    tmpvec=maxdist.Row(row);
    SortAscending(tmpvec);
    this_t_stat=0;
    for(int i=1; i<=tstat_orig.Ncols(); i++)
      for(int j=1; j<=n_perms; j++)
	if (tstat_orig(row,i)>tmpvec(j))
	  this_t_stat(1,i) = float(j)/n_perms;
    output.setmatrix(this_t_stat,mask);
    save_volume4D(output,opts.out_fileroot.value()+"_max_tstat"+num2str(row));

    // vox stat
    this_t_stat=1-vox_numbigger.Row(row)/float(n_perms);
    output.setmatrix(this_t_stat,mask); 
    save_volume4D(output,opts.out_fileroot.value()+"_vox_tstat"+num2str(row));

  }

  if ( opts.cluster_thresh.value() > 0 ) { //cluster thresholding
    tstat4D.setmatrix(tstat_orig,mask);
    tstat4D.binarise(opts.cluster_thresh.value());
    
    for(int row=1; row<=tc.Nrows(); row++){ //for each contrast
      
      //       int cdf=0;
      //       for(int j=1;j<=cluster_counts.Ncols();j++){
      // 	//count smallest (1-p)*100 % clustersizes
      // 	cdf+=int(cluster_counts(row,j));
      // 	if(cdf>=int((1-0.05)*n_perms)){
      // 	  // find clust sizes in real thing
      // 	  tmpvol=connected_components(tstat4D[row],clustersizes,26);
      // 	  delete_clusters(tmpvol,clustersizes,j);
      // 	  tmpvol.binarise(1); //check if this is inclusive
      // 	  output[row]=NEWIMAGE::mask_volume(tstat4D[row],tmpvol); 
      // 	  break; //out of this contrast.
      // 	} 
      //       }

      tmpvol=connected_components(tstat4D[row-1],clustersizes,26);
      output=0;
      tmpvec=maxdistC.Row(row);
      SortAscending(tmpvec);
      for(int i=1; i<=clustersizes.Nrows(); i++)
	if (clustersizes(i)>0) {
	  for(int j=1; j<=n_perms; j++)
	    {
	      if (clustersizes(i)>tmpvec(j)){
	      for(int z=0;z<output.zsize();z++) // is there a nicer way to do this with newimage??
		for(int y=0;y<output.ysize();y++)
		  for(int x=0;x<output.xsize();x++)
		    if (tmpvol(x,y,z)==i)  { output(x,y,z,0) = float(j)/n_perms; }
	      }}
	}
      save_volume4D(output,opts.out_fileroot.value()+"_maxc_tstat"+num2str(row));

    }
  }




 
 if ( opts.clustermass_thresh.value() > 0 ) { //cluster mass thresholding
    tstat4D.setmatrix(tstat_orig,mask);
    tmp_tstat4D=tstat4D;
    tstat4D.binarise(opts.clustermass_thresh.value());
    
    for(int row=1; row<=tc.Nrows(); row++){ //for each contrast
     
      tmpvol=connected_components(tstat4D[row-1],clustersizes,26);
      clustermasses.ReSize(clustersizes.Nrows());
      clustermasses=0;
      
      for(int z=0; z<tmpvol.zsize(); z++)
	for(int y=0; y<tmpvol.ysize(); y++)
	  for(int x=0; x<tmpvol.xsize(); x++)
	    if(tmpvol(x,y,z)>0)
	      clustermasses(tmpvol(x,y,z))=clustermasses(tmpvol(x,y,z))+tmp_tstat4D[row-1](x,y,z);
      
      
      output=0;
      tmpvec=maxdistCmass.Row(row);
      SortAscending(tmpvec);
      for(int i=1; i<=clustermasses.Nrows(); i++)
	if (clustermasses(i)>0) {
	  for(int j=1; j<=n_perms; j++)
	    if (clustermasses(i)>tmpvec(j))
	      for(int z=0;z<output.zsize();z++) // is there a nicer way to do this with newimage??
		for(int y=0;y<output.ysize();y++)
		  for(int x=0;x<output.xsize();x++)
		    if (tmpvol(x,y,z)==i)
		      output(x,y,z,0) = float(j)/n_perms;
	}
      save_volume4D(output,opts.out_fileroot.value()+"_maxcmass_tstat"+num2str(row));

    }
  }


  return 0;
}

