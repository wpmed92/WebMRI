/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
              independent components
    
    melreport.cc - report generation

    Christian F. Beckmann, FMRIB Image Analysis Group
    
    Copyright (C) 1999-2004 University of Oxford */

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

#include "newimage/newimageall.h"
#include "utils/log.h"
#include "melreport.h"

namespace Melodic{
  
  void MelodicReport::IC_rep(MelGMix &mmodel, int cnum, int dim, Matrix ICstats)
  {
    if( bool(opts.genreport.value()) ){
      addlink(mmodel.get_prefix()+".html",num2str(cnum));
      IChtml.setDir(report.getDir(),mmodel.get_prefix()+".html");
      
      {//start IC page

      IChtml << "<HTML> " << endl
	     << "<TITLE>MELODIC Component " << num2str(cnum)
	     << "</TITLE>" << endl
	     << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	     << "/doc/images/fsl-bg.jpg\">" << endl 
	     << "<hr><CENTER><H1>MELODIC Component " << num2str(cnum)
	     << "</H1>"<< endl;
      
      if(cnum>1)
	IChtml << "<a href=\"" << string("IC_")+num2str(cnum-1)
	       <<".html\">previous</a>&nbsp;-&nbsp;";
      
      IChtml << "<a href=\"00index.html\">&nbsp;index&nbsp;</a>" ;

      if(cnum<dim)
	IChtml << "&nbsp;-&nbsp;<a href=\"" << string("IC_")+num2str(cnum+1)
	       <<".html\">next</a><p>";

      IChtml << "<hr><p>" << endl;
      }

      {//output IC stats
	if(ICstats.Storage()>0&&ICstats.Nrows()>=cnum){
	  IChtml << ICstats(cnum,1) << " % of explained variance";
	  if(ICstats.Ncols()>1)
	    IChtml << "; &nbsp; &nbsp; " << ICstats(cnum,2) << " % of total variance";
	  if(ICstats.Ncols()>2){
	    IChtml << "<p>" <<endl;
	    IChtml << " &nbsp; &nbsp; " << ICstats(cnum,3) << " % signal change (pos peak voxel); &nbsp; &nbsp;" << ICstats(cnum,4) << "% signal change (peak neg. voxel)" << endl ;}
	  IChtml << "<hr><p>" << endl;
	}
      }

      volume4D<float> tempVol; 
      
      if(mmodel.get_threshmaps().Storage()>0&&
	 (mmodel.get_threshmaps().Ncols() == mmodel.get_data().Ncols()))
      {//Output thresholded IC map

	
	tempVol.setmatrix(mmodel.get_threshmaps().Row(1),melodat.get_mask());
 	volume<float> map1;
	volume<float> map2;
	map1 = threshold(tempVol[0],float(0.0), tempVol[0].max());
	map2 = threshold(tempVol[0],tempVol[0].min(), float(0.0));
	// 	save_volume(map,report.appendDir(mmodel.get_prefix()+"_thresh")
 	//	    ,melodat.tempInfo);
	volume<float> newvol;
 	// for(int ctr=1; ctr<500; ctr++){
// 	  cerr << ctr << " ";
	miscpic newpic;
	
	float map1min = std::max((map1 + binarise(tempVol[0],tempVol[0].min(), 
		        float(0.0)) * map1.max()).min(),float(0.01));
	float map1max = std::max(map1.max(),float(0.01));
	float map2min = std::min(map2.min(),float(-0.01));
	float map2max = std::min((map2 + binarise(tempVol[0],float(0.0), 
			tempVol[0].max()) * map2.min()).max(),float(-0.01));
	

	newpic.overlay(newvol, melodat.get_mean(), map1, map2, 
		       melodat.get_mean().percentile(0.02),
		       melodat.get_mean().percentile(0.98),
		       map1min, map1max, map2max, map2min, 
		       0, 0, &melodat.tempInfo);
	
 	char instr[10000];
	
	//save_volume(newvol,report.appendDir(mmodel.get_prefix()+"rendered"),
	//      melodat.tempInfo);
 	sprintf(instr," ");
 	strcat(instr," -s 2");
 	strcat(instr," -A 950 ");
 	strcat(instr,string(report.appendDir(mmodel.get_prefix()+
 					     "_thresh.png")).c_str());
 	newpic.set_title(string("Component No. "+num2str(cnum)+
 				" - thresholded IC map  ")+
			 mmodel.get_infstr(1));
 	newpic.set_cbar(string("ysb"));
	//cerr << instr << endl;


	
	if(map1.max()-map1.min()>0.01)
	  newpic.slicer(newvol, instr, &melodat.tempInfo); 
	else
	  newpic.slicer(map1, instr, &melodat.tempInfo);

	}
     //  }
//       exit(2);

      IChtml << "<a href=\""+mmodel.get_prefix()+"_MM.html\">";
      IChtml << "<img BORDER=0 SRC=\""+mmodel.get_prefix()
	+"_thresh.png\" ALT=\"MMfit\"></A><p>" << endl;
      

      {//plot time course
	miscplot newplot;
      
	if(opts.tr.value()>0.0)
	  newplot.timeseries(melodat.get_mix().Column(cnum).t(),
			     report.appendDir(string("t")+num2str(cnum)+".png"),
			     string("Timecourse (in seconds); TR = ")+
			     float2str(opts.tr.value(),0,2,0)+" s", 
			     opts.tr.value(),150,4,1);
	else
	  newplot.timeseries(melodat.get_mix().Column(cnum).t(),
			     report.appendDir(string("t")+num2str(cnum)+".png"),
			     string("Timecourse (in TRs)"));
	write_ascii_matrix(report.appendDir(string("t")
					    +num2str(cnum)+".txt"),
			   melodat.get_mix().Column(cnum));
	IChtml << "<A HREF=\"" << string("t")
	  +num2str(cnum)+".txt" << "\"> ";
	IChtml << "<img BORDER=0 SRC=\"" 
	  +string("t")+num2str(cnum)+".png\"></A><p>" << endl;
      }//time series plot
      
      {//plot frequency 
	miscplot newplot;
	RowVector empty(1);
	empty = 0.0;
	int fact = int(std::pow(10.0,int(std::log10(float(melodat.get_mix().Nrows())))));
      
	if(opts.tr.value()>0.0)
	  newplot.timeseries(empty | melodat.get_fmix().Column(cnum).t(),
			     report.appendDir(string("f")+
					      num2str(cnum)+".png"),
			     string("FFT of timecourse (in Hz / ") +num2str(fact)+")",
			     fact/(opts.tr.value()*melodat.get_mix().Nrows()), 150,
			     0,2);
	else
	  newplot.timeseries(melodat.get_fmix().Column(cnum).t(),
			     report.appendDir(string("f")+
					      num2str(cnum)+".png"),
			     string(string("FFT of timecourse (in cycles); ")
				  +"frequency(Hz)=cycles/("
				  +num2str(melodat.get_mix().Nrows())
				  +"* TR); period(s)=("
				  +num2str(melodat.get_mix().Nrows())
				  +"* TR)/cycles"));
	write_ascii_matrix(report.appendDir(string("f")
					    +num2str(cnum)+".txt"),
			   melodat.get_fmix().Column(cnum));
	IChtml << "<A HREF=\"" << string("f")
	  +num2str(cnum)+".txt" << "\"> ";
	IChtml << "<img BORDER=0 SRC=\"" 
	  +string("f")+num2str(cnum)+".png\"></A><p>" << endl;
      }//frequency plot
      
      
      if(mmodel.get_threshmaps().Storage()>0&&
	 (mmodel.get_threshmaps().Ncols() == mmodel.get_data().Ncols())&&
	 (mmodel.get_threshmaps().Nrows()>1))
	{//Output other thresholded IC map
	  
	  for(int tctr=2; tctr<=mmodel.get_threshmaps().Nrows(); tctr++){
	    tempVol.setmatrix(mmodel.get_threshmaps().Row(tctr),melodat.get_mask());
	    volume<float> map1;
	    volume<float> map2;
	    map1 = threshold(tempVol[0],float(0.0), tempVol[0].max());
	    map2 = threshold(tempVol[0],tempVol[0].min(), float(0.0));
	    // 	save_volume(map,report.appendDir(mmodel.get_prefix()+"_thresh")
	    //	    ,melodat.tempInfo);
	    volume<float> newvol; 
	    miscpic newpic;

	    float map1min = (map1 + binarise(tempVol[0],tempVol[0].min(), 
					     float(0.0)) * map1.max()).min();
	    float map1max = map1.max();
	    float map2min = map2.min();
	    float map2max = (map2 + binarise(tempVol[0],float(0.0), 
					     tempVol[0].max()) * map2.min()).max();
	    
	    //cerr << endl << map1min << " " << map1max << endl
	    //  << map2min << " " << map2max << endl;

	    //	    if(map1.max()-map1.min()>0.01)
	      newpic.overlay(newvol, melodat.get_mean(), map1, map2, 
			     melodat.get_mean().percentile(0.02),
			     melodat.get_mean().percentile(0.98),
			     map1min, map1max, map2max, map2min, 
			     0, 0, &melodat.tempInfo);
	    

	    char instr[10000];
	    
	    //save_volume(newvol,report.appendDir(mmodel.get_prefix()+"rendered"),
	    //      melodat.tempInfo);
	    sprintf(instr," ");
	    strcat(instr," -s 2");
	    strcat(instr," -A 950 ");
	    strcat(instr,string(report.appendDir(mmodel.get_prefix()+"_thresh"+
						 num2str(tctr)+".png")).c_str());
	    newpic.set_title(string("Component No. "+num2str(cnum)+
				    " - thresholded IC map ("+
				    num2str(tctr)+")  ")+
			     mmodel.get_infstr(tctr));
	    newpic.set_cbar(string("ysb"));
	    //cerr << instr << endl;
	    newpic.slicer(newvol, instr, &melodat.tempInfo); 

	    IC_rep_det(mmodel, cnum, dim);
	    IChtml << "<a href=\""+mmodel.get_prefix()+"_MM.html\">";
	    IChtml << "<img BORDER=0 SRC=\""+mmodel.get_prefix()
	      +"_thresh"+num2str(tctr)+".png\" ALT=\"MMfit\"></A><p>" << endl;
	  }
	}


      { //finish IC page
      IChtml<< "<HR><FONT SIZE=1>This page produced automatically by "
	    << "<A HREF=\"http://www.fmrib.ox.ac.uk/fsl/melodic/index.html\"> MELODIC </A>" 
	    << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	    << "FMRIB Software Library</A>.</FONT>" << endl
	    << "</BODY></HTML>" << endl;
      } //finish IC page
      IC_rep_det(mmodel, cnum, dim);
    }
  }

  void MelodicReport::IC_rep_det(MelGMix &mmodel, int cnum, int dim)
  {
    if( bool(opts.genreport.value()) ){

      {//start IC2 page
	IChtml2.setDir(report.getDir(),mmodel.get_prefix()+"_MM.html");
	IChtml2 << "<HTML> " << endl
		<< "<TITLE>Component " << num2str(cnum)
		<< " Mixture Model fit </TITLE>" << endl
		<< "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
		<< "/doc/images/fsl-bg.jpg\">" << endl 
		<< "<hr><CENTER><H1>Component " << num2str(cnum)
		<< " Mixture Model fit </H1>"<< endl;
     
	if(cnum>1)
	  IChtml2 << "<a href=\"" << string("IC_")+num2str(cnum-1)
		 <<"_MM.html\">previous</a>&nbsp;-&nbsp;";
	
	IChtml2 << "<a href=\""+ mmodel.get_prefix() + 
	  ".html\">&nbsp;up&nbsp;to IC report&nbsp;</a>";

	if(cnum<dim)
	  IChtml2 << "&nbsp;-&nbsp;<a href=\"" << string("IC_")+num2str(cnum+1)
		 <<"_MM.html\">next</a><p>";
	IChtml2 << "<hr><p>" << endl;
      }

      volume4D<float> tempVol; 

      if(melodat.get_IC().Storage()>0)
      {//Output raw IC map

	//	tempVol.setmatrix(melodat.get_IC().Row(cnum),
	//	  melodat.get_mask());
 	tempVol.setmatrix(mmodel.get_data(),
			  melodat.get_mask());
 	volume<float> map1;
	volume<float> map2;
	map1 = threshold(tempVol[0],float(0.0), 
			 tempVol[0].max());
	map2 = threshold(tempVol[0],tempVol[0].min(), 
			 float(-0.0));

	volume<float> newvol; 
	miscpic newpic;

	//	float map1min = (map1 + binarise(tempVol[0],tempVol[0].min(), 
	//		float(0.0)) * map1.max()).robustmin();
	float map1max = map1.percentile(0.99);
	float map2min = map2.percentile(0.01);
	//float map2max = (map2 + binarise(tempVol[0],float(0.0), 
	//		tempVol[0].max()) * map2.min()).robustmax();
	
 	newpic.overlay(newvol, melodat.get_mean(), map1, map2, 
 		       float(0.0),
 		       float(0.0),
 		       float(0.01), map1max, float(-0.01), map2min, 
 		       0, 0, &melodat.tempInfo);

 	char instr[10000];
	
	//save_volume(newvol,report.appendDir(mmodel.get_prefix()+"rendered"),
	//      melodat.tempInfo);
 	sprintf(instr," ");
 	strcat(instr," -s 2");
 	strcat(instr," -A 950 ");
 	strcat(instr,string(report.appendDir(mmodel.get_prefix()+
 					     ".png")).c_str());
 	newpic.set_title(string("Component No. "+num2str(cnum)+
 				" - raw Z transformed IC map (1 - 99 percentile)"));
 	newpic.set_cbar(string("ysb"));
	
 	newpic.slicer(newvol, instr, &melodat.tempInfo);
      }
      IChtml2 << "<a href=\""+mmodel.get_prefix()+".html\">";
      IChtml2 << "<img BORDER=0 SRC=\""+ mmodel.get_prefix()+
	".png\"><A><p>" << endl;
	
      

      if(mmodel.get_probmap().Storage()>0&&
	 (mmodel.get_probmap().Ncols() == mmodel.get_data().Ncols())&&
	 (mmodel.get_probmap().Nrows() == mmodel.get_data().Nrows()))
	{//Output probmap  
	  tempVol.setmatrix(mmodel.get_probmap(),melodat.get_mask());
	  
	  volume<float> map;
	  map = tempVol[0];
      
	  volume<float> newvol; 
	  miscpic newpic;

	  newpic.overlay(newvol, melodat.get_mean(), map, map, 
			 melodat.get_mean().percentile(0.02),
			 melodat.get_mean().percentile(0.98),
			 float(0.1), float(1.0), float(0.0), float(0.0),
			 0, 0, &melodat.tempInfo);

	  //newpic.set_minmax(float(0.0),float(0.0),float(0.0),
	  //	    float(1.0),float(0.0),float(0.0));

	  char instr[10000];

	  sprintf(instr," ");
	  strcat(instr,"-l render1 -s 2");
	  strcat(instr," -A 950 ");
	  strcat(instr,string(report.appendDir(mmodel.get_prefix()+
					       "_prob.png")).c_str());      
	  newpic.set_title(string("Component No. "+num2str(cnum)+
 			      " - Mixture Model probability map"));
      
	  newpic.set_cbar(string("y"));
	  newpic.slicer(newvol, instr, &melodat.tempInfo); 

// 	  {
// 	    Log MMdetails;
// 	    MMdetails.setDir(report.getDir(),mmodel.get_prefix()+"_MM.txt");
// 	    MMdetails << mmodel.get_prefix() << " Mixture Model fit " 
// 		      << endl << endl;
// 	    MMdetails << " Means :  " << mmodel.get_means() << endl
// 		      << " Vars  :  " << mmodel.get_vars()  << endl
// 		      << " Prop. :  " << mmodel.get_pi()    << endl;
// 	   // if(mmodel.get_type()=="GGM"){
// 	   //   MMdetails << endl << " Gamma offset :  " 
// 	//		<<  mmodel.get_offset() <<  endl;
// 	 //   }
// 	  }
	  //  IChtml2 << "<A HREF=\"" +mmodel.get_prefix()+"_MM.txt\"> ";
	  IChtml2 << "<a href=\""+mmodel.get_prefix()+".html\">";
	  IChtml2 << "<img BORDER=0 SRC=\""+ mmodel.get_prefix()+
	    "_prob.png\">" << endl;
	  //	  IChtml2 << "</A>";
	  IChtml2 << "</A><p>" << endl;
	}

      {//Output GGM/GMM fit
	miscplot newplot;

// 	cerr << endl << endl;
//  	cerr << mmodel.get_means() << endl;
//  	cerr << mmodel.get_vars()  << endl;
//  	cerr << mmodel.get_pi()    << endl;

	if(mmodel.get_type()=="GGM"){
	  newplot.add_label("IC map histogram");
	  newplot.add_label("full GGM fit");
	  newplot.add_label("background Gaussian");
	  newplot.add_label("Gamma distributions");
	  newplot.gmmfit(mmodel.get_data().Row(1),
			 mmodel.get_means(),
			 mmodel.get_vars(),
			 mmodel.get_pi(),
			 report.appendDir(mmodel.get_prefix()+"_MMfit.png"),
			 string(mmodel.get_prefix() +
				" GGM("+num2str(mmodel.mixtures())+") fit"),
			 true, float(0.0), float(2.0)); 

	  //  newplot.histogram(mmodel.get_data().Row(1),
	  //  report.appendDir(mmodel.get_prefix()+"_MMfit.png"),
	  // 	     string(mmodel.get_prefix() +
	  // 	     " GGM("+num2str(mmodel.mixtures())+") fit"));

	  //, mmodel.get_offset());   
	}
	else{
	  newplot.add_label("IC map histogram");
	  newplot.add_label("full GMM fit");
	  newplot.add_label("individual Gaussians");
	  newplot.gmmfit(mmodel.get_data().Row(1),
			 mmodel.get_means(),
			 mmodel.get_vars(),
			 mmodel.get_pi(),
			 report.appendDir(mmodel.get_prefix()+"_MMfit.png"),
			 string(mmodel.get_prefix() +
				" GMM("+num2str(mmodel.mixtures())+") fit"), 
			 false, float(0.0), float(2.0));
   
	}
	//	cerr << "finish HTML2 " << endl;
	IChtml2 << "<A HREF=\"" +mmodel.get_prefix()+"_MMfit_detail.png\"> ";
	IChtml2 << "<img BORDER=0 SRC=\""+ mmodel.get_prefix()+
	  "_MMfit.png\"></A><p>" << endl;
      } //GGM/GMM plot
      
	  {
	    IChtml2 << "<br> &nbsp;" << mmodel.get_prefix() 
		    << " Mixture Model fit <br>" << endl
		    << "<br> &nbsp; Means :  " << mmodel.get_means() << endl
		    << "<br> &nbsp;  Vars  :  " << mmodel.get_vars()  << endl
		    << "<br> &nbsp;  Prop. :  " << mmodel.get_pi()    << endl;
	    // if(mmodel.get_type()=="GGM"){
	    // 	     IChtml2 << "<br> &nbsp;   Gamma offset :  " 
	    // 		     << mmodel.get_offset() << "<br><p>"<<  endl;
	    // 	    }
	  }

      { //finish IC2 page
	IChtml2<< "<HR><FONT SIZE=1>This page produced automatically by "
	       << "<A HREF=\"http://www.fmrib.ox.ac.uk/fsl/melodic/index.html\"> MELODIC </A>" 
	       << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	       << "FMRIB Software Library</A>.</FONT>" << endl
	       << "</BODY></HTML>" << endl;
      } //finish IC2 page
    }
  }


  void MelodicReport::IC_simplerep(string prefix, int cnum, int dim)
  {
    if( bool(opts.genreport.value()) ){
      addlink(prefix+".html",num2str(cnum));
      IChtml.setDir(report.getDir(),prefix+".html");
      
      {//start IC page
	
	IChtml << "<HTML> " << endl
	       << "<TITLE>MELODIC Component " << num2str(cnum)
	       << "</TITLE>" << endl
	       << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	       << "/doc/images/fsl-bg.jpg\">" << endl 
	       << "<hr><CENTER><H1>MELODIC Component " << num2str(cnum)
	       << "</H1>"<< endl;
	
	if(cnum>1)
	  IChtml << "<a href=\"" << string("IC_")+num2str(cnum-1)
		 <<".html\">previous</a>&nbsp;-&nbsp;";
	
	IChtml << "<a href=\"00index.html\">&nbsp;index&nbsp;</a>" ;
	
	if(cnum<dim)
	  IChtml << "&nbsp;-&nbsp;<a href=\"" << string("IC_")+num2str(cnum+1)
		 <<".html\">next</a><p>";
	
	IChtml << "<hr><p>" << endl;
      }


      volume4D<float> tempVol; 

      if(melodat.get_IC().Storage()>0)
      {//Output raw IC map

	tempVol.setmatrix(melodat.get_IC().Row(cnum),
			  melodat.get_mask());
 	volume<float> map1;
	volume<float> map2;
	map1 = threshold(tempVol[0],float(0.0), 
			 tempVol[0].max());
	map2 = threshold(tempVol[0],tempVol[0].min(), 
			 float(-0.0));

	volume<float> newvol; 
	miscpic newpic;

	//	float map1min = (map1 + binarise(tempVol[0],tempVol[0].min(), 
	//		float(0.0)) * map1.max()).robustmin();
	float map1max = map1.percentile(0.99);
	float map2min = map2.percentile(0.01);
	//float map2max = (map2 + binarise(tempVol[0],float(0.0), 
	//		tempVol[0].max()) * map2.min()).robustmax();
	
 	newpic.overlay(newvol, melodat.get_mean(), map1, map2, 
 		       float(0.0),
 		       float(0.0),
 		       float(0.01), map1max, float(-0.01), map2min, 
 		       0, 0, &melodat.tempInfo);

 	char instr[10000];
	
	//save_volume(newvol,report.appendDir(prefix+"rendered"),
	//      melodat.tempInfo);
 	sprintf(instr," ");
 	strcat(instr," -s 2");
 	strcat(instr," -A 950 ");
 	strcat(instr,string(report.appendDir(prefix+
 					     ".png")).c_str());
 	newpic.set_title(string("Component No. "+num2str(cnum)+
 				" - raw Z transformed IC map (1 - 99 percentile)"));
 	newpic.set_cbar(string("ysb"));
	
 	newpic.slicer(newvol, instr, &melodat.tempInfo);
      }
     
      IChtml << "<img BORDER=0 SRC=\""+ prefix+
	".png\"><p>" << endl;
	
      {//plot time course
	miscplot newplot;
	
	if(opts.tr.value()>0.0)
	  newplot.timeseries(melodat.get_mix().Column(cnum).t(),
			     report.appendDir(string("t")+
					      num2str(cnum)+".png"),
			     string("Timecourse (in seconds); TR = ")+
			     float2str(opts.tr.value(),0,2,0)+" s", 
			     opts.tr.value(),150,4,1);
	else
	  newplot.timeseries(melodat.get_mix().Column(cnum).t(),
			     report.appendDir(string("t")+
					      num2str(cnum)+".png"),
			     string("Timecourse (in TRs)"));
	write_ascii_matrix(report.appendDir(string("t")
					    +num2str(cnum)+".txt"),
			   melodat.get_mix().Column(cnum));
	IChtml << "<A HREF=\"" << string("t")
	  +num2str(cnum)+".txt" << "\"> ";
	IChtml << "<img BORDER=0 SRC=\"" 
	  +string("t")+num2str(cnum)+".png\"></A><p>" << endl;
      }//time series plot
      
      {//plot frequency 
	miscplot newplot;
	int fact = int(std::pow(10.0,
		       int(std::log10(float(melodat.get_mix().Nrows())))));
	
	if(opts.tr.value()>0.0)
	  newplot.timeseries(melodat.get_fmix().Column(cnum).t(),
			     report.appendDir(string("f")+
					      num2str(cnum)+".png"),
			     string("FFT of timecourse (in Hz / ") +
			     num2str(fact)+")",
			     fact/(opts.tr.value()*melodat.get_mix().Nrows()),
			     150,
			     0,2);
	else
	  newplot.timeseries(melodat.get_fmix().Column(cnum).t(),
			     report.appendDir(string("f")+
					      num2str(cnum)+".png"),
			     string(string("FFT of timecourse (in cycles); ")
				    +"frequency(Hz)=cycles/("
				    +num2str(melodat.get_mix().Nrows())
				    +"* TR); period(s)=("
				    +num2str(melodat.get_mix().Nrows())
				    +"* TR)/cycles"));
	write_ascii_matrix(report.appendDir(string("f")
					    +num2str(cnum)+".txt"),
			   melodat.get_mix().Column(cnum));
	IChtml << "<A HREF=\"" << string("f")
	  +num2str(cnum)+".txt" << "\"> ";
	IChtml << "<img BORDER=0 SRC=\"" 
	  +string("f")+num2str(cnum)+".png\"></A><p>" << endl;
      }//frequency plot
 
      { //finish IC page
	IChtml<< "<HR><FONT SIZE=1>This page produced automatically by "
	      << "<A HREF=\"http://www.fmrib.ox.ac.uk/fsl/melodic/index.html\"> MELODIC </A>" 
	      << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	      << "FMRIB Software Library</A>.</FONT>" << endl
	      << "</BODY></HTML>" << endl;
      } //finish IC page 
    } 
  }

  void MelodicReport::PPCA_rep(){
    
    {//plot time course
      report << "<p> <H3>PCA estimates </H3> <p>" << endl;
 
      Matrix what;
      miscplot newplot;

      what  = melodat.get_EV();
      what &= melodat.get_EVP();

      newplot.add_label("ordered Eigenvalues");
      newplot.add_label("% of expl. variance");

      if(melodat.get_PPCA().Storage()>0){
	what = what.Columns(1,melodat.get_PPCA().Nrows());
	what &= melodat.get_PPCA().t();
	newplot.add_label("dim. estimate");
      }

      newplot.timeseries(what,
			 report.appendDir("EVplot.png"),
			 string("Eigenspectrum Analysis"), 
			 0,450,4,0);

      report << "<img BORDER=0 SRC=\"EVplot.png\"><p>" << endl;
    }//time series plot
  }
}
