/*  miscplot.cc  -- miscellaneous plotting functions

    Christian F. Beckmann, FMRIB Image Analysis Group
    
    Copyright (C) 1999-2002 University of Oxford */

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

#include <sstream>
#include <stdio.h>

#include "miscplot.h"
#include "miscmaths/miscprob.h"
#include "gdfonts.h"

extern "C" {
 #include "gd.h"
 #include "gdc.h"
 #include "gdchart.h"
}

namespace MISCPLOT{

//  MATLAB Line colors
unsigned long   sc[64] = {0x8080FF,0x008000,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x404040,0x0000FF,
			  0x008000,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x404040,0x0000FF,0x008000,
			  0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x404040,0x0000FF,0x008000,0xFF0000,
			  0x00BFBF,0xBF00BF,0xBFBF00,0x404040,0x0000FF,0x008000,0xFF0000,0x00BFBF,
			  0xBF00BF,0xBFBF00,0x404040,0x0000FF,0x008000,0xFF0000,0x00BFBF,0xBF00BF,
			  0xBFBF00,0x404040,0x0000FF,0x008000,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,
			  0x404040,0x0000FF,0x008000,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x404040,
			  0x0000FF,0x008000,0xFF0000,0x00BFBF,0xBF00BF,0xBFBF00,0x404040,0x0000FF};

using namespace NEWMAT;
using namespace MISCMATHS;


string int2str(int n)
  {
    ostringstream os;
    os.setf(ios::internal, ios::adjustfield);
    os << n;
    return os.str();
  }

string float2str(float f, int width, int prec, bool scientif)
  {
    ostringstream os;
    int redw = int(std::abs(std::log10(std::abs(f))))+1;
    if(width>0)
      os.width(width);
    if(scientif)
      os.setf(ios::scientific);
    os.precision(redw+std::abs(prec));
    os.setf(ios::internal, ios::adjustfield);
    os << f;
    return os.str();
  }

void miscplot::add_legend(void* ptr, unsigned long cmap[], bool inside){
  
  int xsize = gdImagePtr(ptr)->sx;
  int ysize = gdImagePtr(ptr)->sy;
  int dest_x = 0;

  int linebrk = 5;
  int space = gdFontSmall->w + 2;
  int linelength = 15;

  int rlbllength = 0;
  for(int ctr= 0; ctr < int(labels.size()); ctr++)
    rlbllength = max(rlbllength,int(labels[ctr].length()));

  if(explabel.length()>0)
    ysize += space + 0*gdFontSmall->h;

  if(xlabels.size()>0)
    ysize += xlabels.size()*(linebrk + gdFontSmall->h)+linebrk;

  if(ylabels.size()>0)
    dest_x += ylabels.size()*(linebrk + gdFontSmall->h)+2*linebrk;

  if((!inside)&&(labels.size()>0))
    xsize += 2*space + linelength + rlbllength * gdFontSmall->w;
  
  xsize += dest_x;
  gdImagePtr newim;
  newim = gdImageCreate(xsize,ysize);
  gdImageCopy(newim,gdImagePtr(ptr),dest_x,0,0,0,gdImagePtr(ptr)->sx,gdImagePtr(ptr)->sy);

  int black = gdImageColorResolve(newim, 0, 0, 0);
  int ycoor, xcoor, yoffset=0;

  //make scale info
  if(explabel.length()>0){
    yoffset=space;

    xcoor=gdImagePtr(ptr)->sx - (4 + explabel.length())*gdFontSmall->w;
    ycoor=gdImagePtr(ptr)->sy+linebrk;
    char *s = (char*)string("x10").c_str();
    gdImageString(newim, gdFontSmall,xcoor,ycoor, (unsigned char*) s, black);
    xcoor+=(int)3.5*gdFontSmall->w;
    ycoor-=gdFontSmall->h/2;
    string s2 = string("-")+explabel;
    char *s3 = (char*)s2.c_str();
    gdImageString(newim, gdFontSmall,xcoor,ycoor, (unsigned char*) s3, black);
  }

  // make xlabel
  ycoor=gdImagePtr(ptr)->sy+linebrk+yoffset;
  for(int ctr=0; ctr < int(xlabels.size()); ctr++){
     xcoor = dest_x + gdImagePtr(ptr)->sx/2 - ((int)xlabels[ctr].length())/2*gdFontSmall->w;
     char *s = (char*)xlabels[ctr].c_str();
     gdImageString(newim, gdFontSmall,xcoor,ycoor, (unsigned char*) s, black);
     ycoor+=gdFontSmall->h + linebrk;
  }

  // make legend
  ycoor = 2+2*gdFontSmall->h;
  for(int ctr=0; ctr < int(labels.size()); ctr++){
    xcoor = dest_x+gdImagePtr(ptr)->sx;
    if(inside)
      xcoor -= 2*space + linelength + rlbllength * gdFontSmall->w;

    int col =  gdImageColorResolve(newim, ((cmap[ctr]&PVRED)>>16),((cmap[ctr]&PVGRN)>>8), (cmap[ctr]&0x000000FF));

    gdImageLine(newim, xcoor, ycoor + gdFontSmall->h/2, xcoor + linelength, ycoor + gdFontSmall->h/2,col);
    gdImageLine(newim, xcoor, ycoor + gdFontSmall->h/2+1, xcoor + linelength, ycoor + gdFontSmall->h/2+1,col);
    char *s = (char*)labels[ctr].c_str();
    gdImageString(newim, gdFontSmall, xcoor + linelength + space, ycoor, (unsigned char*) s, black);
    ycoor += gdFontSmall->h + linebrk;
  }

  // make ylabel
  xcoor=space;
  for(int ctr=0; ctr < int(ylabels.size()); ctr++){
    int ycoor = 3*gdImagePtr(ptr)->sy/5 + ((int)ylabels[ctr].length())/2*gdFontSmall->w;
     char *s = (char*)ylabels[ctr].c_str();
     gdImageStringUp(newim, gdFontSmall,xcoor,ycoor, (unsigned char*) s, black);
     xcoor+=gdFontSmall->h + linebrk;
  }

  outim = newim;
}

void miscplot::timeseries(const Matrix& mat, string filename, string title, 
			  float tr, int ysize, int width, int prec, bool sci){
  int numlines=mat.Nrows();
  int numpoint=mat.Ncols();

  float* data = new float[numlines*numpoint];
 
  for(int ctr1=1;ctr1<=numlines; ctr1++){
    for(int ctr2=1;ctr2<=numpoint; ctr2++){
      data[(ctr1-1)*numpoint + ctr2-1] = mat(ctr1,ctr2);
    }
  }

  int spacing = 1;

  if(numpoint>20) spacing = 2;
  if(numpoint>40) spacing = 5;
  if(numpoint>100) spacing = 10;
  if(numpoint>500) spacing = 20;
  if(numpoint>1200) spacing = 50;
  //if(numpoint>2000) spacing = 100;

  char* ctl = new char[numpoint];

  if((tr == 0.0)&&(spacing%2 == 1) || (tr > 0.0))
    ctl[0] = TRUE;
  else
    ctl[0] = FALSE;

  for(int ctr=1;ctr<numpoint;ctr++)
    ctl[ctr]=((((ctr+1)%spacing == 0)&&(tr == 0.0))
	      ||(((ctr)%spacing == 0)&&(tr != 0.0)));
  //ctl[numpoint-1]=TRUE;
 
  char** lbls = new char*[numpoint];
 
  string s;
  lbls[0] = new char[2];
  if(tr == 0.0)
    strcpy(lbls[0],string("1").c_str());
  else
    strcpy(lbls[0],string("0").c_str());
  
  for(int ctr=1;ctr<numpoint;ctr++){
    if(ctl[ctr]){
      if(tr == 0.0)
	s = int2str(ctr+1);
      else
	s = float2str((ctr)*tr,width,prec,sci);
      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
    else{
      lbls[ctr] = new char[1];
      strcpy(lbls[ctr],string("").c_str());
    }
  }


  {
    int tmp = _gdccfoo1;
    tmp++;
    long unsigned int tmp2 = _gdccfoo2;
    tmp2++;
  }

  GDC_xlabel_ctl = ctl;

  GDC_BGColor   = 0xFFFFFFL;                  /* backgound color (white) */
  GDC_LineColor = 0x000000L;                  /* line color      (black) */
  GDC_SetColor  = &(sc[0]);                   /* MATLAB-like colors */

  GDC_title = (char*)title.c_str();

  GDC_ticks = GDC_TICK_LABELS;
  GDC_grid  = GDC_TICK_NONE;
  GDC_yaxis = TRUE;
  GDC_yaxis2 = FALSE;
  GDC_xaxis = TRUE;
  GDC_xaxis_angle = 0;
  GDC_ylabel_density = 50;
  GDC_ylabel_fmt = "%.2f";
  GDC_title_size = GDC_SMALL;

  int i,j;
  GDC_requested_ymax = 1.10*mat.Maximum2(i,j);
  GDC_requested_ymin = 1.10*mat.Minimum2(i,j);

  int xsize = std::max(4*numpoint,750);
  if(xsize>1200) xsize=(int)floor(std::max(xsize/2.0,750.0));

  if(req_xsize>0){
    xsize=req_xsize;
    ysize=req_ysize;
  }

  //GDC_hard_size = TRUE;
  //GDC_hard_xorig = 50;
  //GDC_hard_graphwidth = xsize - 65;

  if(filename.substr(filename.size()-4,filename.size())!=string(".png"))
    filename += string(".png");
  
  FILE  *outpng1 = fopen(filename.c_str(), "wb" );
  GDC_image_type     = GDC_PNG;

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0)
    GDC_hold_img = GDC_EXPOSE_IMAGE;
  
  GDC_out_graph( xsize, ysize, outpng1 , GDC_LINE, numpoint, (char**) lbls , numlines, data, NULL ); 
  //GDC_out_graph( xsize, ysize, outpng1 , GDC_LINE, numpoint, NULL , numlines, data, NULL );
  fclose( outpng1 );
  
  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0){
    outpng1 = fopen(filename.c_str(), "wb" );
    add_legend(GDC_image, sc);
    gdImagePng(outim, outpng1); 
    fclose( outpng1 );
    GDC_destroy_image(GDC_image);
    if(outim) gdImageDestroy(outim);
  }
 
  for(int ctr=0;ctr<numpoint;ctr++){
    delete [] lbls[ctr];
  }
    
  delete [] data;
  delete [] ctl;
  delete [] lbls;
}



void miscplot::histogram(const Matrix& mat, string filename, string title){

  RowVector datam = mat.Row(1);
  int numpoint=datam.Ncols();
  int i,j;
  double scale=std::max(1.0,MISCMATHS::pow((float)10.0,(double)-(std::floor(std::log10(std::min(std::abs(datam.Maximum2(i,j)),std::abs(datam.Minimum2(i,j)))/2.0)))));

  //cerr << datam.Maximum2(i,j) << "  " << datam.Minimum2(i,j) << scale << endl;
  datam *=scale;
  if(scale>100)
    explabel=num2str(std::log10(scale));

  float tmax = datam.Maximum2(i,j);
  float tmin = datam.Minimum2(i,j); 
  float trange = tmax-tmin;
  int bins = (int)floor(MISCMATHS::sqrt(numpoint));
  float intsize = trange / bins;
  int xlint = (int)ceil(trange/6);
  int binperint = (int)ceil(xlint / intsize);

  //cerr << tmax << " " << tmin  << " " << trange  << " " << bins << " " << intsize  << " " << xlint  << " " << binperint << " :" << scale <<endl; 
  intsize = float(xlint) / float(binperint);

  float xmin = ceil(std::abs(1.02*tmin)/intsize)*intsize * sign(tmin);
  float xmax = ceil(std::abs(1.02*tmax)/intsize)*intsize * sign(tmax);

  bins = (int)((xmax-xmin)/intsize);

  Matrix bindata(1,bins);
  bindata = 0.0;

  double binsize = (xmax-xmin)/std::max(bins,1);
  for(int ctr = 1; ctr<= datam.Ncols(); ctr++){
    bindata(1,std::max(std::min(int(floor((datam(ctr) - xmin) / binsize) + 1),bins),1))++;
  }
  //cerr << " calculated bindata " << endl;

  int factor = 1;
  numpoint = factor * bindata.Ncols();
  float* data = new float[numpoint]; 
  for(int ctr1=0;ctr1< numpoint; ctr1++){
      data[ctr1] = bindata(1,(int)floor(float(ctr1 / factor + 1)));
  }
  
  RowVector xax(numpoint);
  for(int ctr1=0;ctr1< numpoint; ctr1++){
    xax(ctr1+1) = xmin + (ctr1+ 0.5) * intsize; 
  }

  char* ctl = new char[numpoint];
  xlint = (int)ceil((xmax-xmin)/9);
  
  ctl[0]=FALSE;
  for(int ctr=1;ctr<numpoint;ctr++){
    int lblpoint = (int)(MISCMATHS::round(abs(xax(ctr))/xlint)*xlint);
    if(xax(ctr)<0) lblpoint *= -1;
     if( (xax(ctr)<lblpoint)&&(xax(ctr+1)>lblpoint) ){
       ctl[ctr]=TRUE;
     }
     else
       ctl[ctr]=FALSE;
  } 

  char** lbls = new char*[numpoint];
 
  string s;

  for(int ctr=0;ctr<numpoint;ctr++){
    if(ctl[ctr]){
      if(scale<2||scale>100)
	s = float2str((int)MISCMATHS::round(xax(ctr+1)/xlint)*xlint,3,2,false);
      else
	s = float2str((float)MISCMATHS::round(xax(ctr+1)/xlint)*xlint/scale,3,2,false);
      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
    else{
      lbls[ctr] = new char[1];
      strcpy(lbls[ctr],string("").c_str());
    }
  }
  
 
  GDC_xlabel_ctl = ctl;

  GDC_BGColor   = 0xFFFFFFL;                  /* backgound color (white) */
  GDC_LineColor = 0x000000L;                  /* line color      (black) */
  GDC_SetColor  = &(sc[0]);                   /* MATLAB-like colors */

  GDC_title = (char*)title.c_str();
  GDC_title_size = GDC_SMALL;

  GDC_ticks = GDC_TICK_LABELS;
  GDC_grid  = GDC_TICK_NONE;
  GDC_yaxis = ylabels.size()>0;
  GDC_yaxis2 = FALSE;
  GDC_xaxis = TRUE;
  GDC_xaxis_angle = 0;
  GDC_bar_width =  75;
  GDC_ylabel_density = 50;
  GDC_requested_ymax = 1.15*bindata.Maximum2(i,j);
  GDC_requested_ymin = 1.15*bindata.Minimum2(i,j);

  //int xsize = std::max(4*numpoint,750);
  //  if(xsize>1200) xsize=(int)floor(std::max(xsize/2.0,750.0));
  int xsize = 600;
  int ysize = 400;

  if(req_xsize>0){
    xsize=req_xsize;
    ysize=req_ysize;
  }


  if(filename.substr(filename.size()-4,filename.size())!=string(".png"))
    filename += string(".png");

  FILE  *outpng1 = fopen(filename.c_str(), "wb" );
  GDC_image_type     = GDC_PNG;

 if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0)
    GDC_hold_img = GDC_EXPOSE_IMAGE;
  
  GDC_out_graph( xsize, ysize, outpng1 , GDC_BAR, numpoint, (char**) lbls , 1, data, NULL ); 
  fclose( outpng1 );

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0){
    outpng1 = fopen(filename.c_str(), "wb" );
    add_legend(GDC_image, sc);
    
    gdImagePng(outim, outpng1); 
    fclose( outpng1 );
    GDC_destroy_image(GDC_image);
    if(outim) gdImageDestroy(outim);
  }

  for(int ctr=0;ctr<numpoint;ctr++){
    delete [] lbls[ctr];
  }
    
  delete [] data;
  delete [] lbls;
  delete [] ctl;
}

void miscplot::gmmfit(const Matrix& mat, Matrix& mu, Matrix& sig, 
		      Matrix& pi, string filename, string title, 
		      bool gammamix, float meanoffset, float detailfactor)
{
  RowVector datam = mat.Row(1);
  int numpoint=datam.Ncols();
  int i,j;
  double scale=std::max(1.0,MISCMATHS::pow((float)10.0,(double)-(std::floor(std::log10(std::min(std::abs(datam.Maximum2(i,j)),std::abs(datam.Minimum2(i,j)))/2.0)))));

  //cerr << datam.Maximum2(i,j) << "  " << datam.Minimum2(i,j) << scale << endl;
  datam *=scale; mu *= scale; sig *= scale;
  if(scale>100)
    explabel=num2str(std::log10(scale));

  float tmax = datam.Maximum2(i,j);
  float tmin = datam.Minimum2(i,j); 
  float trange = tmax-tmin;
  int bins = (int)floor(MISCMATHS::sqrt(numpoint));
  float intsize = trange / bins;
  int xlint = (int)ceil(trange/6);
  int binperint = (int)ceil(xlint / intsize);

  //cerr << tmax << " " << tmin  << " " << trange  << " " << bins << " " << intsize  << " " << xlint  << " " << binperint << " :" << scale <<endl; 
  intsize = float(xlint) / float(binperint);

  float xmin = ceil(std::abs(1.02*tmin)/intsize)*intsize * sign(tmin);
  float xmax = ceil(std::abs(1.02*tmax)/intsize)*intsize * sign(tmax);

  bins = (int)((xmax-xmin)/intsize);

  Matrix bindata(1,bins);
  bindata = 0.0;

  double binsize = (xmax-xmin)/std::max(bins,1);
  for(int ctr = 1; ctr<= datam.Ncols(); ctr++){
    bindata(1,std::max(std::min(int(floor((datam(ctr) - xmin) / binsize) + 1),bins),1))++;
  }
  //cerr << " calculated bindata " << endl;

  int factor = 5;

  bindata = bindata / max(double(factor*bindata.SumAbsoluteValue()),double(1.0));

  numpoint = factor * bindata.Ncols();
  float* histdata = new float[numpoint]; 
  for(int ctr1=0;ctr1< numpoint; ctr1++){
      histdata[ctr1] = bindata(1,(int)floor(float(ctr1 / factor + 1)));
  }
  
  RowVector xax(numpoint);
  for(int ctr1=0;ctr1< numpoint; ctr1++){
    xax(ctr1+1) = xmin + (ctr1+ 0.5) * intsize/factor; 
  }

  char* ctl = new char[numpoint];
  xlint = (int)ceil((xmax-xmin)/(factor*5));
  
  ctl[0]=FALSE;
  for(int ctr=1;ctr<numpoint;ctr++){
    int lblpoint = (int)(MISCMATHS::round(abs(xax(ctr))/xlint)*xlint);
    if(xax(ctr)<0) lblpoint *= -1;
     if( (xax(ctr)<lblpoint)&&(xax(ctr+1)>lblpoint) ){
       ctl[ctr]=TRUE;
     }
     else
       ctl[ctr]=FALSE;
  } 

  char** lbls = new char*[numpoint];
 
  string s;

  for(int ctr=0;ctr<numpoint;ctr++){
    if(ctl[ctr]){
      if(scale<2||scale>100)
	s = float2str((int)MISCMATHS::round(xax(ctr+1)/xlint)*xlint,3,2,false);
      else
	s = float2str((float)MISCMATHS::round(xax(ctr+1)/xlint)*xlint/scale,3,2,false);
      lbls[ctr] = new char[s.length()+1];
      strcpy(lbls[ctr],s.c_str());
    }
    else{
      lbls[ctr] = new char[1];
      strcpy(lbls[ctr],string("").c_str());
    }
  }
  

  // Calculate lines
  int numlines=mu.Ncols()+1;
  for(int ctr1=1;ctr1< sig.Ncols(); ctr1++)
    if(sig(1,ctr1)<0.000001){
      sig(1,ctr1) = 0.000001;
      pi(1,ctr1) = 0.0;
    }

  Matrix fit;
  fit = SP(normpdf(xax,mu,sig),pi.t()*ones(1,numpoint));

  if(null_shift!=0.0)
    fit.Row(1) = pi(1,1)*gammapdf(-(xax + null_shift) + meanoffset, std::abs(mu(1,1) + null_shift), sig(1,1));

  if(gammamix){
    if(pi(1,2)>0.0001)
      fit.Row(2)  = pi(1,2)*gammapdf( xax - meanoffset, std::abs(mu(1,2)), sig(1,2));
    if(pi(1,3)>0.0001)
      fit.Row(3)  = pi(1,3)*gammapdf(-xax + meanoffset, std::abs(mu(1,3)), sig(1,3));
  }
  fit = sum(fit,1) & fit;
  fit = fit / fit.Row(1).SumAbsoluteValue();

  float* linesdata = new float[(numlines)*numpoint];
  for(int ctr1=1;ctr1<=numlines; ctr1++){
    for(int ctr2=1;ctr2<=numpoint; ctr2++){
      linesdata[(ctr1-1)*numpoint + ctr2-1] = fit(ctr1,ctr2);
    }
  }	

  GDC_xlabel_ctl = ctl;

  unsigned long*   sc2 = new unsigned long[numlines];
  for(int ctr=0;ctr<numlines;ctr++) sc2[ctr] = 0xFFDD00;
  sc2[0]=0xFF0000;
  if(gammamix)
    sc2[1]=0x00aa00;

  unsigned long  sc3[64];
  for(int ctr=0;ctr<64;ctr++) sc3[ctr] = sc[ctr];
  for(int ctr=0;ctr<numlines;ctr++){
    sc3[ctr+1] = sc2[ctr];
  }

  GDC_BGColor   = 0xFFFFFFL;                  /* backgound color (white) */
  GDC_LineColor = 0x000000L;                  /* line color      (black) */
  GDC_SetColor  = &(sc2[0]);                  /* MATLAB-like colors */

  GDC_title = (char*)title.c_str();
  GDC_title_size = GDC_SMALL;

  GDC_ticks = GDC_TICK_LABELS;
  GDC_grid  = GDC_TICK_NONE;
  GDC_yaxis = ylabels.size()>0;
  GDC_yaxis2 = FALSE;
  GDC_xaxis = TRUE;
  GDC_xaxis_angle = 0;
  GDC_bar_width =  75;
  GDC_ylabel_density = 50;
  GDC_requested_ymax = 1.15*bindata.Maximum2(i,j);
  GDC_requested_ymin = 1.15*bindata.Minimum2(i,j);

  //int xsize = std::max(4*numpoint,750);
  //  if(xsize>1200) xsize=(int)floor(std::max(xsize/2.0,750.0));
  int xsize = 600;
  int ysize = 400;

  if(req_xsize>0){
    xsize=req_xsize;
    ysize=req_ysize;
  }


  if(filename.substr(filename.size()-4,filename.size())!=string(".png"))
    filename += string(".png");

  FILE  *outpng1 = fopen(filename.c_str(), "wb" );
  GDC_image_type     = GDC_PNG;

 if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0)
    GDC_hold_img = GDC_EXPOSE_IMAGE;
  
  GDC_out_graph( xsize, ysize, outpng1 , GDC_COMBO_LINE_AREA, numpoint, (char**) lbls , numlines, linesdata, histdata ); 
  fclose( outpng1 );

  if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0){
    outpng1 = fopen(filename.c_str(), "wb" );
    add_legend(GDC_image, sc3, TRUE);
    
    gdImagePng(outim, outpng1); 
    fclose( outpng1 );
    GDC_destroy_image(GDC_image);
    if(outim) gdImageDestroy(outim);
  }

 if((detailfactor > 0.0)&&(fit.Nrows()>2)){
  
    Matrix tmp(1,fit.Ncols()-2);
    tmp=sum(fit.Rows(3,fit.Nrows()),1);
   
    float req_max = detailfactor*tmp.Maximum2(i,j);

    if(req_max<0.0000001)
      req_max= 1.05*std::max(bindata.Maximum2(i,j),fit.Maximum2(i,j));
    
    for(int ctr1=1;ctr1<=numlines; ctr1++){
      for(int ctr2=1;ctr2<=numpoint; ctr2++){
	linesdata[(ctr1-1)*numpoint + ctr2-1] = std::min(req_max,linesdata[(ctr1-1)*numpoint + ctr2-1]);
      }
    }	

    for(int ctr1=0;ctr1< numpoint; ctr1++){
     histdata[ctr1] =  std::min(req_max, histdata[ctr1]);
    }

    GDC_requested_ymax = req_max;
    GDC_requested_ymin = float(0.0);
    GDC_title = (char*)title.append("  (detail)").c_str();
    GDC_title_size = GDC_SMALL;

    if(filename.substr(filename.size()-4,filename.size())==string(".png"))
      filename = filename.substr(0,filename.size()-4);

    FILE  *outpng1 = fopen((filename.append("_detail.png")).c_str(), "wb" );
    GDC_image_type = GDC_PNG;
    
    if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0)
      GDC_hold_img = GDC_EXPOSE_IMAGE;
    
    GDC_out_graph( xsize, ysize, outpng1 , GDC_COMBO_LINE_BAR, numpoint, (char**) lbls , numlines, linesdata, histdata); 
    //GDC_out_graph( xsize, ysize, outpng1 , GDC_COMBO_LINE_BAR, numpoint, NULL , numlines, data, histdata); 
    fclose( outpng1 );
    
    if (labels.size()>0||ylabels.size()>0||xlabels.size()>0||explabel.length()>0){
      outpng1 = fopen(filename.c_str(), "wb" );
      add_legend(GDC_image, sc3,TRUE);
      
      gdImagePng(outim, outpng1); 
      fclose( outpng1 );
      GDC_destroy_image(GDC_image);
      if(outim) gdImageDestroy(outim);
    }
 }

  for(int ctr=0;ctr<numpoint;ctr++){
    delete [] lbls[ctr];
  }
    
  delete [] histdata;
  delete [] linesdata;
  delete [] lbls;
  delete [] ctl;
}


}
