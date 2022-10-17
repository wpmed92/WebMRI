/* {{{ Copyright etc. */

/*  tsplot - FMRI time series and model plotting

    Stephen Smith and Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2006 University of Oxford  */

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

/* }}} */
/* {{{ background theory */

/*

GLM : Y = Xb + e

The "partial model fit" shows, in the case of a contrast which selects
a single EV, the part of the full model fit explained by that EV. In
the case of a more complex contrast, it is basically a plot of the sum
of each EV weighted by the product of that EV's PE and that EV's
contrast weighting.

i.e. the partial model fit is X*diag(c)*b, where X is the design and c
is the contrast vector (renormalised to unit length and then turned
into a diagonal matrix) and b is the parameter estimate vector.

Thus we plot this versus the "reduced data" ie plot

Y - Xb + X*diag(c)*b  vs  X*diag(c)*b
i.e.
residuals + partial fit    vs   partial fit

NOTE: this plot cannot simply be used to generate the t/Z score
associated with this contrast (eg by straight correlation) - this
would not be correct. In order to do that you would need to correlate
the Hansen projection c'*pinv(X) with Y instead.

*/

/* }}} */
/* {{{ defines, includes and typedefs */

#include "featlib.h"

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("Usage: tsplot <feat_directory.feat> [options]\n");
  printf("[-f <4D_data>] input main filtered data, in case it's not <feat_directory.feat>/filtered_func_data\n");
  printf("[-c <X Y Z>] : use X,Y,Z instead of max Z stat position\n");
  printf("[-C <X Y Z output_file.txt>] : use X,Y,Z to output time series only - no stats or modelling\n");
  printf("[-m <mask>] : use mask image instead of thresholded activation images\n");
  printf("[-o <output_directory>] change output directory from default of input feat directory\n");
  printf("[-n] don't weight cluster averaging with Z stats\n");
  printf("[-p] prewhiten data and model timeseries before plotting\n");
  exit(1);
}

/* }}} */

int main(int argc, char **argv)
{
  /* {{{ variables */

FILE         *ifp, *ofp, *rofp;
int          argi=1, prewhiten=0, ev, i, j, t, v, X=0, Y=0, Z=0, nevs, npts,
  ncon=1, nftests=0, size, ymin, ymax, coordset=0, dataonly=0, modelfree=0, zweight=1,
  use_triggers, level, custommask=0;
double       *model, *contrasts=NULL, *norm_contrasts=NULL, tsmean, tmpf, maxz,
  *TS_model, *TS_copemodel, *TS_pemodel, *TS_data, *TS_residuals;
float        *triggers;
char         rofpM[100000], rofpF[100000], rofpP[100000], 
  fmridata[10000], fsldir[10000], featdir[10000], outputdir[10000], gpstart[10000], gpname[10000], gprootname[10000],
  thestring[10000], datafile[10000], statname[100], gnuplotstr[10000], convertstr[10000], vname[100];
ColumnVector pwts;

volume<float> immask;

rofpM[0]='\0';

/* }}} */

  /* {{{ process arguments */

if (argc<2) usage();

strcpy(featdir,argv[argi++]);

strcpy(outputdir,featdir);

sprintf(fmridata,"%s/filtered_func_data",featdir);

/* use environment variables to work out where "gnuplot" and "convert" are, for later system calls */
sprintf(fsldir,getenv("FSLDIR"));
sprintf(gnuplotstr,getenv("FSLGNUPLOT"));
sprintf(convertstr,getenv("FSLCONVERT"));

for (;argi<argc;argi++)
{
  if (!strcmp(argv[argi], "-f"))
    /* {{{ alternative fmri data */

{
  argi++;
  if (argc<argi+1)
    {
      printf("Error: no value given following -f\n");
      usage();
    }
  strcpy(fmridata,argv[argi]);
}

/* }}} */
  else if (!strcmp(argv[argi], "-c"))
    /* {{{ alternative voxel position */

{
  coordset=1;
  argi++;
  if (argc<argi+3) /* options following c haven't been given */
    {
      printf("Error: incomplete values given following -c\n");
      usage();
    }
  X=atoi(argv[argi++]);
  Y=atoi(argv[argi++]);
  Z=atoi(argv[argi]);
}

/* }}} */
  else if (!strcmp(argv[argi], "-C"))
    /* {{{ output data only */

{
  coordset=1;
  dataonly=1;
  argi++;
  if (argc<argi+4)
    {
      printf("Error: incomplete values given following -C\n");
      usage();
    }
  X=atoi(argv[argi++]);
  Y=atoi(argv[argi++]);
  Z=atoi(argv[argi++]);
  strcpy(datafile,argv[argi]);
}

/* }}} */
  else if (!strcmp(argv[argi], "-m"))
    /* {{{ alternative mask image */

{
  custommask=1;

  argi++;
  if (argc<argi+1)
    {
      printf("Error: no mask image given following -m\n");
      usage();
    }

  if ( read_volume(immask,argv[argi]) )
    {
      printf("Error: mask image chosen doesn't exist\n");
      usage();
    }
}

/* }}} */
  else if (!strcmp(argv[argi], "-o"))
    /* {{{ output dir */

{
  argi++;
  if (argc<argi+1)
    {
      printf("Error: no value given following -o\n");
      usage();
    }
  strcpy(outputdir,argv[argi]);
}

/* }}} */
  else if (!strcmp(argv[argi], "-n"))
    /* {{{ zweight clusters? */

{
  zweight=0;
}

/* }}} */
  else if (!strcmp(argv[argi], "-p"))
    /* {{{ turn on prewhitening */

{
  prewhiten=1;
}

/* }}} */
}

/* }}} */
  /* {{{ read filtered_func_data */

volume4D<int> im;
read_volume4D(im, fmridata);

size=im.nvoxels();

if (dataonly)
  /* {{{ output raw data and exit */

{
  if((ofp=fopen(datafile,"wb"))==NULL)
    {
      fprintf(stderr,"Can't open output data file %s\n",datafile);
      exit(1);
    }

  for(t=0; t<im.tsize(); t++)
    fprintf(ofp,"%e\n",(float)im(X,Y,Z,t));

  fclose(ofp);

  return 0;
}

/* }}} */

/* }}} */
  /* {{{ read design.mat */

sprintf(thestring,"%s/design.mat",featdir);
model=read_model(thestring,&nevs,&npts);

if (npts==0)
{
  modelfree=1;
  nftests=1;
  ncon=0;
}

npts=im.tsize();

TS_model     = (double *)malloc(sizeof(double)*npts);
TS_copemodel = (double *)malloc(sizeof(double)*npts);
TS_pemodel   = (double *)malloc(sizeof(double)*npts*nevs);
TS_data      = (double *)malloc(sizeof(double)*npts);
TS_residuals = (double *)malloc(sizeof(double)*npts);

/* }}} */
  /* {{{ read auto correlation estimates for prewhitening */

double* pwmodel=0;
volume4D<float> acs;

if ( prewhiten ) {

  prewhiten=0;

  sprintf(thestring,"%s/stats/threshac1",featdir);
  if (fsl_imageexists(string(thestring)))
    {
      read_volume4D(acs, thestring);

      if (acs[1].max()!=0) /* hacky test for whether prewhitening was actually carried out */
	{
	  pwmodel=new double[nevs*npts];
	  prewhiten=1;
	}
    }

}

/* }}} */
  /* {{{ read design.con and PEs */

vector< volume<float> > impe(nevs);
if (!modelfree)
{
  sprintf(thestring,"%s/design.con",featdir);
  contrasts=read_contrasts(thestring,&nevs,&ncon);

  /* create normalised contrasts */
  norm_contrasts=(double*)malloc(sizeof(double)* nevs * ncon);
  for(i=0; i<ncon; i++)
    {
      double norm_factor=0;

      for(ev=0; ev<nevs; ev++)
	norm_factor += contrasts[i*nevs+ev] * contrasts[i*nevs+ev];

      for(ev=0; ev<nevs; ev++)
	norm_contrasts[i*nevs+ev] = contrasts[i*nevs+ev] / sqrt(norm_factor);
    }

  for(i=1;i<=nevs;i++)
    {
      sprintf(thestring,"%s/stats/pe%d",featdir,i);
      read_volume(impe[i-1], thestring);
    }
}

/* }}} */
  /* {{{ read design.fts */

if (!modelfree)
{
  sprintf(thestring,"%s/design.fts",featdir);
  read_ftests(thestring,&nftests);
}

/* }}} */
  /* {{{ read triggers */

use_triggers=read_triggers(featdir,&triggers,nevs,npts);

/*
for(ev=0;ev<nevs;ev++)
{
  for(i=0;i<=triggers[ev]+1;i++) printf("%f ",triggers[i*nevs+ev]);
  printf("\n");
}
*/

/* }}} */
  /* {{{ check analysis level */

level=1;

sprintf(thestring,"%s/design.lev",featdir);
if((ifp=fopen(thestring,"rb"))!=NULL)
{
  fclose(ifp);
  level=2;
}

/* }}} */
  /* {{{ create plot(s) for each contrast */

for(j=0;j<2;j++)
{
  /* {{{ setup stats type */

int maxi;

if (j==0) { sprintf(statname,"zstat");  maxi=ncon; }
else      { sprintf(statname,"zfstat"); maxi=nftests; }

/* }}} */

  for(i=0; i<maxi; i++)
    {
      /* {{{ vars */

volume<float> imcope, imz, imweight;
bool haveclusters=false;

rofpF[0]='\0';
rofpP[0]='\0';

/* }}} */
      /* {{{ read COPE and derived stats; test for f-test output */

/* load zstat or zfstat */
sprintf(thestring,"%s/stats/%s%d",featdir,statname,i+1);
if (fsl_imageexists(string(thestring))) 
{
  read_volume(imz,thestring);
  imweight=imz;
  if (!zweight)
    imweight=1;
}
else
  continue; /* f-test i wasn't valid - no zfstat image */

/* load cope */
if ( (j==0) && (!modelfree) )
{
  sprintf(thestring,"%s/stats/cope%d",featdir,i+1);
  read_volume(imcope,thestring);
}

/* load cluster mask */
if (!coordset)
{
  if (!custommask)
    {
      sprintf(thestring,"%s/cluster_mask_%s%d",featdir,statname,i+1);
      if (fsl_imageexists(string(thestring)))
	read_volume(immask,thestring);
    }
  haveclusters=(immask.max()>0);
}

/* }}} */
      /* {{{ find max Z and X,Y,Z */

if (!coordset)
{
  X=Y=Z=0;
  maxz=-1000;

  for(int z=0; z<im.zsize(); z++)
    for(int y=0; y<im.ysize(); y++)
      for(int x=0; x<im.xsize(); x++)
	if ( (imz(x,y,z)>maxz) &&
	     ( (!haveclusters) ||   /* make max Z be inside a cluster if we found a cluster map */
	       (immask(x,y,z)>0) && (!prewhiten || acs(x,y,z,1)!=0 || acs(x,y,z,2)!=0) ) )
	  {
	    maxz=imz(x,y,z);
	    X=x; Y=y; Z=z;
	  }
}
else
  maxz=imz(X,Y,Z);

/* }}} */

      /* first do peak voxel plotting then do mask-averaged plotting */
      for(v=0;v<2;v++)
	{
          /* {{{ setup */

double wtotal=0;
int count=0;
volume<float> tmp_imweight=imweight, tmp_immask=immask;

if (v==0) {
  vname[0]='\0';
  tmp_imweight=0;
  tmp_imweight(X,Y,Z)=1;
  tmp_immask=tmp_imweight;
} else {
  sprintf(vname,"c");
}

/* }}} */
          /* {{{ create model and data time series */

for(t=0; t<npts; t++)
{
  TS_data[t]=TS_model[t]=TS_copemodel[t]=0;
  for(ev=0; ev<nevs; ev++)
    TS_pemodel[ev*npts+t]=0;
}

for(int x=0; x<im.xsize(); x++) for(int y=0; y<im.ysize(); y++) for(int z=0; z<im.zsize(); z++)
  if (tmp_immask(x,y,z)>0 && (!prewhiten || acs(x,y,z,1)!=0 || acs(x,y,z,2)!=0))
  {
    count++;
    wtotal+=tmp_imweight(x,y,z);

    if(prewhiten)
      prewhiten_timeseries(acs.voxelts(x,y,z), im.voxelts(x,y,z), pwts, npts);
    else
      pwts = im.voxelts(x,y,z);
    for(t=0; t<npts; t++)
      TS_data[t]+=pwts(t+1)*tmp_imweight(x,y,z);

    if (!modelfree) {
      if (prewhiten)
	prewhiten_model(acs.voxelts(x,y,z), model, pwmodel, nevs, npts);
      else
	pwmodel=model;
      for(t=0; t<npts; t++)
	for(ev=0; ev<nevs; ev++)
	  {
	    tmpf=pwmodel[t*nevs+ev]*impe[ev](x,y,z)*tmp_imweight(x,y,z);
	    TS_model[t]           += tmpf;
	    TS_copemodel[t]       += tmpf*norm_contrasts[i*nevs+ev];
	    TS_pemodel[ev*npts+t] += tmpf;
	  }
    }
  }

tsmean=0;
for(t=0; t<npts; t++)
{
  TS_data[t]/=wtotal;
  tsmean+=TS_data[t];
}
tsmean/=npts;

if (level==2) tsmean=0;

if (!modelfree)
  for(t=0; t<npts; t++)
    {
      TS_model[t] = TS_model[t]/wtotal + tsmean;
      TS_copemodel[t] = TS_copemodel[t]/wtotal + tsmean;
      TS_residuals[t]=TS_data[t]-TS_model[t];
      for(ev=0; ev<nevs; ev++)
	TS_pemodel[ev*npts+t] = TS_pemodel[ev*npts+t]/wtotal + tsmean;
    }

/* }}} */
	  /* {{{ output data text files */

sprintf(thestring,"%s/tsplot%s_%s%d.txt",outputdir,vname,statname,i+1);
ofp=fopen(thestring,"wb");

ymin=ymax=(int)TS_data[0];

for(t=0; t<npts; t++)
{
  fprintf(ofp,"%e ",TS_data[t]); ymin=(int)MISCMATHS::Min(TS_data[t],ymin); ymax=(int)MISCMATHS::Max(TS_data[t],ymax);

  if (!modelfree)
    {
      if (j==0)
	{
	  fprintf(ofp,"%e ",TS_copemodel[t]); ymin=(int)MISCMATHS::Min(TS_copemodel[t],ymin); ymax=(int)MISCMATHS::Max(TS_copemodel[t],ymax);
	}

      fprintf(ofp,"%e ",TS_model[t]); ymin=(int)MISCMATHS::Min(TS_model[t],ymin); ymax=(int)MISCMATHS::Max(TS_model[t],ymax);
	
      if (j==0)
	fprintf(ofp,"%e ",TS_residuals[t]+TS_copemodel[t]);
    }

  fprintf(ofp,"\n");
}

fclose(ofp);

ymax+=(ymax-ymin)/5;
ymin-=(ymax-ymin)/20;

/* }}} */
	  /* {{{ run gnuplot */

sprintf(gprootname,"tsplot%s_%s%d",vname,statname,i+1);
sprintf(gpname,"%s/%s",outputdir,gprootname);

sprintf(gpstart,"echo \"set term pbm color\nset size %f,0.4\nset title '%s%d :", MISCMATHS::Min(MISCMATHS::Max(npts*4,600),3000)/640.0,statname,i+1);

if (v==0) {
  if (!coordset)
    sprintf(gpstart,"%s max Z stat of %.1f at voxel (%d,%d,%d)",gpstart,maxz,X,Y,Z);
  else
    sprintf(gpstart,"%s Z stat of %.1f at selected voxel (%d,%d,%d)",gpstart,maxz,X,Y,Z);
} else {
    sprintf(gpstart,"%s averaged over %d voxels",gpstart,count);
}

sprintf(gpstart,"%s'\nset yrange [%d:%d]\nplot '%s.txt' using",gpstart,ymin,ymax,gpname);

if (!modelfree)
{
  if (j==0)
    {
      sprintf(thestring,"%s 1 title 'data' with lines , '%s.txt' using 2 title 'cope partial model fit' with lines , '%s.txt' using 3 title 'full model fit' with lines \" | %s | %s - %s.gif",gpstart,gpname,gpname,gnuplotstr,convertstr,gpname);
      system(thestring);
      sprintf(thestring,"%s 4 title 'reduced data' with lines , '%s.txt' using 2 title 'cope partial model fit' with lines \" | %s | %s - %sp.gif",gpstart,gpname,gnuplotstr,convertstr,gpname);
      system(thestring);
      sprintf(rofpF,"%sFull model fit - <a href=\"%sp.gif\">Partial model fit</a> - <a href=\"%s.txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\"%s.gif\"><br><br>\n",rofpF,gprootname,gprootname,gprootname);
    }
  else
    {
      sprintf(thestring,"%s 1 title 'data' with lines , '%s.txt' using 2 title 'full model fit' with lines \" | %s | %s - %s.gif",gpstart,gpname,gnuplotstr,convertstr,gpname);
      system(thestring);
      sprintf(rofpF,"%sFull model fit - <a href=\"%s.txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\"%s.gif\"><br><br>\n",rofpF,gprootname,gprootname);
    }
}
else
{
  sprintf(thestring,"%s 1 title 'data' with lines \" | %s | %s - %s.gif",gpstart,gnuplotstr,convertstr,gpname);
  system(thestring);
  sprintf(rofpF,"%sData plot - <a href=\"%s.txt\">Raw data</a>\n<IMG BORDER=0 SRC=\"%s.gif\"><br><br>\n",rofpF,gprootname,gprootname);
}

/* picture for main web index page */
if (v==0)
     sprintf(rofpM,"%s<a href=\"%s.html\"><IMG BORDER=0 SRC=\"%s.gif\"></a><br><br>\n",rofpM,gprootname,gprootname);

/* }}} */
  	  /* {{{ peri-stimulus: output the data & run gnuplot */

if (use_triggers)
{
  if (!modelfree)
    sprintf(rofpP,"%s<table><tr>\n",rofpP);

  for(ev=0; ev<nevs; ev++)
    if (triggers[ev]>0.5)
{
  int which_event;
  float ps_period=triggers[((int)triggers[ev]+1)*nevs+ev];

  sprintf(gprootname,"ps_tsplot%s_%s%d_ev%d",vname,statname,i+1,ev+1);

  sprintf(thestring,"%s/%s.txt",outputdir,gprootname);
  ofp=fopen(thestring,"wb");

  for(which_event=1;which_event<=triggers[ev];which_event++)
    {
      double min_t=triggers[which_event*nevs+ev];
      int int_min_t=(int)ceil(min_t),
	max_t=MISCMATHS::Min(npts-1,int_min_t+(int)ps_period);

      for(t=int_min_t;t<max_t;t++)
	{
	  fprintf(ofp,"%.1f ",t-min_t); /* time (restricted temporal accuracy (0.1*TR) to make plotting slightly smoother ie force binning) */

	  fprintf(ofp,"%e ",TS_residuals[t]+TS_model[t]);                                    /* full data (full model + residuals) */

	  if (!modelfree)
	    {
	      fprintf(ofp,"%e ",TS_pemodel[ev*npts+t]);                  /* model - this EV only */
	      fprintf(ofp,"%e ",TS_model[t]);                                                /* full model */
	      fprintf(ofp,"%e ",TS_residuals[t]+TS_pemodel[ev*npts+t]);  /* data for this EV only (EV + residuals) */
	    }

	  fprintf(ofp,"\n");
	}
    }
  
  fclose(ofp);

  sprintf(gpname,"%s/%s",outputdir,gprootname);
 
  sprintf(gpstart,"echo \"set term pbm color\nset xlabel 'peristimulus time (TRs)'\nset size %f,0.4\nset title '%s%d ev%d :", MISCMATHS::Min(MISCMATHS::Max(ps_period*3,400),3000)/640.0,statname,i+1,ev+1);

  if (v==0) {
    if (!coordset)
      sprintf(gpstart,"%s max Z stat of %.1f at voxel (%d,%d,%d)",gpstart,maxz,X,Y,Z);
    else
      sprintf(gpstart,"%s Z stat of %.1f at selected voxel (%d,%d,%d)",gpstart,maxz,X,Y,Z);
  } else {
    sprintf(gpstart,"%s averaged over %d voxels",gpstart,count);
  }

  sprintf(gpstart,"%s'\nset yrange [%d:%d]\nplot '%s.txt' using",gpstart,ymin,ymax,gpname);

  if (!modelfree)
    {
      sprintf(thestring,"%s 1:2 title 'data' , '%s.txt' using 1:3 smooth unique title 'EV %d model fit' , '%s.txt' using 1:4 smooth unique title 'full model fit' \" | %s | %s - %s.gif",gpstart,gpname,ev+1,gpname,gnuplotstr,convertstr,gpname);
      system(thestring);
      sprintf(thestring,"%s 1:5 title 'reduced data' , '%s.txt' using 1:3 smooth unique title 'EV %d model fit' \" | %s | %s - %sp.gif",gpstart,gpname,ev+1,gnuplotstr,convertstr,gpname);
      system(thestring);
      sprintf(rofpP,"%s<td>Full model fit - <a href=\"%sp.gif\">Partial model fit</a> - <a href=\"%s.txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\"%s.gif\">\n",rofpP,gprootname,gprootname,gprootname);
    }
  else
    {
      sprintf(thestring,"%s 1:2 title 'data' with lines \" | %s | %s - %s.gif",gpstart,gnuplotstr,convertstr,gpname);
      system(thestring);
      sprintf(rofpP,"%sData plot - <a href=\"%s.txt\">Raw data</a>\n<IMG BORDER=0 SRC=\"%s.gif\"><br><br>\n",rofpP,gprootname,gprootname);
    }
}

  if (!modelfree)
    sprintf(rofpP,"%s</tr></table><br><br>\n",rofpP);

}

/* }}} */

  	  if (!haveclusters) v=10;
        }

      /* {{{ web output */

sprintf(thestring,"%s/tsplot_%s%d.html",outputdir,statname,i+1);

if((rofp=fopen(thestring,"wb"))==NULL)
{
  fprintf(stderr,"Can't open output report file %s\n",outputdir);
  exit(1);
}

fprintf(rofp,"<HTML>\n<TITLE>%s%d</TITLE>\n<BODY BACKGROUND=\"file:%s/doc/images/fsl-bg.jpg\">\n<hr><CENTER>\n<H1>FEAT Time Series Report - %s%d</H1>\n</CENTER>\n<hr><b>Full plots</b><p>\n%s\n<hr><b>Peristimulus plots</b><p>\n%s\n<HR></BODY></HTML>\n\n",statname,i+1,fsldir,statname,i+1,rofpF,rofpP);

fclose(rofp);

/* }}} */
    }
}

/* {{{ main web index page output */

/* first output full index page (eg for use by featquery) */

sprintf(thestring,"%s/tsplot_index.html",outputdir);

if((rofp=fopen(thestring,"wb"))==NULL)
{
  fprintf(stderr,"Can't open output report file %s\n",outputdir);
  exit(1);
}

fprintf(rofp,"<HTML>\n<TITLE>FEAT Time Series Report</TITLE>\n<BODY BACKGROUND=\"file:%s/doc/images/fsl-bg.jpg\">\n<hr><CENTER>\n<H1>FEAT Time Series Report</H1>\n</CENTER>\n<hr>%s<HR></BODY></HTML>\n\n",fsldir,rofpM);

fclose(rofp);


/* now output same thing without start and end, for inclusion in feat report */

sprintf(thestring,"%s/tsplot_index",outputdir);

if((rofp=fopen(thestring,"wb"))==NULL)
{
  fprintf(stderr,"Can't open output report file %s\n",outputdir);
  exit(1);
}

fprintf(rofp,"%s\n\n",rofpM);

fclose(rofp);

/* }}} */

/* }}} */

  exit(0);
}

