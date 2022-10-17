/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

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

#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include "mriseg.h"
#include "imageio.h"
#include "imalgorithm.h"
#include "classification.h"
using namespace std;

bool	b2d=false;
int		type=0;
int		nclass=3;
float	neighbor=0.3f;
int		nbiter=8;
int		nblowpass=100;
int	bapriori = 0;
string  basepathAPriori = "";
char	fnmean[1024];
bool	autopara = true;
string  outbase;
bool	outseg = true;
bool	outbias = false;
bool	outallbias = false;
bool	outBYTEbias = false;
int	nItersSpeardBias = 5;
bool	outrestore = false;
bool	outsegs = false;
bool	outpbmap = false;
bool	outpve = false;

bool	pve = false;
bool	pvepara = true;
float	pveb = 1.5f;

string	VERSION = "3.53";

int		verbose = 1;

void usage(bool full=false)
{
	cerr << "\nFAST - FMRIB's Automated Segmentation Tool\n";
	cerr << "FAST Version " << VERSION << "\n\n";

	cerr << "Usage: fast [-options] imagefile \n";

	cerr << "options:\n";
	cerr << "    -t<n>:\tinput image type. 1-T1, 2-T2, 3-PD (default T1)\n";
	cerr << "    -c <n>:\tnumber of classes (default 3)\n";
#ifndef WIN32
	cerr << "    -a: \tuse apriori probability maps (for initialisation only)\n";
	cerr << "    -A: \tuse apriori probability maps (for initialisation and for posteriors)\n";
	cerr << "    -ap <prefix> | -Ap <prefix>: \tpath+filename prefix to pre-registered (to input data) apriori probability maps, <prefix>_gm _wm and _csf (uncompressed Analyze format)\n";
#endif
	cerr << "    -od <base>:\toutput basename (default <input>)\n";
	cerr << "    -os:\toutput segmentation with one class per image\n";
	cerr << "    -op:\toutput probability maps (one for each class)\n";
	cerr << "    -or:\toutput restored image\n";
	cerr << "    -ob:\toutput bias correction field\n";
	cerr << "    -n:\t\tdo not output segmentation\n";

	if(full)
	{
		cerr << "    -m <file>:\tuse file with starting values for the mean of each tissue\n";
		cerr << "    -2:\t\t2D segmentation (only for 3D images)\n";
		cerr << "    -i <n>:\tnumber of iterations (default 8)\n";
		cerr << "    -l <n>:\titerations for bias correction field smoothing (default 100)\n";
		cerr << "    -p:\t\tdisable automatic parameter updating\n";
		cerr << "    -b <v>:\tMRF neighbourhood beta value (default 0.3)\n";
		cerr << "    -oba <n>:\toutput dilated bias correction field (n iterations)\n";
		cerr << "    -e:\t\tenable partial volume classification\n";
		cerr << "    -ov:\toutput partial volume images (one for each class)\n";
		cerr << "    --b <v>:\tMRF beta value for PVE classification (default 1.5)\n";
		cerr << "    --p:\tdisable parameter updating during PVE classification\n";
		cerr << "    -v<n>:\tverbose level (0-5; default 1)\n\n";
	}
	else 
	{
		cerr << "    -v<n>:\tverbose level (0-5; default 1)\n";
		cerr << "    -h:\t\tadvanced options\n\n";
	}

	exit(0);
}

void ParameterAnalyze(int argc, char** argv)
{
	if(argc < 2) usage(false);
	if(argv[argc-1][0] == '-')
		if(argv[argc-1][1] == 'h') usage(true);
		else usage(false);

	fnmean[0] = 0;
	int argindex = 1;
	while (argindex < argc-1)
	{
		char* tcp = argv[argindex];
		if (*tcp == '-')
		{
			switch (*++tcp)
			{
			case 'h': // full usage
				usage(true);
			  	break;

			case 'e': // pve
				pve = true;
			  	break;

			case 'c': // number of classes
				nclass = atoi(argv[++argindex]);
				if(nclass < 2 || nclass > 6) 
				{
					nclass = 3;
					cerr << "number of classes can only be between 2 and 6 and is set to 3\n";
				}
			  	break;

			case '2': // 2D segmentation
				b2d = true;
				break;

#ifndef WIN32
			case 'a': // apriori 
				bapriori = 1;
				if ( (*++tcp)=='p' )
 				  basepathAPriori = argv[++argindex];
				break;
			case 'A': // apriori including posteriors
				bapriori = 2;
				if ( (*++tcp)=='p' )
 				  basepathAPriori = argv[++argindex];
				break;
#endif

			case 't':
				type = *++tcp - '1';
				if(type < 0 || type > 2) type = 0;
				break;

			case 'b': // MRF neighbourhood value
				neighbor = float(atof(argv[++argindex]));
				if(neighbor > 3 || neighbor < 0) neighbor = 0.3f;
				break;

			case 'i': // number of iterations
				nbiter = atoi(argv[++argindex]);
				break;

			case 'l': // number of lowpass filtering
				nblowpass = atoi(argv[++argindex]);
				if(nblowpass < 0) nblowpass = 100;
				break;

			case 'm': // mean parameter file
				strcpy(fnmean, argv[++argindex]);
				break;
			
			case 'p': 
				autopara = false;
				break;

			case 'n': // No segmentation output
				outseg = false;
				break;

			case 'o': // output result
				switch(*++tcp)
				{
				case 'd':		// output basename
					outbase = argv[++argindex];
					break;
				case 's':		//individual class segmenation image
					outsegs = true;
					break;
				case 'b':		// bias field
					++tcp;
					if(*tcp == 'a')
					{
						outallbias = true;	//spreaded bias field
						++argindex;
						if(argv[argindex][0] < '0' || argv[argindex][0] > '9')
						{
							cerr << "number of iterations for spreading bias field is required!\n"; 
							usage();
						}
						nItersSpeardBias = atoi(argv[argindex]);	
					}
					if(*tcp == 'b') outBYTEbias = true;
					else outbias = true;
					break;
				case 'r':		//restored image
					outrestore = true;
					break;
				case 'p':		//probability map
					outpbmap = true;
					break;
				case 'v':		//pve
					outpve = true;
					break;
				}
				break;
			
			case '-': // pve options
				if(*++tcp == 'p') pvepara = false;
				else if(*tcp == 'b') pveb = float(atof(argv[++argindex]));
				if(pveb > 3 || pveb < 0) pveb = 1.0f;
			  	break;

			case 'v': // density of the extra uniform distributed class
				verbose = *++tcp - '0';
				if(verbose < 0 || verbose > 5) verbose = 1;
				break;

			default: // unknown argument
				cerr << "Unknow argument: -" << *tcp << "\n\n";
				usage(true);
				break;
			}	    
		}
		argindex++;
	}

	if(fnmean[0] != 0 && nclass == 2)
	{
		cerr << "Using mean input requires 3 or more classes!\n";
		exit(0);
	}

	if(!pve && outpve) pve = true;
}

void SpreadBiasField(ZGrayFloatImage& Bias, int nIter)
{
	ZGrayFloatImage tmp = Bias;
	int width = Bias.Width();
	int height = Bias.Height();
	int depth = Bias.Depth();
	int slicesize = Bias.PixelsPerSlice();
	int offset3d = -slicesize - width - 1;
	int offset2d = - width - 1;
	bool m_b3D = false;
	
	if(depth > 2)
	{
		Bias.SetROI(1, width-2, 1, height-2, 1, depth-2);
		tmp.SetROI(1, width-2, 1, height-2, 1, depth-2);
		m_b3D = true;
	}
	else
	{
		Bias.SetROI(1, width-2, 1, height-2, 0, depth-1);
		tmp.SetROI(1, width-2, 1, height-2, 0, depth-1);
		m_b3D = false;
	}

	ZImage<float>::iterator	p_bias(Bias), p_tmp(tmp);

	for(int iter=0; iter<nIter; iter++)
	{
		for (p_bias.reset(), p_tmp.reset(); p_bias.more(); p_bias++, p_tmp++)
		{
			if(*p_bias != 1) continue;

			float average=0;
			int num=0;

			if(m_b3D)
			{
				float* src = p_bias + offset3d;
				for(int k=0; k<3; k++, src+=slicesize)
				{
					float* ptr = src;
					for(int j=0; j<3; j++, ptr+=width)
					for(int i=0; i<3; i++)
					{
						if(ptr[i] != 1) average += ptr[i], num++;
					}
				}
			}
			else 
			{
				float* ptr = p_bias + offset2d;
				for(int j=0; j<3; j++, ptr+=width)
				for(int i=0; i<3; i++)
				{
					if(ptr[i] != 1) average += ptr[i], num++;
				}
			}

			if(num==0) continue;

			average /= num;

			*p_tmp = average;
		}
		Bias.Replace(tmp);
	}
	Bias.FullROI();
}

int main(int argc, char* argv[])
{
	ParameterAnalyze(argc, argv);

	if(verbose) 
		cerr << "\nFAST - FMRIB's Automated Segmentation Tool\t\tVersion " << VERSION << "\n\n";

	float pixdim[3];
	ZImageReader imagereader;
	ZMRISegmentation	mri;
	ZImageBase* inputimage;

	if(verbose) cerr << "Image: " << argv[argc-1] << endl << endl;

	string ext = ".hdr"; 
	string avwname = string(argv[argc-1])+ext;

	fstream file(avwname.c_str(), ios::in);
	if(!file)
	{
		char* cext = strrchr(argv[argc-1], '.');
		if(outbase.empty()) outbase = BaseName(argv[argc-1]);

		if(cext==0 || cext > argv[argc]-5 ) 
		{ 
			cerr << "Unknown file format!\n";
			return -1;
		}
		imagereader.SetFile(argv[argc-1]); ext = cext;

		if((inputimage = imagereader.ReadFile()) == 0) return -1;
	}
	else
	{
		file.close();
		if(outbase.empty()) outbase = argv[argc-1];
		imagereader.SetFile(argv[argc-1], AVW);
		if((inputimage = imagereader.ReadFile()) == 0)
		{
			cerr << "Can not open file " << avwname << "!\n";
			return -1;
		}
	}

	ZVector<float> mean(nclass);
	if(fnmean[0] != 0)
	{
		fstream file(fnmean, ios::in);
		if(!file) 
		{
			cerr << "Can not open file " << fnmean << " to read!\n";
			return -1;
		}

		ZVector<float>::iterator p_m = mean.pbegin();
		//p_m = (float*)mean.begin();
		for(; p_m != mean.pend(); p_m++) file >> *p_m;
	}

	inputimage->GetPixelDim(pixdim[0], pixdim[1], pixdim[2]);

	int width = inputimage->Width();
	int height = inputimage->Height();
	int depth = inputimage->Depth();
	int size = width * height * depth;

	if(verbose) 
	{
		switch(type)
		{
		default:
		case 0:
			cerr << "T1-weighted image" << endl;
			break;
		case 1:
			cerr << "T2-weighted image" << endl;
			break;
		case 2:
			cerr << "PD-weighted image" << endl;
			break;
		}

		cerr << "Imagesize : " << width << " x " << height << " x " << depth << endl;
		cerr << "Pixelsize : " << pixdim[0] << " x " << pixdim[1] << " x "  << pixdim[2] << endl << endl;
	}

	mri.Create(*inputimage, type, nclass, b2d, nbiter, nblowpass, neighbor, 2, autopara, false);

	time_t T1, T2;
	T1 = time(NULL);
	mri.SetVerboseMode(verbose);

	ZGrayByteImage *pGM=0, *pWM=0, *pCSF=0;

#ifndef WIN32
	if(bapriori>0 && nclass>3)
	{
		cerr << "Apriori can only be used for 3-class segmentation\n";
		bapriori = 0;
	}

	if(bapriori>0 && depth==1)
	{
		cerr << "Apriori can only be used for 3-D images\n";
		bapriori = 0;
	}

	string brainpath, GMpath, WMpath, CSFpath;
 	if((bapriori>0)&&basepathAPriori.empty())
	{
		bool found = false;
		string basepath(getenv("FSLDIR"));
		if(!basepath.empty())
		{
			brainpath = basepath + "/etc/standard/avg152T1_brain.hdr";
			fstream file(brainpath.c_str(), ios::in);
			if(!(!file))
			{
				found = true;
				basepath += "/etc/standard/";
			}
		}
		if(!found)
		{
			basepath = ""; brainpath = basepath + "avg152T1_brain.hdr";
			fstream file(brainpath.c_str(), ios::in);
			if(!file)
			{
				if(verbose) cerr << "prior images avg152T1_brain are not found! priors are not used!\n";
				bapriori = 0;
			}
		}
			
		if(bapriori>0)
		{
			GMpath = basepath + "avg152T1_gray.hdr";
			fstream gray(GMpath.c_str(), ios::in);
			WMpath = basepath + "avg152T1_white.hdr";
			fstream white(WMpath.c_str(), ios::in);
			CSFpath = basepath + "avg152T1_csf.hdr";
			fstream csf(CSFpath.c_str(), ios::in);
			if(!gray || !white || !csf)
			{
				if(verbose) cerr << "prior images avg152T1_brain are not found! priors are not used!\n";
				bapriori = 0;
			}
		}
	}

	if(bapriori>0)
	{
	  if (!basepathAPriori.empty())
	    {
	      if(verbose) cerr << "reading pre-registered probability maps...\n";
	      imagereader.SetFile((basepathAPriori+"_gm.hdr").c_str());
	      
	      pGM = new ZGrayByteImage;
	      pWM = new ZGrayByteImage;
	      pCSF = new ZGrayByteImage;
	      if(!imagereader.ReadFile(*pGM))
 		{
		  cerr << "Can not read apriori GM probability image gray.tif" << endl;
		  delete inputimage;	return -1;
 		}
	      imagereader.SetFile((basepathAPriori+"_wm.hdr").c_str());
	      if(!imagereader.ReadFile(*pWM))
 		{
		  cerr << "Can not read apriori WM probability image white.tif" << endl;
		  delete inputimage;	return -1;
 		}
	      imagereader.SetFile((basepathAPriori+"_csf.hdr").c_str());
	      if(!imagereader.ReadFile(*pCSF))
 		{
		  cerr << "Can not read apriori CSF probability image csf.tif" << endl;
		  delete inputimage;	return -1;
 		}
	    }
	  else
	    {
	        if(verbose) cerr << "registering the probability maps to the image...\n";
		char reg[1024];
		sprintf(reg, "%s/bin/flirt -omat %s -in %s -ref %s", getenv("FSLDIR"), (outbase+"in2st.mat").c_str(), argv[argc-1], brainpath.c_str());
		if(verbose) cerr << reg << endl;
		system(reg);
		sprintf(reg, "%s/bin/convert_xfm -inverse -omat %s %s", getenv("FSLDIR"), (outbase+"st2in.mat").c_str(), (outbase+"in2st.mat").c_str());
		if(verbose) cerr << reg << endl;
		system(reg);
		sprintf(reg, "%s/bin/flirt -init %s -applyxfm -o %s -in %s -ref %s", getenv("FSLDIR"), (outbase+"st2in.mat").c_str(), (outbase+"gm.hdr").c_str(), GMpath.c_str(), argv[argc-1]);
		if(verbose) cerr << reg << endl;
		system(reg);
		sprintf(reg, "%s/bin/flirt -init %s -applyxfm -o %s -in %s -ref %s", getenv("FSLDIR"), (outbase+"st2in.mat").c_str(), (outbase+"wm.hdr").c_str(), WMpath.c_str(), argv[argc-1]);
		if(verbose) cerr << reg << endl;
		system(reg);
		sprintf(reg, "%s/bin/flirt -init %s -applyxfm -o %s -in %s -ref %s", getenv("FSLDIR"), (outbase+"st2in.mat").c_str(), (outbase+"csf.hdr").c_str(), CSFpath.c_str(), argv[argc-1]);
		if(verbose) cerr << reg << endl;
		system(reg);

		imagereader.SetFile((outbase+"gm.hdr").c_str());

		pGM = new ZGrayByteImage;
		pWM = new ZGrayByteImage;
		pCSF = new ZGrayByteImage;
		if(!imagereader.ReadFile(*pGM))
		{
			cerr << "Can not read apriori GM probability image gray.tif" << endl;
			delete inputimage;	return -1;
		}
		imagereader.SetFile((outbase+"wm.hdr").c_str());
		if(!imagereader.ReadFile(*pWM))
		{
			cerr << "Can not read apriori WM probability image white.tif" << endl;
			delete inputimage;	return -1;
		}
		imagereader.SetFile((outbase+"csf.hdr").c_str());
		if(!imagereader.ReadFile(*pCSF))
		{
			cerr << "Can not read apriori CSF probability image csf.tif" << endl;
			delete inputimage;	return -1;
		}

		remove((outbase+"in2st.mat").c_str());
		remove((outbase+"st2in.mat").c_str());
		remove((outbase+"gm.hdr").c_str());
		remove((outbase+"wm.hdr").c_str());
		remove((outbase+"csf.hdr").c_str());
		remove((outbase+"gm.img").c_str());
		remove((outbase+"wm.img").c_str());
		remove((outbase+"csf.img").c_str());
	    }
	}
#endif
	
	if(fnmean[0] != 0) mri.SetMeans(mean.pbegin());
	
	if(pve) mri.SetPVEPara(pvepara, pveb);
	if(mri.Segment(pve, pGM, pWM, pCSF, bapriori) == false) { delete inputimage; return -1; }

	T2 = time(NULL);

	if(verbose) 
	{
		cerr << "\nSegmentation done successfully!\n" << endl;
		cerr << "Calculation time " << difftime(T2, T1) << " seconds!\n" << endl;
	}

	bool ret = false;

	if(outseg)
	{
		string outext = depth>1 ? ext : ".png";
		string fname = outbase + "_seg" + outext;

		imagereader.SetFile(fname.c_str());
	
		// new code by SS to fix numbering of outputs for T2
		ZGrayByteImage gboutput(width, height, depth);
		gboutput.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
		ZImage<BYTE>::iterator p_out(gboutput), p_seg(mri.m_Segment);
		for(p_out.reset(), p_seg.reset(); p_seg.more(); p_out++, p_seg++)
		  {
		    *p_out=*p_seg;
		    if ((type==1)&&(nclass>2))
		      { if (*p_seg==1) *p_out=3; if (*p_seg==3) *p_out=1; }
		  }
		ret = imagereader.WriteFile(gboutput);

		//ret = imagereader.WriteFile(mri.m_Segment);
		if(verbose)
			if(ret) cerr << "Write segmentation image " << fname << " Successfully!" << endl;
			else cerr << "Can not write segmentation image " << fname << endl;
	}

	if(outrestore)
	{
		ZImageBase* pRes = mri.GetRestored();
		string fname = outbase + "_restore" + ext;
		imagereader.SetFile(fname.c_str());

		ret = imagereader.WriteFile(*pRes);
		if(verbose)
			if(ret) cerr << "Write restored image " << fname << " Successfully!" << endl;
			else cerr << "Can not write restored image " << fname << endl;
	}

	if(outbias)
	{
		string fname = outbase + "_bias.hdr";

		imagereader.SetFile(fname.c_str());

		ret = imagereader.WriteFile(mri.m_FinalBias);
		if(verbose)
			if(ret) cerr << "Write bias field image " << fname << " Successfully!" << endl;
			else cerr << "Can not write bias field image " << fname << endl;
	}

	if(outBYTEbias)
	{
		ZGrayByteImage gboutput;
		mri.MakeBiasDisplay(gboutput);
		string fname = outbase + "_bbias" + ext;

		imagereader.SetFile(fname.c_str());
		ret = imagereader.WriteFile(gboutput);
		if(verbose)
			if(ret) cerr << "Write 8-bit bias field image " << fname << " Successfully!" << endl;
			else cerr << "Can not write 8-bit bias field image " << fname << endl;
	}

	if(outallbias)
	{
		string fname = outbase + "_abias.hdr";

		imagereader.SetFile(fname.c_str());

		SpreadBiasField(mri.m_FinalBias, nItersSpeardBias);

		ret = imagereader.WriteFile(mri.m_FinalBias);
		if(verbose)
			if(ret) cerr << "Write spreaded bias field image " << fname << " Successfully!" << endl;
			else cerr << "Can not write spreaded bias field image " << fname << endl;
	}

	if(outsegs)
	{
		ZGrayByteImage gboutput(width, height, depth);
		gboutput.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
		
		ZImage<BYTE>::iterator p_out(gboutput), p_seg(mri.m_Segment);
		string outext = depth>1 ? ext : ".png";

		for(int i=0; i<nclass; i++)
		{
			char tmp[10];
			sprintf(tmp, "_seg_%d", i);
			string fname = outbase + tmp + outext;

			imagereader.SetFile(fname.c_str());

			int j=i; if ((type==1)&&(nclass>2)) { if (i==0) j=2; if (i==2) j=0; } // fix numbering of outputs for T2

			for(p_out.reset(), p_seg.reset(); p_seg.more(); p_out++, p_seg++)
			{
				if(*p_seg == j+1) *p_out = 1;
				else *p_out = 0;
			}

			ret = imagereader.WriteFile(gboutput);
			if(verbose)
				if(ret) cerr << "Write segmentation image of tissue " << i << " Successfully!" << endl;
				else cerr << "Can not write segmentation image of tissue " << i << ' ' << fname << endl;
		}
	}

	if(outpbmap)
	{
		for(int i=0; i<nclass; i++)
		{
			char tmp[10];
			sprintf(tmp, "_pbmap_%d", i);
			string fname = outbase + tmp + ext;

			imagereader.SetFile(fname.c_str());

			ZGrayFloatImage pbmap(width, height, depth);

			int j=i; if ((type==1)&&(nclass>2)) { if (i==1) j=2; if (i==2) j=1; } // fix numbering of outputs for T2

			memcpy(pbmap.GetBuffer(), mri.m_post[j], pbmap.MemorySize());
			pbmap.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	
			ret = imagereader.WriteFile(pbmap);
			if(verbose)
				if(ret) cerr << "Write probability map of tissue " << i << " Successfully!" << endl;
				else cerr << "Can not write probability map of tissue " << i << ' ' << fname << endl;
		}
	}

	if(outpve)
	{
		for(int i=0; i<nclass; i++)
		{
			char tmp[10];
			sprintf(tmp, "_pve_%d", i);
			string fname = outbase + tmp + ext;

			imagereader.SetFile(fname.c_str());

			int j=i; if ((type==1)&&(nclass>2)) { if (i==1) j=2; if (i==2) j=1; } // fix numbering of outputs for T2

			ret = imagereader.WriteFile(mri.GetPVE(j));
			if(verbose)
				if(ret) cerr << "Write PVE image of tissue " << i << " Successfully!" << endl;
				else cerr << "Can not write PVE image of tissue " << i << ' ' << fname << endl;
		}
	}

	if(verbose && (nbiter>0 || pve))
	{
		if(nclass < 4)
		{
			cerr << endl;
			cerr << "Class:\t\tCSF\t";
			for(int i=1; i<nclass; i++) cerr << "\ttissue " << i;
			cerr << "\tbrain percentage\nVolumes:\t";

			float vol[3] = {0, 0, 0}, total=0;

			int c;
			for(c=0; c<nclass; c++)
			{
				if(pve)	for(int i=0; i<size; i++) vol[c] += mri.m_PVE(i, c);
				else for(int i=0; i<size; i++) vol[c] += mri.m_post[c][i];
			}

			for(c=0; c<nclass; c++) vol[c] *= pixdim[0] * pixdim[1] * pixdim[2];

			if(nclass == 2)
			{
				total = vol[0] + vol[1];
				printf("%-10.1f\t%-10.1f\t", vol[0], vol[1]);
			}
			else
			{
				total = vol[0] + vol[1] + vol[2];
				printf("%-10.1f\t%-10.1f\t%-10.1f\t", vol[0], vol[1], vol[2]);
			}
                        fflush(stdout);
			if(nclass==2) cerr << vol[1] / total << endl;
			else cerr << (vol[1]+vol[2]) / total << endl;
		}
	}

	delete inputimage;

	return 0;
}


