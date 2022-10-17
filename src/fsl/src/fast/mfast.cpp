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
#include <strstream>
#include <string.h>
#include "imageio.h"
#include "imalgorithm.h"
#include "vclassification.h"
#include "mmriseg.h"
using namespace std;

bool	b2d=false;
int		nchannel=2;
int		nclass=3;
float	beta=0.3f;
int		nbiter=4;
int		sampleratio=2;
int		nblowpass=100;
bool	bbias = true;
int	bapriori = 0;
string  basepathAPriori = "";
string	outbase;
bool	outseg = true;
bool	outbias = false;
bool	outrestore = false;
bool	outsegs = false;
bool	outpbmap = false;
bool	outpve = false;

bool	pve = false;
bool	pvepara = true;
float	pveb = 1.5f;

string	VERSION = "3.53";

int		verbose = 1;

void usage()
{
	cerr << "\nMFAST - Multispectral FAST\n";
	cerr << "MFAST Version " << VERSION << "\n\n";

	cerr << "Usage: mfast [-options] channel1 channel2 ... \n";

	cerr << "options:\n";
	cerr << "    -s <n>\tnumber of input channels (default 2)\n";
	cerr << "    -c <n>\tnumber of classes (default 3)\n";
	cerr << "    -i <n>\tmaximum number of iterations (default 4)\n";
#ifndef WIN32
	cerr << "    -a: \tuse apriori probability maps (for initialisation only)\n";
	cerr << "    -A: \tuse apriori probability maps (for initialisation and for posteriors)\n";
        cerr << "    -ap <prefix> | -Ap <prefix>: \tpath+filename prefix to pre-registered (to input data) apriori probability maps, <prefix>_gm _wm and _csf (uncompressed Analyze format)\n";
#endif
	cerr << "    -b <v>\tMRF neighbourhood beta value (default 0.3)\n";
	cerr << "    -f\t\tdisable bias field correction\n";
	cerr << "    -n:\t\tdo not output segmentation\n";
	cerr << "    -e:\t\tenable partial volume classification\n";
	cerr << "    --b <v>:\tMRF beta value for PVE classification (default 1.5)\n";
	cerr << "    --p:\tdisable parameter updating during PVE classification\n";
	cerr << "    -od <base>:\toutput basename (default the current)\n";
	cerr << "    -os\t\toutput segmentation with one class per image\n";
	cerr << "    -ob\t\toutput bias correction field images\n";
	cerr << "    -or\t\toutput restored images\n";
	cerr << "    -op\t\toutput probability maps (one for each class)\n";
	cerr << "    -ov:\toutput partial volume images (one for each class)\n";
	cerr << "    -2:\t\t2D segmentation (only for 3D images)\n";
	cerr << "    -v<n>:\tverbose level (0-5; default 1)\n\n";

	exit(0);
}

void ParameterAnalyze(int argc, char** argv)
{
	if(argc < 2 || argv[argc-1][0] == '-') usage();

	int argindex = 1;
	while (argindex < argc-1)
	{
		char* tcp = argv[argindex];
		if (*tcp == '-')
		{
			switch (*++tcp)
			{
			case 's':
				nchannel = atoi(argv[++argindex]);
			  	break;

			case 'e': // pve
				pve = true;
			  	break;

			case 'c':
				nclass = atoi(argv[++argindex]);
				if(nclass < 2 || nclass > 6) 
				{
					nclass = 3;
					cerr << "number of classes can only be between 2 and 6 and is set to 3\n";
				}
			  	break;

			case 'f': // disable bias field estimation
				bbias = false;
				break;

			case '2': // 2D or 3D segmentation
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

			case 'b': // MRF neighbourhood value
				beta = float(atof(argv[++argindex]));
				if(beta > 3 || beta < 0) beta = 0.3f;
				break;

			case 'i': // number of iterations
				nbiter = atoi(argv[++argindex]);
				break;

			case 'l': // number of lowpass filtering
				nblowpass = atoi(argv[++argindex]);
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
					outbias = true;
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
				usage();
				break;
			}	    
		}
		argindex++;
	}

	if(bbias == false) outbias = false;

	if(outpve) pve = true;
}

int main(int argc, char* argv[])
{
	ParameterAnalyze(argc, argv);

	if(verbose) 
		cerr << "\nMFAST - Multispectral FAST\t\tVersion " << VERSION << "\n\n";

	if(nchannel < 2) 
	{
		cerr << "There must be at least 2 channels\n";
		exit(0);
	}

	ZImageReader imagereader;
	float	pixdim[3];
	int width, height, depth, SliceSize, i;

	vector<ZImageBase*>	vimages(nchannel);
	ZImageBase* inputimage = 0;

	string ext = ".hdr"; 
	string avwname = string(argv[argc-1])+ext;

	bool avwinput = true;
	fstream file(avwname.c_str(), ios::in);
	if(!file) avwinput = false, ext = strrchr(argv[argc-1], '.');
	else file.close();

	vector<string> basename(nchannel);
	for(i=nchannel; i>0; i--)
	{
		if(verbose) cerr << "Channel " << nchannel-i << ": " << argv[argc-i] << endl;
		
		if(avwinput) basename[nchannel-i] = FileName(argv[argc-i]);
		else basename[nchannel-i] = BaseName(FileName(argv[argc-i]));

		if(avwinput || ExtName(argv[argc-i])[0] == 0) 
			imagereader.SetFile(argv[argc-i], AVW);
		else 
			imagereader.SetFile(argv[argc-i]);

		if((inputimage = imagereader.ReadFile()) == 0) return -1;
		
		vimages[nchannel-i] = inputimage;
	}
	
	width = vimages[0]->Width();
	height = vimages[0]->Height();
	depth = vimages[0]->Depth();
	
	for(i=1; i<nchannel; i++)
	{
		if(vimages[i]->Width() != UINT(width) || vimages[i]->Height() != UINT(height) || vimages[i]->Depth() != UINT(depth))
		{
			cerr << "Channel " << i << " has a different size from channel 0\n";
			for(i=0; i<nchannel; i++) delete vimages[i];
			return -1;
		}
	}

	SliceSize = width * height;

	vimages[0]->GetPixelDim(pixdim[0], pixdim[1], pixdim[2]);

	time_t T1, T2;
	T1 = time(NULL);

	ZMultiMRISegmentation mri;
	mri.Create(vimages, nclass, b2d, nbiter, nblowpass, beta, 2, bbias);

	mri.SetVerboseMode(verbose);

	if(verbose) cerr << "Imagesize : " << width << " x " << height << " x " << depth << endl;
	if(verbose) cerr << "Pixelsize : " << pixdim[0] << " x " << pixdim[1] << " x "  << pixdim[2] << "\n\n";

	mri.SetVerboseMode(verbose);

	if(outbase.empty()) outbase = BaseName(argv[argc-1]);
	string outpath = Path(outbase.c_str());

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
		sprintf(reg, "%s/bin/flirt -omat %s -in %s -ref %s", getenv("FSLDIR"), (outpath+"in2st.mat").c_str(), argv[argc-1], brainpath.c_str());
		if(verbose) cerr << reg << endl;
		system(reg);
		sprintf(reg, "%s/bin/convert_xfm -inverse -omat %s %s", getenv("FSLDIR"), (outpath+"st2in.mat").c_str(), (outpath+"in2st.mat").c_str());
		if(verbose) cerr << reg << endl;
		system(reg);
		sprintf(reg, "%s/bin/flirt -init %s -applyxfm -o %s -in %s -ref %s", getenv("FSLDIR"), (outpath+"st2in.mat").c_str(), (outpath+"gm.hdr").c_str(), GMpath.c_str(), argv[argc-1]);
		if(verbose) cerr << reg << endl;
		system(reg);
		sprintf(reg, "%s/bin/flirt -init %s -applyxfm -o %s -in %s -ref %s", getenv("FSLDIR"), (outpath+"st2in.mat").c_str(), (outpath+"wm.hdr").c_str(), WMpath.c_str(), argv[argc-1]);
		if(verbose) cerr << reg << endl;
		system(reg);
		sprintf(reg, "%s/bin/flirt -init %s -applyxfm -o %s -in %s -ref %s", getenv("FSLDIR"), (outpath+"st2in.mat").c_str(), (outpath+"csf.hdr").c_str(), CSFpath.c_str(), argv[argc-1]);
		if(verbose) cerr << reg << endl;
		system(reg);

		imagereader.SetFile((outpath+"gm.hdr").c_str());

		pGM = new ZGrayByteImage;
		pWM = new ZGrayByteImage;
		pCSF = new ZGrayByteImage;
		if(!imagereader.ReadFile(*pGM))
		{
			cerr << "Can not read apriori GM probability image gray.tif" << endl;
			delete inputimage;	return -1;
		}
		imagereader.SetFile((outpath+"wm.hdr").c_str());
		if(!imagereader.ReadFile(*pWM))
		{
			cerr << "Can not read apriori WM probability image white.tif" << endl;
			delete inputimage;	return -1;
		}
		imagereader.SetFile((outpath+"csf.hdr").c_str());
		if(!imagereader.ReadFile(*pCSF))
		{
			cerr << "Can not read apriori CSF probability image csf.tif" << endl;
			delete inputimage;	return -1;
		}

		remove((outpath+"in2st.mat").c_str());
		remove((outpath+"st2in.mat").c_str());
		remove((outpath+"gm.hdr").c_str());
		remove((outpath+"wm.hdr").c_str());
		remove((outpath+"csf.hdr").c_str());
		remove((outpath+"gm.img").c_str());
		remove((outpath+"wm.img").c_str());
		remove((outpath+"csf.img").c_str());
	    }
	}
#endif
	
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
		string fname = outbase + "_seg" + ext;
		imagereader.SetFile(fname.c_str());
	
		ret = imagereader.WriteFile(mri.m_Segment);
		if(verbose)
			if(ret) cerr << "Write segmentation image " << fname << " Successfully!" << endl;
			else cerr << "Can not write segmentation image " << fname << endl;
	}

	if(outrestore)
	{
		for(i=0; i<nchannel; i++)
		{
			ZImageBase* pRes = mri.GetRestored(i);
			string fname = outpath + basename[i] + "_restore" + ext;
			imagereader.SetFile(fname.c_str());

			ret = imagereader.WriteFile(*pRes);
			if(verbose)
				if(ret) cerr << "Write restored image " << fname << " of channel " << i << " Successfully!" << endl;
				else cerr << "Can not write restored image " << fname << endl;
		}
	}

	if(outbias)
	{
		for(i=0; i<nchannel; i++)
		{
			ZImageBase* pBias = mri.GetBiasField(i);
			string fname = outpath + basename[i] + "_bias" + ext;
			imagereader.SetFile(fname.c_str());

			ret = imagereader.WriteFile(*pBias);
			if(verbose)
				if(ret) cerr << "Write bias field image " << fname << " of channel " << i << " Successfully!" << endl;
				else cerr << "Can not write bias field image " << fname << endl;
		}
	}

	if(outpbmap)
	{
		for(i=0; i<nclass; i++)
		{
			char tmp[10];
			sprintf(tmp, "_pbmap_%d", i);
			string fname = outbase + tmp + ext;
			imagereader.SetFile(fname.c_str());

			ZImageBase* pProb = mri.GetProbability(i);
			ret = imagereader.WriteFile(*pProb);
			if(verbose)
				if(ret) cerr << "Write probability map image " << fname << " of class " << i << " Successfully!" << endl;
				else cerr << "Can not write probability map image " << fname << endl;
		}
	}

	if(outsegs)
	{
		ZGrayByteImage gboutput(width, height, depth);
		gboutput.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
		
		ZImage<BYTE>::iterator p_out(gboutput), p_seg(mri.m_Segment);

		for(int i=0; i<nclass; i++)
		{
			char tmp[10];
			sprintf(tmp, "_seg_%d", i);
			string fname = outbase + tmp + ext;

			imagereader.SetFile(fname.c_str());

			for(p_out.reset(), p_seg.reset(); p_seg.more(); p_out++, p_seg++)
			{
				if(*p_seg == i+1) *p_out = 1;
				else *p_out = 0;
			}

			ret = imagereader.WriteFile(gboutput);
			if(verbose)
				if(ret) cerr << "Write segmentation image " << fname << " of tissue " << i << " Successfully!" << endl;
				else cerr << "Can not write segmentation image of tissue " << i << ' ' << fname << endl;
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

			ret = imagereader.WriteFile(*(mri.GetPVE(i)));
			if(verbose)
				if(ret) cerr << "Write PVE image of tissue " << i << " Successfully!" << endl;
				else cerr << "Can not write PVE image of tissue " << i << ' ' << fname << endl;
		}
	}

	if(verbose && (nbiter>0 || pve))
	{
		cerr << "\nClass:\t";
		for(int i=0; i<nclass; i++) cerr << "\ttissue " << i;
		cerr << "\ttotal\nVolume:\t";

		vector<float> vol(nclass);

		int c;
		for(c=0; c<nclass; c++)
		{
			if(pve)	for(UINT i=0; i<mri.m_PVE.NRows(); i++) vol[c] += mri.m_PVE(i, c);
			else for(UINT i=0; i<mri.m_post.NRows(); i++) vol[c] += mri.m_post(i, c);
		}

		float total=0;
		for(c=0; c<nclass; c++) 
		{
			vol[c] *= pixdim[0] * pixdim[1] * pixdim[2];
			total+=vol[c];
			printf("\t%-10.1f", vol[c]);
		}
		printf("\t%-10.1f\n\n", total);
	}

	for(i=0; i<nchannel; i++) delete vimages[i];

	return 0;
}
