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
#include <vector>
#include <fstream>
#include <iomanip>   // setprecision(),...
#include <strstream>

#include "imalgorithm.h"
#include "classification.h"
#include "mriseg.h"
#include "zimage/common/commfunc.cpp"
#include "zimage/common/handler.cpp"
#include "zimage/common/mstring.cpp"
#include "zimage/common/operator.cpp"
#include "zimage/common/random.cpp"
#include "zimage/imagebase/image.cpp"
#include "zimage/imagebase/imagebase.cpp"
#include "zimage/imagebase/imagehelper.cpp"
#include "zimage/imagebase/palette.cpp"
#include "zimage/imageio/fileavw.cpp"
#include "zimage/imageio/imageiocommon.cpp"
#include "zimage/imageio/imagereader.cpp"
#include "zimage/imageprocess/classification.cpp"
#include "zimage/imageprocess/filter.cpp"
#include "zimage/imageprocess/imalgorithm.cpp"
#include "zimage/imageprocess/thresholding.cpp"
#include "zimage/imageprocess/vclassification.cpp"
#include "zimage/math/eigensystem.cpp"
#include "zimage/math/jacobi.cpp"
#include "zimage/math/lu.cpp"
#include "zimage/math/mathcommon.cpp"
#include "zimage/math/matrix.cpp"
#include "zimage/math/matrixD.cpp"
#include "zimage/math/svd.cpp"
#include "imageio.h"

using namespace std;

ZMRISegmentation::ZMRISegmentation()
{
	m_bAutoIni = m_bAutoPara = true;
	m_b2c = false;
	m_type = 0;
	m_nVerbose = MAIN;
	m_b3DLowpass = true;
	m_pImage = m_pRestore = 0;

	m_bPVEPara = false;
	PVEB = 1.0f;
}

void ZMRISegmentation::Create(const ZImageBase& image, int type, int nclass, bool b2d, int nbiter,
								int nblowpass, float neighbor, UINT sampleratio, bool autopara, bool scaleseg)
{
	m_pImage = &image;
	image.GetPixelDim(pixdim[0], pixdim[1], pixdim[2]);

	m_type = type;
	m_nbClasses = nclass;
	m_b3D = !b2d && (image.Depth() > 3) && (pixdim[2] > 0);

	m_nbIter = nbiter;
	m_nbLowpass = nblowpass;

	m_nSampleRatio = sampleratio;

	B = neighbor;

	m_bAutoPara = autopara;

	m_bScaleSeg = scaleseg;

	if(m_pImage->ImageWidth() < 100 && m_pImage->ImageWidth()<100)
	{
		m_nSampleRatio = 1;
		B = 0.1f; 
		m_type = 1;
	}
	diaB = 0.8f; ddiaB = 0.7f;
	if(m_b3D) ddiaB *= Max(pixdim[0], pixdim[1]) / pixdim[2];

	if(m_nbClasses==2)		// brain/CSF segmentation
	  //{ m_b2c = true; m_nbClasses = 3; if(m_nbIter>4) m_nbIter=4;} // old buggy code
	  { m_b2c = true; m_nbClasses = 3; }     // new fix SS 14/7/5

	m_mean.resize(m_nbClasses);		m_inimean.resize(m_nbClasses);
	m_stddev.resize(m_nbClasses);	m_inistddev.resize(m_nbClasses); 
	m_variance.resize(m_nbClasses);m_inivariance.resize(m_nbClasses);

	m_prob.resize(m_nbClasses); m_post.resize(m_nbClasses); m_prior.resize(m_nbClasses);
}

void ZMRISegmentation::KMeansVolume()
{
	vector<float> mean(4), stddev(4), prior(4);
	TreeKMean(m_Mri, m_Segment, 3, true, &(*mean.begin()), &(*stddev.begin()), &(*prior.begin()), m_nVerbose>=PARA?&cerr:0);
	// We expect GM and WM to be the two largest classes and to be next to each other
	// We also expect GM and WM to have a similar size otherwise it means the classification is failed.
	if(m_type == 1) // T2
	{
		if(prior[0] > prior[1]*2 || prior[1] > prior[0]*2) 
		{
			if(m_nVerbose >= SLICE) cerr << "3-class K-means failed. Apply 4-class K-means ...\n";
			TreeKMean(m_Mri, m_Segment, 4, true, &(*mean.begin()), &(*stddev.begin()), &(*prior.begin()), m_nVerbose>=PARA?&cerr:0);
			if(prior[2] > prior[0]) // 1 and 2 are GM and WM
			{
				m_inimean[1] = mean[1], m_inimean[2] = mean[2];
				m_inistddev[1] = stddev[1], m_inistddev[2] = stddev[2];
			}
			else					// 0 and 1 are GM and WM
			{
				m_inimean[1] = mean[0], m_inimean[2] = mean[1];
				m_inistddev[1] = stddev[0], m_inistddev[2] = stddev[1];
			}
			m_inimean[0] = mean[3]; m_inistddev[0] = stddev[3];  // CSF is always the brightest
			Classification(m_Mri, m_Segment, 3, true, m_inimean.pbegin(), m_inistddev.pbegin());
		}
		else
		{
			m_inimean[0] = mean[2]; m_inistddev[0] = stddev[2];  // CSF is always the brightest
			m_inimean[1] = mean[0]; m_inistddev[1] = stddev[0];
			m_inimean[2] = mean[1]; m_inistddev[2] = stddev[1];
		}
	}
	else
	{
		if(prior[2] > prior[1]*2 || prior[1] > prior[2]*2) 
		{
			if(m_nVerbose >= SLICE) cerr << "3-class K-means failed. Apply 4-class K-means ...\n";
			TreeKMean(m_Mri, m_Segment, 4, true, &(*mean.begin()), &(*stddev.begin()), &(*prior.begin()), m_nVerbose>=PARA?&cerr:0);
			if(prior[1] > prior[3]) // 1 and 2 are GM and WM
			{
				m_inimean[1] = mean[1], m_inimean[2] = mean[2];
				m_inistddev[1] = stddev[1], m_inistddev[2] = stddev[2];
				m_inimean[0] = mean[0]; m_inistddev[0] = stddev[0];  // CSF
			}
			else					// 2 and 3 are GM and WM
			{
				m_inimean[1] = mean[2], m_inimean[2] = mean[3];
				m_inistddev[1] = stddev[2], m_inistddev[2] = stddev[3];
				if(prior[1] > prior[0]) 
					m_inimean[0] = mean[1], m_inistddev[0] = stddev[1];  // CSF
				else
					m_inimean[0] = mean[0], m_inistddev[0] = stddev[0];  // CSF

			}
			Classification(m_Mri, m_Segment, 3, true, m_inimean.pbegin(), m_inistddev.pbegin());
		}
		else
		{
			m_inimean[0] = mean[0]; m_inistddev[0] = stddev[0];  // CSF is always the brightest
			m_inimean[1] = mean[1]; m_inistddev[1] = stddev[1];
			m_inimean[2] = mean[2]; m_inistddev[2] = stddev[2];
		}
	}
}

void ZMRISegmentation::KMeansSlices()
{
	// including those slices that have a number of pixels no less than
	// half of the maximum number of pixels of a slices for parameter
	// estimation

	int k, max_pixels_a_slice = 0, mslice=m_nDepth/2;
	vector<int> nbpixels(m_nDepth);

	float* ptr=m_Mri.ImageOrigin(); 
	for(k=0; k<m_nDepth; k++, ptr+=m_nSliceSize)
	{
		int m = 0; float* p = ptr;
		for(int i=0; i<m_nSliceSize; i++, p++) if(*p != 0) m ++;

		nbpixels[k] = m;

		if(max_pixels_a_slice < m) max_pixels_a_slice = m, mslice = k;
	}

	for(k=0; k<mslice; k++) if(nbpixels[k] > max_pixels_a_slice/2) break;
	m_nBegin = k;
	for(k=m_nDepth-1; k>mslice; k--) if(nbpixels[k] > max_pixels_a_slice/2) break;
	m_nEnd = k;

	int nslice = 0;
	ZGrayByteImage seg;
	vector<float> mean(3), stddev(3), prior(3);
	if(m_nVerbose>=SLICE) cerr << "Slice-by-slice estiamtion...";
	for(k=m_nBegin; k<=m_nEnd; k++)
	{
		if(m_nVerbose>=SLICE) cerr << "\n\tSlice " << k << "...";
		m_Mri.GotoSlice(k);
		TreeKMean(m_Mri, seg, 3, true, &(*mean.begin()), &(*stddev.begin()), &(*prior.begin()));
		if(m_type == 1) // T2
		{
			if(prior[0] > prior[1]*2 || prior[1] > prior[0]*2) continue;
		}
		else			// T1 and PD
		{
			if(prior[2] > prior[1]*2 || prior[1] > prior[2]*2) continue;
		}
		if(m_nVerbose>=SLICE) cerr << "pass!";
		nslice ++; 
		
		int gm=1, wm=2, csf=0;
		if(m_type == 1) gm=0, wm=1, csf = 2;

		m_inimean[0] += mean[csf], m_inistddev[0] += stddev[csf];
		m_inimean[1] += mean[gm], m_inistddev[1] += stddev[gm];
		m_inimean[2] += mean[wm], m_inistddev[2] += stddev[wm];
	}
	if(nslice) m_inimean /= nslice, m_inistddev /= nslice;
	else
	{
		m_Mri.GotoSlice(mslice);
		KMeansVolume();
	}
	m_Mri.FullROI();
	Classification(m_Mri, m_Segment, 3, true, m_inimean.pbegin(), m_inistddev.pbegin());
}

void ZMRISegmentation::Initialize()
{
	if(m_apriori.size())
	{
		if(m_nVerbose!=QUIET) cerr << "Initial estimation using apriori probability maps...\n";
		ZImage<float>::iterator p_mri(m_Mri);
		int c, pos; float n[3]; n[0] = n[1] = n[2] = 0;
		for(pos=0; p_mri.more(); p_mri++, pos++) 
		{
			if(*p_mri)
			for(c=0; c<3; c++) 
			{
				n[c] += m_apriori[c][pos];
				m_inimean[c] += m_apriori[c][pos] * *p_mri;
			}
		}
		for(c=0; c<3; c++) m_inimean[c] /= n[c];
		for(pos=0, p_mri.reset(); p_mri.more(); p_mri++, pos++) 
		{
			if(*p_mri)
			for(c=0; c<3; c++)
			{
				float v = *p_mri - m_inimean[c];
				m_inistddev[c] += m_apriori[c][pos] * v * v;
			}
		}
		for(c=0; c<3; c++) m_inistddev[c] = float(sqrt(m_inistddev[c]/n[c]));

		m_Segment.Create(m_nWidth, m_nHeight, m_nDepth);
		ZImage<BYTE>::iterator p_seg(m_Segment);
		for(pos=0, p_mri.reset(); p_mri.more(); p_mri++, pos++, p_seg++) 
		{
			if(!*p_mri) continue;

			int clas = 0; float max = 0;
			for(c=0; c<3; c++) 
				if(max < m_apriori[c][pos])	max = m_apriori[c][pos], clas = c;

			*p_seg = clas;
		}
	}
	else
	{
		if(m_bAutoIni)
		{
			if(m_nVerbose!=QUIET) cerr << "Initial K-means segmentation...\n";
			if(m_nbClasses == 3)
			{
				// 0 - CSF; 1, 2 -- GM/WM
				if(m_b3D || m_b2DImage) KMeansVolume();
				else KMeansSlices();
				m_Segment -= 1;
			}
			else
			{
				vector<float> prior(m_nbClasses);
				TreeKMean(m_Mri, m_Segment, m_nbClasses, true, &(*m_inimean.begin()), &(*m_inistddev.begin()), &(*prior.begin()), m_nVerbose>=PARA?&cerr:0);
				if(m_type==1)		// T2
				{
					m_inimean = m_inimean.cshift(1); 
					m_inistddev = m_inistddev.cshift(1);
					ZGrayByteImage::iterator p_seg(m_Segment);
					for(; p_seg.more(); p_seg++) if(*p_seg == m_nbClasses) *p_seg = 0;
				}
				else m_Segment -= 1;
			}
		}
		else
		{
			for(int c=0; c<m_nbClasses; c++) m_inistddev[c] = m_inimean[c] / 5;
			Classification(m_Mri, m_Segment, m_nbClasses, true, m_inimean.pbegin(), m_inistddev.pbegin());
			if(m_type==1)		// T2
			{
				m_inimean = m_inimean.cshift(1); 
				m_inistddev = m_inistddev.cshift(1);
				ZGrayByteImage::iterator p_seg(m_Segment);
				for(; p_seg.more(); p_seg++) if(*p_seg == m_nbClasses) *p_seg = 0;
			}
			else m_Segment -= 1;
		}
	}
	m_mean = m_inimean; m_stddev = m_inistddev;

	ZImage<float>::iterator p_mri(m_Mri);
	for(; p_mri.more(); p_mri++)
	{
		if(*p_mri > 0) 
			*p_mri = float(log(*p_mri));
		else 
			*p_mri = 0;
	}

	if(m_nVerbose>=ITER)
	{
		cerr << "\n    Initial Parameter Estimates:\n";
		for(int c=0; c<m_nbClasses; c++)
			cerr << "\t" << m_mean[c] << "\t" << m_stddev[c] << endl;
		cerr << endl;
	}

	int c;
	for(c=0; c<m_nbClasses; c++)
	{
		m_stddev[c] = m_stddev[c] / m_mean[c];
		m_variance[c] = m_stddev[c] * m_stddev[c];
	}

	for(c=0; c<m_nbClasses; c++) m_mean[c] = float(log(m_mean[c]));
}

void ZMRISegmentation::EMloop(void)
{
	const int WM = m_type==1 ? 1 : m_nbClasses - 1;
	ZGrayFloatImage res;
	ZImage<float>::iterator p_mri(m_Mri);

	ZVector<float> m(m_nbClasses);
	res.Create(m_nWidth, m_nHeight, m_nDepth);

	for(int pos=0; p_mri.more(); p_mri++, pos++)
	{
		if(!*p_mri) continue;

		float norm = 0; int c;
		for(c = 0; c < m_nbClasses; c++)
		{
			float val = *p_mri - m_mean[c];
			val = 0.5f * val * val / m_variance[c];
			val = float(exp(-val) / m_stddev[c] * M_1_SQRT_2PI);
			if(m_bapriori==2 && m_apriori.size() && m_nbClasses == 3) // change to make apriori affect posteriors
			  val *= m_apriori[c][pos];
			norm +=	m_post[c][pos] = val;
		}

		if(norm < 1e-10) norm = 1;
		for(c = 0; c < m_nbClasses; c++) m_post[c][pos] /= norm;
	}

	int iter;
	for (iter = 0; iter < m_nbIter; iter++)
	{
		int pos;

		if(m_nVerbose>=ITER) cerr << "Main iteration, No. " << iter << endl;
		
		memset(m_prob[1], 0, m_nSize * sizeof(float));
		memset(m_prob[2], 0, m_nSize * sizeof(float));

		// MAP estimator of the p_bias field 

		ZImage<float>::iterator p_mri(m_Mri);
		float* p_resmean = m_prob[1];
		float* p_meaninvcov = m_prob[2];

		for(pos=0; pos<m_nSize; p_meaninvcov++, p_resmean++, pos++, p_mri++)
		{
			if(!*p_mri) continue;

			// compute the residual mean   and mean inverse variance 
			for (int c = 1; c < m_nbClasses; c++)
			{
				float tempf = m_post[c][pos] / m_variance[c];
				*p_meaninvcov += tempf;
				*p_resmean += tempf * (*p_mri - m_mean[c]);
			} 
		}

		// filter residual mean

		if(m_nVerbose == ITER) cerr << "    Lowpass filtering...";
		Lowpass(m_prob[1]);
		if(m_nVerbose>=PARA) 
			cerr << "    Mean Residual lowpass filter done!" << endl;
		Lowpass(m_prob[2]);
		if(m_nVerbose>=PARA) 
			cerr << "    Mean inverse covariance lowpass filter done!" << endl;
		if(m_nVerbose == ITER) cerr << "done!\n";

		// computation of the p_bias field
		ZImage<float>::iterator p_bias(m_Bias);
		p_resmean = m_prob[1];
		p_meaninvcov = m_prob[2];
		ZGrayFloatImage::iterator p_res(res);
		for(p_mri.reset(); p_mri.more(); p_mri++, p_res++, p_bias++, p_meaninvcov++, p_resmean++)
		{
			if(!*p_mri) continue;

			// deal with unfortunate division by 0
			if ( Abs(*p_meaninvcov) < TINY ) *p_bias = 0;
			else *p_bias = *p_resmean / *p_meaninvcov; 
			*p_res = *p_mri - *p_bias;
		}

		// segmentation and parameter estimation

		if(m_nVerbose >= PARA && m_nVerbose < PROGRESS) cerr << "    HMRF-EM iteration...\n";

		for(p_mri.reset(), pos=0; p_mri.more(); p_mri++, pos++)
		{
			if(!*p_mri) continue;

			for(int c = 0; c < m_nbClasses; c++)
			{
				if(m_type==1)		// T2
				{
					if(c==1) if(*p_mri > m_mean[0]) { m_prob[c][pos] = 20; continue; }
					if(c==0) if(*p_mri < m_mean[WM]) { m_prob[c][pos] = 20; continue; }
				}
				else
				{
					if(c==0) if(*p_mri > m_mean[WM]) { m_prob[c][pos] = 20; continue; }
				}

				float val = *p_mri - m_mean[c];
				m_prob[c][pos] = 0.5f * val * val / m_variance[c] + log(m_stddev[c]);
			}
		}

		m = m_mean; 			//backup

		float BB = 0.1f;
		if(m_nVerbose == ITER) cerr << "    HMRF-EM Classification...";
		ZProgressIndicator progress("HMRF-EM", 40, 20, m_nVerbose==PROGRESS);
		for(int mrfiter=0; mrfiter<20; mrfiter++)
		{
			int c, change = ICM(res, BB);
			if(BB < B) BB += 0.03f;

			ZGrayFloatImage::iterator p_mri(res);
			for(pos=0; p_mri.more(); p_mri++, pos++)
			{
				if(!*p_mri) continue;

				float norm = 0;

				if(m_bapriori==2 && m_apriori.size() && m_nbClasses == 3)
				  for(c=0; c<m_nbClasses; c++) norm += m_post[c][pos] = float(exp(-m_post[c][pos]) * m_apriori[c][pos]);
				else
				  for(c=0; c<m_nbClasses; c++) norm += m_post[c][pos] = float(exp(-m_post[c][pos]));

				if(norm < 1e-10) { *p_mri = 0;  continue; }
				for(c=0; c<m_nbClasses; c++) m_post[c][pos] /= norm;
			}

			if(change < m_nSize / (10000*BB)) break;

			if(m_bAutoPara)
			{
				m_mean = 0; 
				ZVector<float> nor(m_nbClasses);
				int start = m_nBegin * m_nSliceSize, end = (m_nEnd+1) * m_nSliceSize;

				p_mri.reset(), p_mri += start;
				for(pos=start; pos<end; p_mri+=m_nSampleRatio, pos+=m_nSampleRatio)
				{
					if(!*p_mri) continue;
					for(c=0; c<m_nbClasses; c++)
					{
						float ppp = m_post[c][pos];
						m_mean[c] += ppp * *p_mri;
						nor[c] += ppp;
					}
				}

				for(c=0; c<m_nbClasses; c++) 
				{
					if(nor[c] != 0) m_mean[c] /= nor[c];
					else { m_mean = m; break; }
				}
				if(c<m_nbClasses) break;

				m_variance = 0; 
				p_mri.reset(), p_mri += start;
				for(pos=start; pos<end; p_mri+=m_nSampleRatio, pos+=m_nSampleRatio)
				{
					if(!*p_mri) continue;

					for(c=0; c<m_nbClasses; c++)
					{
						float var = *p_mri - m_mean[c];
						m_variance[c] += m_post[c][pos] * var * var;
					}
				}

				for(c=0; c<m_nbClasses; c++) 
					m_variance[c] /= nor[c], m_stddev[c] = float(sqrt(m_variance[c]));

				for(pos=0, p_mri.reset(); pos<m_nSize; p_mri++, pos++)
				{
					if(!*p_mri) continue;

					for (c = 0; c < m_nbClasses; c++)
					{
						if(m_type==1)		// T2
						{
							if(c==1) if(*p_mri > m_mean[0]) { m_prob[c][pos] = 20; continue; }
							if(c==0) if(*p_mri < m_mean[WM]) { m_prob[c][pos] = 20; continue; }
						}
						else
						{
							if(c==0) if(*p_mri > m_mean[WM]) { m_prob[c][pos] = 20; continue; }
						}
						float val = *p_mri - m_mean[c];
						m_prob[c][pos] = float(val * val / m_variance[c] / 2 + log(m_stddev[c]));
					}
				}
			}
			progress.Step(mrfiter+1);
		}
		if(m_nVerbose == ITER) cerr << "done!\n";
		if(m_nVerbose == PROGRESS) cerr << endl;

		if(m_nVerbose >= PARA) 
		{
			cerr << "    estimated parameters:\n";
			for(int c=0; c<m_nbClasses; c++)
				cerr << "\t" << exp(m_mean[c]) << "\t" << m_stddev[c]*exp(m_mean[c]) << endl;
			cerr << endl;
		}

		if(iter % 4 == 0 && iter != 0)
		{
			for(p_mri.reset(), p_bias.reset(); p_mri.more(); p_mri++, p_bias++)
				if(*p_mri) *p_mri -= *p_bias;

			m_FinalBias += m_Bias;
		}

		// if(m_b2c && m_type == 0) if(m_stddev[1] > m_stddev[2]) break; // old buggy code, taken out SS 14/7/5
	}

	if(iter == 1 || (iter % 4) != 1) 
	{
		ZImage<float>::iterator p_bias(m_Bias);
		for(p_mri.reset(); p_mri.more(); p_mri++, p_bias++)
			if(*p_mri) *p_mri -= *p_bias;

		m_FinalBias += m_Bias;
	}
}

bool ZMRISegmentation::Segment(bool pve, const ZGrayByteImage* gm, const ZGrayByteImage* wm, const ZGrayByteImage* csf, int bapriori)
{
        m_bapriori=bapriori;

	if(m_pImage == 0)
	{
		ZError("ZMRISegmentation::Segment", "No input image was given!");
		return false;
	}

#ifndef _WINDOWS
	if(*PBYTE(m_pImage->GetBuffer()) != 0)
		cerr << "Warning - it looks like BET has not already been used to remove the non-brain parts of the image!\n";
#endif

	if(m_nbClasses > 6)
	{
		ZError("ZMRISegmentation", "The maximum number of classes is 6!");
		return false;
	}

	int width = m_pImage->Width();			// the original width of the image
	int height = m_pImage->Height();		// the original height of the image
	int depth = m_pImage->Depth();			// the original depth of the image
	UINT size = m_pImage->Size();			// the original size of the image
	m_b2DImage = depth == 1;					// if 2D image

	int _left=m_pImage->Left(), _right=m_pImage->Right();
	int _top=m_pImage->Top(), _bottom=m_pImage->Bottom();
	int _front=m_pImage->Front(), _back=m_pImage->Back();

	int left=0, right=width-1, top=0, bottom=height-1, front=0, back=depth-1;
	AutoRange(*m_pImage, left, right, top, bottom, front, back);
	
	if(left == right || top == bottom)
	{
		ZError("ZMRISegmentation::Segment", "Empty input image. Nothing to do!");
		return false;
	}

	ImageRect roi = m_pImage->GetROI();		// backup the original ROI
	m_pImage->SetROI(left, right, top, bottom, front, back);
	
	m_pImage->SendPixelsTo(m_Mri);

	m_nWidth = m_pImage->Width();		// width of the new image with zeor boundaries removed
	m_nHeight = m_pImage->Height();		// height of the new image with zeor boundaries removed
	m_nSliceSize = m_nWidth * m_nHeight;// slicesize of the new image with zeor boundaries removed
	m_nDepth = m_pImage->Depth();		// depth of the new image with zeor boundaries removed
	m_nSize = m_nDepth * m_nSliceSize;	// size of the new image with zeor boundaries removed
	
	offset3d = - m_nSliceSize - m_nWidth - 1;
	offset2d = - m_nWidth - 1;

	if(depth>1) { if(DetermineValidSlices()==false) return false; }
	else m_nBegin = m_nEnd = 0;

	bool inflate = m_pImage->Size() != size;// Black boundaries are removed

	if(m_b3D && (back-front+1 < 10 || pixdim[2] > 6)) 
	{
		cerr << "Slice is too thick, 2D segmentation will be performed!\n";
		m_b3D = false;
	}

	// to remove dark pixels outside brain
	if(m_type == 1)
	{
		ZGrayByteImage seg;
		vector<float> mean(5), stddev(5), prior(5);
		TreeKMean(m_Mri, seg, 5, true, &(*mean.begin()), &(*stddev.begin()), &(*prior.begin()));
		if(prior[0] < 0.1) { seg -= 1; 	m_Mri.Mask(seg); }
	}
	else
	{
		ZImage<float>::iterator p_mri(m_Mri);
		for(; p_mri.more(); p_mri++) if(*p_mri < 10) *p_mri = 0;
	}

	if(gm && wm && csf && m_b3D && m_nbClasses==3)
	{
		gm->SetROI(left, right, top, bottom, front, back);
		wm->SetROI(left, right, top, bottom, front, back);
		csf->SetROI(left, right, top, bottom, front, back);
		int c, pos;
		m_apriori.resize(3);
		for (c = 0; c < 3; c++) m_apriori[c] = new float[m_nSize];
		ZGrayByteImage::const_iterator p_csf(*csf);
		for(pos=0; p_csf.more(); p_csf++, pos++) m_apriori[0][pos] = float(*p_csf)/255.0f;

		int g=1, w=2; if(m_type == 1) g=2, w=1;
		ZGrayByteImage::const_iterator p_gm(*gm);
		for(pos=0; p_gm.more(); p_gm++, pos++) m_apriori[g][pos] = float(*p_gm)/255.0f;
		ZGrayByteImage::const_iterator p_wm(*wm);
		for(pos=0; p_wm.more(); p_wm++, pos++) m_apriori[w][pos] = float(*p_wm)/255.0f;
	}

	if(pixdim[2] > 3 || !m_b3D) m_b3DLowpass = false;
	
	TB = m_b3D ? 4 + diaB*4 + diaB*2 + ddiaB*8 : 4 + diaB*4;

	Initialize();						// Initial segmentation

	if(!m_b2DImage && !m_b3D) {	if(m_nVerbose > SLICE) m_nVerbose = SLICE; }

	if(m_nbIter) 
	{
		m_Bias.Create(m_nWidth, m_nHeight, m_nDepth);
		m_FinalBias.Create(m_nWidth, m_nHeight, m_nDepth);
	}
	
	int _size = m_nSize;
	if(!m_b3D && !m_b2DImage) _size = m_nSliceSize;

	vector<float*>  _post(m_nbClasses);
	for (int c = 0; c < m_nbClasses; c++)
	{
		m_prob[c] = new float[_size];
		m_post[c] = new float[_size];
		memset(m_prob[c], 0, _size*sizeof(float));
		memset(m_post[c], 0, _size*sizeof(float));
		if(!m_b3D && !m_b2DImage) _post[c] = new float[m_nSize];
	}

	int ns = m_b3D ? 1 : m_nDepth;	// number slices
	if(m_nVerbose>QUIET && m_nVerbose<ITER) 
	{
		if(ns==1) cerr << m_nbIter << " main iterations ...\n";
		else cerr << ns << " slice-by-slice segmentation ...\n";
	}

	for(int ii=0; ii<ns; ii++)
	{
		if(!m_b3D && ns>1)
		{
			if(m_nVerbose >= SLICE) cerr << "\nProcessing slice " << ii << " ...";
			m_Mri.GotoSlice(ii);	m_Segment.GotoSlice(ii);
			m_Bias.GotoSlice(ii);	m_FinalBias.GotoSlice(ii);
			m_nSize = m_nSliceSize; m_nDepth = 1;
			m_nBegin = m_nEnd = 0;
		}

		if(m_nbIter) EMloop();

		if(!m_b3D && !m_b2DImage)
		{
			for(int c=0; c<m_nbClasses; c++) 
			{
				memcpy(_post[c]+ii*m_nSliceSize, m_post[c], m_nSliceSize*sizeof(float));
				memset(m_post[c], 0, m_nSliceSize*sizeof(float));
			}
		}
	}

	if(!m_b3D && ns>1) 
	{ 
		if(m_nVerbose >= SLICE) cerr << endl << endl;
		m_Mri.FullROI();	m_Segment.FullROI();
		m_Bias.FullROI();	m_FinalBias.FullROI();
		m_nSize = m_Mri.Size();  m_nDepth = m_Mri.Depth();
	}
		
	if(!m_b3D && ns>1)
	{
		for(int c=0; c<m_nbClasses; c++) delete [](m_post[c]), m_post[c] = _post[c];
	}
	
	if(pve) PVE();

	if(m_nbIter) 
	{
		ZGrayFloatImage::iterator p_bias(m_FinalBias);
		for (; p_bias.more(); p_bias++) *p_bias = float(exp(- *p_bias));
	}

	ZGrayByteImage::iterator p_seg(m_Segment);
	ZImage<float>::iterator p_mri(m_Mri);
	if(m_b2c)
	{
		p_mri.reset();
		for(int i=0; p_seg.more(); i++, p_seg++, p_mri++) 
		{
			if(*p_mri)
			{
				if(*p_seg >= 1) *p_seg = 2; else *p_seg = 1;
				m_post[1][i] += m_post[2][i];
			}
		}
		if(pve)
		{
			ZMatrix<float>::iterator p_pve = m_PVE.begin();
			for(; p_pve != m_PVE.end(); p_pve++) (*p_pve)[1] += (*p_pve)[2];
		}
	}
	else
	{
		if(m_type==1)
			for(p_mri.reset(); p_seg.more(); p_seg++, p_mri++) { if(*p_mri) if(*p_seg == 0) *p_seg = m_nbClasses; }
		else
			for(p_mri.reset(); p_seg.more(); p_seg++, p_mri++) { if(*p_mri) (*p_seg)++; }
	}

	if(inflate)
	{
		m_Segment.Inflate(left-_left, _right-right, top-_top, _bottom-bottom, front-_front, _back-back, true, 0);
		m_FinalBias.Inflate(left-_left, _right-right, top-_top, _bottom-bottom, front-_front, _back-back, true, 1);

		int xysize = width*height;
		for(int c=0; c<m_nbClasses; c++) 
		{
			_post[c] = new float[xysize*depth];
			memset(_post[c], 0, xysize*depth*sizeof(float));
			float* p1 = _post[c]+(front-_front)*xysize+(top-_top)*width+(left-_left);
			float* d1 = m_post[c];
			for(int k=0; k<=back-front; k++, p1+=xysize)
			{
				float* p2 = p1; 
				for(int j=0; j<=bottom-top; j++, p2+=width, d1+=m_nWidth)
				{
					memcpy(p2, d1, m_nWidth*sizeof(float));
				}
			}
			delete [](m_post[c]), m_post[c] = _post[c];
		}

		if(pve)
		{
			ZMatrix<float> opve = m_PVE;
			m_PVE.resize(xysize*depth, m_nbClasses);
			m_PVE = 0;
			ZVector<float> *p1 = m_PVE.pbegin() + (front-_front)*xysize+(top-_top)*width+(left-_left);
			ZVector<float> *d1 = opve.pbegin();
			for(int k=0; k<=back-front; k++, p1+=xysize)
			{
				ZVector<float>* p2 = p1; 
				for(int j=0; j<=bottom-top; j++, p2+=width, d1+=m_nWidth)
				{
					copy(d1, d1+m_nWidth, p2);
				}
			}
		}
	}

	if(m_bScaleSeg)
	{
		for (p_seg.reset(); p_seg.more(); p_seg++)
			*p_seg = BYTE(255 * (*p_seg) / m_nbClasses);
	}

	if(m_nVerbose>=ITER) 
	{
		cerr << "\nThe final statistics are:\n";
		for(int c=0; c<m_nbClasses; c++)
		{
			m_mean[c] = float(exp(m_mean[c]));
			cerr << "\tTissue " << c << ":\t" << m_mean[c] << '\t' << m_mean[c] * m_stddev[c] << endl;
		}
	}

	if(m_b2c) m_nbClasses = 2;

	m_Mri.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	m_FinalBias.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	m_Segment.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	m_pImage->SetROI(roi);
	return true;
}

void ZMRISegmentation::Lowpass(float* imagebuf)
{
	if(m_b3DLowpass)
	{
        float lowz = Max(pixdim[2] / Min(pixdim[0], pixdim[1]), 1.5f);
        lowz *= lowz;
        ZProgressIndicator progress("Lowpass Filtering (Z)", 40, int(m_nbLowpass/lowz), m_nVerbose==PROGRESS);

		for(int niter = 0 ; niter < m_nbLowpass/lowz; niter++)
		{
			float *p_mri = m_Mri.ImageOrigin(), *p_img = imagebuf;

 			for(int i=0; i<m_nSliceSize; i++, p_img++, p_mri++)
			{
				float *pp_mri = p_mri, *pp_img = p_img;
				
				int k;
				for(k=0; k<m_nDepth && !*pp_mri; k++, pp_mri+=m_nSliceSize, pp_img+=m_nSliceSize);
				if(k==m_nDepth) continue;

				bool m1=false, m2=*pp_mri!=0, m3=m2;
				float v1=0, v2=*pp_img, v3=v2;
				for(; k<m_nDepth-1; k++, pp_mri+=m_nSliceSize, pp_img+=m_nSliceSize)
				{
					if(!*pp_mri)  {m1=false; continue;}

					m2 = *pp_mri!=0; 				v2 = *pp_img;
					m3 = *(pp_mri+m_nSliceSize)!=0;	v3 = *(pp_img+m_nSliceSize);
					
					if(m1)
					{
						if(m3)
							*pp_img = (v1 + v2 + v2 + v3) / 4;
						else 
							*pp_img = (v1 + v2 + v2) / 3;
					}
					else if(m3)
					{
						*pp_img = (v3 + v2 + v2) / 3;
					}
					
					v1 = v2; m1 = m2;
				}
				if(*pp_mri && m1) *pp_img = (v1 + *pp_img + *pp_img) / 3;
			}
			
			progress.Step(niter+1);
		}
		if(m_nVerbose==PROGRESS) cerr << endl;
	}

	ZProgressIndicator progress("Lowpass Filtering (XY)", 40, m_nbLowpass, m_nVerbose==PROGRESS);

	for(int niter = 0 ; niter < m_nbLowpass; niter++)
	{
		float* p_mri = m_Mri.ImageOrigin();
		float* p_img = imagebuf;

		for(int k=0; k<m_nDepth; k++)
		{
			for(int j=0; j<m_nHeight; j++)
			{
				int i;
				for(i=0; i<m_nWidth && !*p_mri; i++, p_mri++, p_img++);
				if(i==m_nWidth) continue;
				
				bool m1=false, m2=(*p_mri!=0), m3=m2;
				float v1=0, v2=*p_img, v3=v2;
				for(;i<m_nWidth-1; i++, p_mri++, p_img++)
				{
					if(!*p_mri) {m1=false; continue;}
					
					m2 = *p_mri!=0;		v2 = *p_img;
					m3 = *(p_mri+1)!=0;	v3 = *(p_img+1);

					if(m1)
					{
						if(m3)
							*p_img = (v1 + v2 + v2 + v3) / 4;
						else 
							*p_img = (v1 + v2 + v2) / 3;
					}
					else if(m3)
					{
						*p_img = (v3 + v2 + v2) / 3;
					}
					v1 = v2;
					m1 = m2;
				}

				if(*p_mri && m1) *p_img = (v1 + *p_img + *p_img) / 3;
				p_mri ++; p_img ++;
			}

			p_img = imagebuf + k * m_nSliceSize;
			p_mri = m_Mri.ImageOrigin() + k * m_nSliceSize;

			for(int i=0; i<m_nWidth; i++, p_mri++, p_img++)
			{
				int j;

				float *pp_mri = p_mri, *pp_img = p_img;
				
				for(j=0; j<m_nHeight && !*pp_mri; j++, pp_mri+=m_nWidth, pp_img+=m_nWidth);
				if(j==m_nHeight) continue;

				bool	m1 = false, m2 = (*pp_mri!=0), m3=m2;
				float	v1 = 0,	 v2 = *pp_img,    v3=v2;
				for(; j<m_nHeight-1; j++, pp_mri+=m_nWidth, pp_img+=m_nWidth)
				{
					if(!*pp_mri) {m1=false; continue;}

					m2 = *pp_mri!=0;			v2 = *pp_img;
					m3 = *(pp_mri+m_nWidth)!=0;	v3 = *(pp_img+m_nWidth);
					
					if(m1)
					{
						if(m3)
							*pp_img = (v1 + v2 + v2 + v3) / 4;
						else 
							*pp_img = (v1 + v2 + v2) / 3;
					}
					else if(m3)
					{
						*pp_img = (v3 + v2 + v2) / 3;
					}
					
					v1 = v2;  m1 = m2;
				}
				if(*pp_mri && m1) *pp_img = (v1 + *pp_img + *pp_img) / 3;
			}
		}
		
		progress.Step(niter+1);
	}

	if(m_nVerbose==PROGRESS) cerr << endl;
}

void ZMRISegmentation::MakeBiasDisplay(ZGrayByteImage& gbBias)
{
	if(gbBias.Size() != m_FinalBias.Size())
		gbBias.Create(m_FinalBias.Width(), m_FinalBias.Height(), m_FinalBias.Depth());
	gbBias.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	FillPattern(gbBias);

	unsigned long	hist[256];

	memset(hist, 0, 256 * sizeof(long));

	int nbpixel = 0;

	ZImage<float>::iterator	p_bias(m_FinalBias);
	for (; p_bias.more(); p_bias++) if(*p_bias != 1) *p_bias = 1.0f / *p_bias;
	float min=0, max=256; m_FinalBias.MinMax(min, max);
	float nor = 255.0f / (max - min);

	for (p_bias.reset(); p_bias.more(); p_bias++)
	{
		if(*p_bias == 1) continue;
		int l = int((*p_bias - min) * nor);
		if(l<0) l = 0; if(l > 255) l = 255;
		hist[l] ++;
		nbpixel ++;
	}

	float low=0, high = 0, count=0;
	int i;
	for(i=0; i<256; i++)
	{
		count += hist[i];
		if(count < float(nbpixel) * 0.1) 
			low = float(i)/255 * (max - min) + min;
		if(count > float(nbpixel) * 0.98)
		{
			high = float(i)/255 * (max - min) + min;
			break;
		}
	}

	nor = 255.0f / (high-low);
	ZImage<BYTE>::iterator	gbias(gbBias);
	for(p_bias.reset(); p_bias.more(); gbias++, p_bias++)
	{
		if(*p_bias == 1) continue;

		if(*p_bias > high) *gbias = 255;
		else if(*p_bias < low) *gbias = 0;
		else *gbias = BYTE((*p_bias-low) * nor);
	}

	for (p_bias.reset(); p_bias.more(); p_bias++) 
		if(*p_bias != 1) *p_bias = 1.0f / *p_bias;
}

int ZMRISegmentation::ICM(const ZGrayFloatImage& mri, float BB)
{
	int pos=0, change=0;

	ZImage<BYTE>::iterator			p_seg(m_Segment);
	ZImage<float>::const_iterator	p_mri(mri);
	for (; p_seg.more(); p_seg++, p_mri++, pos++)
	{
		if(!*p_mri) continue;

		PBYTE	neighbor1 = neighbour, src;

		if(m_b3D)
		{
			src = p_seg + offset3d;
			for(int k=0; k<3; k++, src+=m_nSliceSize)
			{
				PBYTE ptr = src;
				for(int j=0; j<3; j++, ptr+=m_nWidth)
				for(int i=0; i<3; i++, neighbor1++)
				{
					if(ptr+i < m_Segment.GetBuffer() || ptr+i >= m_Segment.GetLimit())
						*neighbor1 = 0;
					else 
						*neighbor1 = ptr[i];
				}
			}
		}
		else 
		{
			PBYTE ptr = p_seg + offset2d;
			for(int j=0; j<3; j++, ptr+=m_nWidth)
			for(int i=0; i<3; i++, neighbor1++)
			{
				if(ptr+i < m_Segment.GetBuffer() || ptr+i >= m_Segment.GetLimit())
					*neighbor1 = 0;
				else *neighbor1 = ptr[i]; 
			}
		}

		int c, clas=0; 
		m_prior = 0;

		PBYTE center = m_b3D ? &neighbour[13] : &neighbour[4];
		
		m_prior[center[-1]] ++;
		m_prior[center[1]] ++;
		m_prior[center[-3]] ++;
		m_prior[center[3]] ++;

		m_prior[center[-4]] += diaB;
		m_prior[center[4]] += diaB;
		m_prior[center[-2]] += diaB;
		m_prior[center[2]] += diaB;

		if(m_b3D)
		{
			m_prior[center[-9]] += diaB;
			m_prior[center[9]] += diaB;

			m_prior[center[-12]] += ddiaB;
			m_prior[center[12]] += ddiaB;
			m_prior[center[-10]] += ddiaB;
			m_prior[center[10]] += ddiaB;
			m_prior[center[-8]] += ddiaB;
			m_prior[center[8]] += ddiaB;
			m_prior[center[-6]] += ddiaB;
			m_prior[center[6]] += ddiaB;
		}

		for(c=0; c<m_nbClasses; c++) 
		{
			if(!m_prior[c]) m_prior[c] = -5;
			m_prior[c] = TB - m_prior[c];
		}

		float min = 1e10;
		for(c=0; c<m_nbClasses; c++)
		{
			float v = BB * m_prior[c];
			m_post[c][pos] = m_prob[c][pos] + v;
		
			if(min > m_post[c][pos])
			{
				min = m_post[c][pos];
				clas = c;
			}
		}

		if(*p_seg != clas) change++, *p_seg = BYTE(clas);
	}
	return change;
}

bool ZMRISegmentation::DetermineValidSlices()
{
	// including those slices that have a number of pixels no less than
	// half of the maximum number of pixels of a slices for parameter
	// estimation

	int k, max_pixels_a_slice = 0;
	vector<int> nbpixels(m_nDepth);

	float* ptr=m_Mri.ImageOrigin(); 
	for(k=0; k<m_nDepth; k++, ptr+=m_nSliceSize)
	{
		int m = 0; float* p = ptr;
		for(int i=0; i<m_nSliceSize; i++, p++) if(*p != 0) m ++;

		nbpixels[k] = m;

		if(max_pixels_a_slice < m) max_pixels_a_slice = m;
	}

	for(k=0; k<m_nDepth; k++) if(nbpixels[k] > max_pixels_a_slice/3) break;
	m_nBegin = k;
	for(k=m_nDepth-1; k>=0; k--) if(nbpixels[k] > max_pixels_a_slice/3) break;
	m_nEnd = k;

	if(m_nBegin > m_nEnd)
	{
		ZError("ZMRISegmentation::DetermineValidSlices", "Invalid image!");
		return false;
	}

	if(m_nVerbose>=SLICE) cerr << "Parameter estimation will be performed from slice " 
						<< m_nBegin << " to " << m_nEnd << endl;

	return true;
}

void ZMRISegmentation::PVE()
{
	static long seed=-1;

	if(m_nVerbose > QUIET) cerr << "Estimating PVE...\n";
	
// free some memory
	UINT c;
	for (c = 0; c < m_prob.size(); c++) delete []m_prob[c]; m_prob.clear(); 
	for (c = 0; c < m_apriori.size(); c++) delete []m_apriori[c]; m_apriori.clear();
	m_Bias.CleanUp();

	InitializePVE();

	ZVector<float> npve(m_nbClasses);

	float total;
	ZProgressIndicator progress("HMRF-EM for PVE", 40, 100, m_nVerbose==PROGRESS);
	for(int mrfiter=0; mrfiter<100; mrfiter++)
	{
		int c;
		total = 0;
		ZImage<float>::iterator p_mri(m_Mri);
		ZMatrix<float>::iterator p_pve = m_PVE.begin();
		ZVector<float>::iterator p_eng = m_energy.pbegin();
		for(; p_mri.more(); p_mri++, p_eng++, p_pve++)
		{
			if(!*p_mri) continue;

			npve = 0;

			int which = rand() % (m_nbClasses-1);

			npve[which] = rand1(seed);
			npve[which+1] = 1.0f - npve[which];

			float cond = ComputeConditional(*p_mri, npve);
			float post = ComputePosterior(m_Mri, &(*p_mri), &(*p_pve), npve);
			float eng = cond + post;

			float change = *p_eng - eng;
			if(change>0) *p_eng = eng, *p_pve = npve;

			total += *p_eng;
		}

		if(m_bPVEPara)
		{
			m_mean = 0; m_variance = 0;
			ZVector<float> nor(m_nbClasses);
			ZMatrix<float>::iterator p_post = m_PVE.begin();
			for(p_mri.reset(); p_mri.more(); p_mri+=m_nSampleRatio, p_post+=m_nSampleRatio)
			{
				if(!*p_mri) continue;

				for(c=0; c<m_nbClasses; c++)
				{
					m_mean[c] += (*p_post)[c] * *p_mri;
					nor[c] += (*p_post)[c];
				}
			}

			for(c=0; c<m_nbClasses; c++) m_mean[c] /= nor[c];

			p_post = m_PVE.begin();
			for(p_mri.reset(); p_mri.more(); p_mri+=m_nSampleRatio, p_post+=m_nSampleRatio)
			{
				if(!*p_mri) continue;

				for(c=0; c<m_nbClasses; c++)
				{
					float var = *p_mri - m_mean[c];
					m_variance[c] += (*p_post)[c] * var * var;
				}
			}

			for(c=0; c<m_nbClasses; c++)
			{
				m_variance[c] /= nor[c];
				
				m_stddev[c] = float(sqrt(m_variance[c]));
			}
		}
		progress.Step(mrfiter+1);
	}
	if(m_nVerbose >= PARA) cerr << "Final energy " << total << endl; 

	if(m_nVerbose >= PARA) 
	{
		cerr << "    estimated parameters:\n";
		for(int c=0; c<m_nbClasses; c++)
			cerr << "\t" << exp(m_mean[c]) << "\t" << m_stddev[c]*exp(m_mean[c]) << endl;
		cerr << endl;
	}
}

void ZMRISegmentation::InitializePVE()
{
	m_PVE.resize(m_nSize, m_nbClasses);
	m_energy.resize(m_nSize);

	ZGrayFloatImage::iterator p_mri(m_Mri);
	ZMatrix<float>::iterator p_pve = m_PVE.begin();

	/* old 3.3 code */
	/*	for(; p_mri.more(); p_mri++, p_pve++)
	{
		if(!*p_mri) continue;

		double	norm = 0; 
		for (int c = 0; c < m_nbClasses; c++)
		{
			float val = *p_mri - m_mean[c];
			val = val * val / m_variance[c] / 2;
			(*p_pve)[c] = val = float(exp(-val) / m_stddev[c] * M_1_SQRT_2PI);
			norm += val;
		}

		if(norm < TINY) *p_mri = 0; else *p_pve /= norm;
		}*/

	/* new 3.4 code */
	int pos=0;
	for(pos=0; p_mri.more();pos++, p_mri++, p_pve++)
	  {
	    if(!*p_mri) continue;
	    for (int c = 0; c < m_nbClasses; c++)
	      (*p_pve)[c]=m_post[c][pos];//JV fix of initialisation of pve histogram
	  }

	float total=0; p_pve=m_PVE.begin();
	ZVector<float>::iterator p_eng = m_energy.pbegin();
	for(p_mri.reset(); p_mri.more(); p_mri++, p_pve++, p_eng++)
	{
		if(!*p_mri) continue;

		float cond = ComputeConditional(*p_mri, *p_pve);
		float post = ComputePosterior(m_Mri, &(*p_mri), &(*p_pve), *p_pve);
		total += *p_eng = cond + post;
	}
	if(m_nVerbose >= PARA) cerr << "initial energy " << total << endl; 
}

float ZMRISegmentation::ComputeConditional(float mri, const ZVector<float>& pve)
{
	int c; float m=0, v=0;
	for(c=0; c<m_nbClasses; c++) 
	{
		m += pve[c] * m_mean[c];
		float temp = pve[c] * pve[c];
		v += temp * m_variance[c];
	}

	float temp = mri - m;
	float cond = (temp * temp / v + log(v)) / 2;

	return cond;
}

float ZMRISegmentation::ComputePosterior(ZGrayFloatImage& mri, float* p_mri, ZVector<float>* p_pve, ZVector<float>& npve)
{

#define _POST(offset) \
	if(p_mri[offset] > 0) \
	{ ZVector<float>& pve = *(p_pve+offset); \
	  for(int c=0; c<m_nbClasses; c++) post += (pve[c] - npve[c]) * (pve[c] - npve[c]); }

	float *pmri, post=0;

	if(m_b3D)
	{
		if(PBYTE(pmri = p_mri - m_nSliceSize) > mri.GetBuffer() && *pmri>0) 
		{
			ZVector<float>& pve = *(p_pve - m_nSliceSize);
			float dif = 0;
			for(int c=0; c<m_nbClasses; c++) dif += (pve[c] - npve[c]) * (pve[c] - npve[c]);
			post += dif * diaB;
		}
		if(PBYTE(pmri = p_mri + m_nSliceSize) < mri.GetLimit() && *pmri>0) 
		{
			ZVector<float>& pve = *(p_pve + m_nSliceSize);
			float dif = 0;
			for(int c=0; c<m_nbClasses; c++) dif += (pve[c] - npve[c]) * (pve[c] - npve[c]);
			post += dif * diaB;
		}
	}

	if(PBYTE(p_mri - 1) > mri.GetBuffer()) _POST(-1); 
	if(PBYTE(p_mri + 1) < mri.GetLimit()) _POST(1); 

	if(PBYTE(p_mri - m_nWidth) > mri.GetBuffer()) _POST(-m_nWidth); 
	if(PBYTE(p_mri + m_nWidth) < mri.GetLimit()) _POST(m_nWidth); 

	if(PBYTE(p_mri-1-m_nWidth) > mri.GetBuffer()) _POST(-1-m_nWidth); 
	if(PBYTE(p_mri+1-m_nWidth) > mri.GetBuffer()) _POST(1-m_nWidth); 

	if(PBYTE(p_mri-1+m_nWidth) < mri.GetLimit()) _POST(-1+m_nWidth); 
	if(PBYTE(p_mri+1+m_nWidth) < mri.GetLimit()) _POST(1+m_nWidth); 

	return post * PVEB;
}

ZGrayFloatImage& ZMRISegmentation::GetPVE(int nclass)
{
	if(nclass < 0) nclass = 0;
	if(nclass >= m_nbClasses) nclass = m_nbClasses - 1;
	m_PVEimage.Create(m_pImage->Width(), m_pImage->Height(), m_pImage->Depth());
	m_PVEimage.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	ZGrayFloatImage::iterator p_ipve(m_PVEimage);
	ZMatrix<float>::const_iterator p_pve = m_PVE.begin();
	for(; p_ipve.more(); p_ipve++, p_pve++) *p_ipve = (*p_pve)[nclass];

	return m_PVEimage;
}
