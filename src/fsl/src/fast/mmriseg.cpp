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

#include "imalgorithm.h"
#include "common.h"
#include "vclassification.h"
#include "mmriseg.h"
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
#include "zmath.h"

using namespace std;

static ZVector<float>	tv, tvc;
static ZMatrix<float>	m1, m2;

ZMultiMRISegmentation::ZMultiMRISegmentation()
{
	m_nVerbose = MAIN;
	m_b3DLowpass = true;
	m_images = 0;
	m_bBias = true;
	m_FloatImage = 0;
}

void ZMultiMRISegmentation::Create(const std::vector<ZImageBase*>& images, int nclass, bool b2d, 
					int nbiter, int nblowpass, float neighbor, UINT sampleratio, bool bBias)
{
	m_images = &images;
	m_nChannel = images.size();

	images[0]->GetPixelDim(pixdim[0], pixdim[1], pixdim[2]);

	m_nbClasses = nclass;
	m_b3D = !b2d && (images[0]->Depth() > 3) && (pixdim[2] > 0);

	m_nbIter = nbiter;
	m_nbLowpass = nblowpass;

	m_nSampleRatio = sampleratio;

	B = neighbor;

	m_bBias = bBias;

	if(images[0]->ImageWidth() < 100)
	{
		m_nSampleRatio = 1;
		m_b3D = false;
		B = 0.1f;
	}
	diaB = 0.8f; ddiaB = 0.7f;
	if(m_b3D) ddiaB *= Max(pixdim[0], pixdim[1]) / pixdim[2];

	m_mean.resize(m_nbClasses, m_nChannel);		
	m_inimean.resize(m_nbClasses, m_nChannel);
	m_stddev.resize(m_nbClasses);
	m_inicovariance.resize(m_nbClasses); 
	int c;
	for(c=0; c<m_nbClasses; c++) m_inicovariance[c].resize(m_nChannel, m_nChannel);
	m_icovariance.resize(m_nbClasses); 
	for(c=0; c<m_nbClasses; c++) m_icovariance[c].resize(m_nChannel, m_nChannel);

	m_prior.resize(m_nbClasses);

	tv.resize(m_nbClasses); tvc.resize(m_nChannel);
}

void ZMultiMRISegmentation::Initialize()
{
	m_mask.resize(m_nSize);
	ZMatrix<float>::row_iterator p_mri = m_Mri.pbegin();
	vector<BYTE>::iterator p_mask = m_mask.begin();
	for(; p_mri != m_Mri.pend(); p_mri++, p_mask++) 
		if((*p_mri)[0]) { *p_mri = log(*p_mri); *p_mask=1; }
		else *p_mask = 0;

	m_Segment.Create(m_nWidth, m_nHeight, m_nDepth);
	if(m_apriori.size())
	{
		if(m_nVerbose!=QUIET) cerr << "Initial estimation using apriori probability maps...\n";
		ZMatrix<float>::row_iterator p_mri = m_Mri.pbegin();
		ZMatrix<BYTE>::row_iterator p_pri = m_apriori.pbegin();
		int c; ZVector<float> n(3);
		for(; p_mri != m_Mri.pend(); p_mri++, p_pri++) 
		{
			if((*p_mri)[0])
			for(c=0; c<3; c++) 
			{
				float ppp = float((*p_pri)[c]) / 255.0f;
				n[c] += ppp;
				m_inimean[c] += *p_mri * ppp;
			}
		}

		for(c=0; c<3; c++) 
		{
			if(n[c] == 0) ZError("Initialization", "Initial estimation has been failed!");
			m_inimean[c] /= n[c];
		}

		ZVector<float> dif(m_nChannel);
		
		for(p_mri=m_Mri.pbegin(), p_pri=m_apriori.pbegin(); p_mri != m_Mri.pend(); p_mri++, p_pri++)
		{
			if((*p_mri)[0])
			for(c=0; c<3; c++)
			{
				dif = *p_mri; dif -= m_inimean[c];
				float ppp = float((*p_pri)[c]) / 255.0f;
				for(int i=0; i<m_nChannel; i++)
				for(int j=0; j<m_nChannel; j++)
					m_inicovariance[c](i,j) += dif[i] * dif[j] * ppp;
			}
		}
		for(c=0; c<3; c++) if(n[c] > 1) m_inicovariance[c] /= n[c];

		ZImage<BYTE>::iterator p_seg(m_Segment);
		for(p_mri=m_Mri.pbegin(), p_pri=m_apriori.pbegin(); p_mri != m_Mri.pend(); p_mri++, p_pri++, p_seg++)
		{
			if(!(*p_mri)[0]) continue;

			int clas = 0; float max = 0;
			for(c=0; c<3; c++) 
				if(max < (*p_pri)[c]) max = (*p_pri)[c], clas = c;

			*p_seg = clas;
		}
	}
	else 
	{
		if(m_nVerbose!=QUIET) cerr << "initial segmentation by KMeans...." << endl;
		ZVector<float> prior(m_nbClasses);
		TreeKMean(m_Mri, m_Segment, m_nbClasses, true, m_inimean, m_inicovariance, prior.pbegin(), m_nVerbose>=PARA?&cerr:0);
		m_Segment -= 1;
	}

	m_mean = m_inimean;  m_covariance = m_inicovariance; 
	for(int c=0; c<m_nbClasses; c++)
	{
		m_icovariance[c] = inv(m_covariance[c]);
		m_stddev[c] = sqrt(Abs(det(m_covariance[c])));
	}

	if(m_nVerbose>=ITER)
	{
		cerr << "\n    Initial Parameter Estimates:\n";
		for(int c=0; c<m_nbClasses; c++)
			cerr << "\t" << exp(m_mean[c]) << endl;
		cerr << endl;
	}
}

void ZMultiMRISegmentation::EMloop(void)
{
	int iter;

	ZMatrix<float>			meanres;
	vector< ZMatrix<float> >		meaninvcov;
	ZMatrix<float>			res;
	ZMatrix<float>			mm(m_nChannel, m_nChannel);
	ZVector<float>			tvc(m_nChannel);
	ZVector<float>			tv(m_nbClasses);

	ZMatrix<float>::row_iterator p_mri = m_Mri.pbegin(), p_post=m_post.pbegin();
        ZMatrix<BYTE>::row_iterator p_pri = m_apriori.pbegin();
	for(; p_mri != m_Mri.pend(); p_mri++, p_post++)
	{
		if(!(*p_mri)[0]) continue;

		float norm = 0; int c;
		for(c = 0; c < m_nbClasses; c++)
		{
			tvc = *p_mri; tvc -= m_mean[c];
			float v = 0;
			for(int i=0; i<m_nChannel; i++)
			{
				float u = 0;
				for(int j=0; j<m_nChannel; j++) u += tvc[j] * m_icovariance[c](j, i);
				v += u * tvc[i];
			}

			v = float(exp(-v/2) / m_stddev[c] * M_1_SQRT_2PI);
                        if(m_bapriori==2 && m_apriori.size() && m_nbClasses == 3)
			  v *= (*p_pri)[c];
                        norm += (*p_post)[c] = v;
		}

		if(norm < 1e-10) norm = 1;
		for(c = 0; c < m_nbClasses; c++) (*p_post)[c] /= norm;
	}

	if(m_bBias) 
	{
		m_Bias.resize(m_nSize, m_nChannel);
		meanres.resize(m_nSize, m_nChannel);
		meaninvcov.resize(m_nSize);
		for(int i=0; i<m_nSize; i++) meaninvcov[i].resize(m_nChannel, m_nChannel);
		res.resize(m_nSize, m_nChannel);
	}

	for (iter = 0; iter < m_nbIter; iter++)
	{
		if(m_nVerbose>=ITER) cerr << "Main iteration, No. " << iter << endl;

		///// M-Step : MAP estimator of the p_bias field 

		p_mri = m_Mri.pbegin(), p_post=m_post.pbegin();
		vector<BYTE>::iterator p_mask = m_mask.begin();
		if(m_bBias)
		{
			meanres = 0;
			fill(meaninvcov.begin(), meaninvcov.end(), 0);

			ZMatrix<float>::row_iterator p_mres = meanres.pbegin();
			ZMatrix<float>* p_minv = &(*meaninvcov.begin());
			for(; p_mri != m_Mri.pend(); p_post++, p_mres++, p_minv++, p_mri++, p_mask++)
			{
				if(!*p_mask) continue;

				for (int c = 0; c < m_nbClasses; c++)
				{
					mm = m_icovariance[c]; mm *= (*p_post)[c];
					*p_minv += mm;

					tvc = *p_mri; tvc -= m_mean[c];

					for(int i=0; i<m_nChannel; i++)
					for(int j=0; j<m_nChannel; j++)
						(*p_mres)[i] += mm(i, j) * tvc[j];
				} 
			} 
			if(m_nVerbose == ITER) cerr << "    Lowpass filtering...";
			Lowpass(meanres, meaninvcov);
			if(m_nVerbose == ITER) cerr << "done!\n";

			// computation of the p_bias field
			p_mask = m_mask.begin();
			p_mres = meanres.pbegin();
			p_minv = &(*meaninvcov.begin());
			ZMatrix<float>::row_iterator p_bias = m_Bias.pbegin();
			for(; p_bias!=m_Bias.pend(); p_mask++, p_mres++, p_minv++, p_bias++)
			{
				if(!*p_mask) continue;

				float v = 0;
				for(int i=0; i<m_nChannel; i++) v += (*p_minv)(i,i);
				if(Abs(v) < 0.00001f) *p_bias = 0;
				else
				{
					//*p_bias = p_minv->Inverse() * *p_mres;

					mm = inv(*p_minv);
					for(int i=0; i<m_nChannel; i++)
					{
						float v=0;
						for(int j=0; j<m_nChannel; j++) v += mm(i, j) * (*p_mres)[j];
						(*p_bias)[i] = v;
					}
				}
			}

			p_bias = m_Bias.pbegin(); p_mri = m_Mri.pbegin();
			p_mask = m_mask.begin();
			ZMatrix<float>::row_iterator p_res = res.pbegin();
			for(; p_mri!=m_Mri.pend(); p_mri++, p_bias++, p_res++, p_mask++)
			{
				if(!*p_mask) continue;

				*p_res = *p_mri; *p_res -= *p_bias;
			}
		}

		p_mri = m_bBias ? res.pbegin() : m_Mri.pbegin();
		ZMatrix<float>::row_iterator p_prob = m_prob.pbegin();
		p_mask = m_mask.begin();
		for(; p_prob!=m_prob.pend(); p_mri++, p_prob++, p_mask++)
		{
			if(!*p_mask) continue;

			for (int c = 0; c < m_nbClasses; c++)
			{
				double v = 0; 
				tvc = *p_mri; tvc -= m_mean[c];

				for(int i=0; i<m_nChannel; i++)
				{
					float u = 0;
					for(int j=0; j<m_nChannel; j++) u += tvc[j] * m_icovariance[c](j, i);
					v += u * tvc[i];
				}

				(*p_prob)[c] = float(v / 2 + log(m_stddev[c]));
			}
		}

		float BB = 0.1f;
		ZProgressIndicator progress("HMRF-EM", 40, 20, m_nVerbose==PROGRESS);
		for(int mrfiter=0; mrfiter<20; mrfiter++)
		{
			int c, change = ICM(BB);
			if(BB < B) BB += 0.03f;

			ZMatrix<float>::row_iterator p_post=m_post.pbegin();
			p_mask = m_mask.begin();

                        ZMatrix<BYTE>::row_iterator p_pri=m_apriori.pbegin();
                        for(; p_post!=m_post.pend(); p_post++, p_mask++, p_pri++)
			{
				if(!*p_mask) continue;

				float norm = 0;

				if(m_bapriori==2 && m_apriori.size() && m_nbClasses == 3)
				  for(c=0; c<m_nbClasses; c++) norm += (*p_post)[c] = float(exp(-(*p_post)[c]) * (*p_pri)[c]);
				else
				  for(c=0; c<m_nbClasses; c++) norm += (*p_post)[c] = float(exp(-(*p_post)[c]));

				if(norm < 1e-10) { *p_mask = 0;  continue; }
				*p_post /= norm;
			}

			if(change < m_nSize / (10000*BB)) break;

			m_mean = 0; 
			fill(m_covariance.begin(), m_covariance.end(), 0);
			tv = 0;
			int start = m_nBegin * m_nSliceSize;
			ZMatrix<float>::row_iterator end = m_post.pbegin() + (m_nEnd+1) * m_nSliceSize;
			
			p_mri = m_bBias ? res.pbegin() : m_Mri.pbegin(); p_mri += start;
			p_post = m_post.pbegin() + start;
			p_mask = m_mask.begin() + start;
			for(; p_post<end; p_mri+=m_nSampleRatio, p_post+=m_nSampleRatio, p_mask+=m_nSampleRatio)
			{
				if(!*p_mask) continue;

				for(c=0; c<m_nbClasses; c++)
				{
					float ppp = (*p_post)[c];
					for(int i=0; i<m_nChannel; i++) m_mean[c][i] += (*p_mri)[i] * ppp;
					tv[c] += ppp;
				}
			}

			for(c=0; c<m_nbClasses; c++) m_mean[c] /= tv[c];

			p_mri = m_bBias ? res.pbegin() : m_Mri.pbegin(); p_mri += start;
			p_post = m_post.pbegin() + start;
			p_mask = m_mask.begin() + start;
			for(; p_post<end; p_mri+=m_nSampleRatio, p_post+=m_nSampleRatio, p_mask+=m_nSampleRatio)
			{
				if(!*p_mask) continue;

				for(c=0; c<m_nbClasses; c++)
				{
					tvc = *p_mri; tvc -= m_mean[c];

					float ppp = (*p_post)[c];
					for(int i=0; i<m_nChannel; i++)
					for(int j=0; j<m_nChannel; j++)
						m_covariance[c](i,j) += tvc[i] * tvc[j] * ppp;
				}
			}
				
			for(c=0; c<m_nbClasses; c++)
			{
				m_covariance[c] /= tv[c];
				m_icovariance[c] = inv(m_covariance[c]);
				m_stddev[c] = sqrt(Abs(det(m_covariance[c])));
			}

			p_mri = m_bBias ? res.pbegin() : m_Mri.pbegin(); 
			p_prob = m_prob.pbegin();
			p_mask = m_mask.begin();
			for(; p_prob<m_prob.pend(); p_mri++, p_prob++, p_mask++)
			{
				if(!*p_mask) continue;

				for (c = 0; c < m_nbClasses; c++)
				{
					tvc = *p_mri; tvc -= m_mean[c];
					float v = 0;
					for(int i=0; i<m_nChannel; i++)
					{
						float u = 0;
						for(int j=0; j<m_nChannel; j++) u += tvc[j] * m_icovariance[c](j, i);
						v += u * tvc[i];
					}

					(*p_prob)[c] = float(v/2 + log(m_stddev[c]));
				}
			}
			progress.Step(mrfiter+1);
		}

		if(m_nVerbose==PROGRESS) cerr << endl;

		if(m_nVerbose >= PARA) 
		{
			cerr << "    estimated parameters:\n";
			for(int c=0; c<m_nbClasses; c++)
				cerr << "\t" << exp(m_mean[c]) << endl;
			cerr << endl;
		}
	}
}

bool ZMultiMRISegmentation::Segment(bool pve, const ZGrayByteImage* gm, const ZGrayByteImage* wm, const ZGrayByteImage* csf, int bapriori)
{
        m_bapriori=bapriori;

	if(m_images == 0)
	{
		ZError("ZMultiMRISegmentation::Segment", "No input image was given!");
		return false;
	}

	if(m_nbClasses > 6)
	{
		ZError("ZMultiMRISegmentation", "The maximum number of classes is 6!");
		return false;
	}

	int width = (*m_images)[0]->Width();			// the original width of the image
	int height = (*m_images)[0]->Height();		// the original height of the image
	int depth = (*m_images)[0]->Depth();			// the original depth of the image
	UINT size = (*m_images)[0]->Size();			// the original size of the image
	m_b2DImage = depth == 1;					// if 2D image

	int _left=(*m_images)[0]->Left(), _right=(*m_images)[0]->Right();
	int _top=(*m_images)[0]->Top(), _bottom=(*m_images)[0]->Bottom();
	int _front=(*m_images)[0]->Front(), _back=(*m_images)[0]->Back();

	m_left=0, m_right=width-1, m_top=0, m_bottom=height-1, m_front=0, m_back=depth-1;
	if(!AutoRange(*m_images, m_left, m_right, m_top, m_bottom, m_front, m_back)) return false;
	
	if(m_left == m_right || m_top == m_bottom)
	{
		ZError("ZMultiMRISegmentation::Segment", "Empty input image. Nothing to do!");
		return false;
	}

	int i;
	ImageRect roi = (*m_images)[0]->GetROI();		// backup the original ROI
	for(i=0; i<m_nChannel; i++) (*m_images)[i]->SetROI(m_left, m_right, m_top, m_bottom, m_front, m_back);
	Reform(*m_images, m_Mri);

	if(m_b2DImage) m_front = m_back = _front = _back = 0;
	
	m_nWidth = (*m_images)[0]->Width();		// width of the new image with zeor boundaries removed
	m_nHeight = (*m_images)[0]->Height();		// height of the new image with zeor boundaries removed
	m_nSliceSize = m_nWidth * m_nHeight;// slicesize of the new image with zeor boundaries removed
	m_nDepth = (*m_images)[0]->Depth();		// depth of the new image with zeor boundaries removed
	m_nSize = m_nDepth * m_nSliceSize;	// size of the new image with zeor boundaries removed
	
	m_prob.resize(m_nSize, m_nbClasses); 
	m_post.resize(m_nSize, m_nbClasses); 

	offset3d = - m_nSliceSize - m_nWidth - 1;
	offset2d = - m_nWidth - 1;

	if(depth>1) { if(DetermineValidSlices()==false) return false; }
	else m_nBegin = m_nEnd = 0;

	ZMatrix<float>::row_iterator p_mri = m_Mri.pbegin() + 1;
	for(; p_mri < m_Mri.pend()-1; p_mri++) 
	{
		for(i=0; i<m_nChannel; i++)
		{
			if((*p_mri)[i] < 10) {*p_mri = 0; break; }

			if((*(p_mri-1))[i]==0 && (*(p_mri+1))[i]==0) { *p_mri = 0; break; }
		}
	}

	bool inflate = m_nSize != int(size);// Black boundaries are removed

	if(m_b3D && (m_back-m_front+1 < 10 || pixdim[2] > 6)) 
	{
		cerr << "Slice is too thick, 2D segmentation will be performed!\n";
		m_b3D = false;
	}

	if(gm && wm && csf && m_b3D && m_nbClasses==3)
	{
		gm->SetROI(m_left, m_right, m_top, m_bottom, m_front, m_back);
		wm->SetROI(m_left, m_right, m_top, m_bottom, m_front, m_back);
		csf->SetROI(m_left, m_right, m_top, m_bottom, m_front, m_back);
		m_apriori.resize(m_nSize, 3);
		ZMatrix<BYTE>::row_iterator p_pri = m_apriori.pbegin();
		ZGrayByteImage::const_iterator p_csf(*csf);
		ZGrayByteImage::const_iterator p_gm(*gm);
		ZGrayByteImage::const_iterator p_wm(*wm);
		for(; p_csf.more(); p_csf++, p_gm++, p_wm++, p_pri++) 
		{
			(*p_pri)[0] = *p_csf;
			(*p_pri)[1] = *p_gm;
			(*p_pri)[2] = *p_wm;
		}
	}

	if(pixdim[2] > 3 || !m_b3D) m_b3DLowpass = false;
	
	TB = m_b3D ? 4 + diaB*4 + diaB*2 + ddiaB*8 : 4 + diaB*4;

	Initialize();						// Initial segmentation

	m_prob.resize(m_nSize, m_nbClasses);
	m_post.resize(m_nSize, m_nbClasses);

	if(m_nVerbose==MAIN) cerr << m_nbIter << " main iterations ...\n";

	EMloop();

	ZGrayByteImage::iterator p_seg(m_Segment);
	for(p_mri=m_Mri.pbegin(); p_seg.more(); p_seg++, p_mri++) { if((*p_mri)[0]) (*p_seg)++; }

	if(m_bBias)
	{
		ZMatrix<float>::row_iterator p_bias = m_Bias.pbegin();
		for(; p_bias!=m_Bias.pend(); p_bias++) 
		for(i=0; i<m_nChannel; i++)
		{
			if((*p_bias)[i] < -0.7) { *p_bias = 0; break; }
		}

		m_Mri -= m_Bias;
	}
	if(pve) PVE();

	if(inflate)
	{
		m_Segment.Inflate(m_left-_left, _right-m_right, m_top-_top, _bottom-m_bottom, m_front-_front, _back-m_back, true, 0);
	}

	if(m_nVerbose>=ITER) 
	{
		cerr << "\nThe final statistics are:\n";
		for(int c=0; c<m_nbClasses; c++)
		{
			m_mean[c] = exp(m_mean[c]);
			cerr << "\tTissue " << c << ":\t" << m_mean[c] << endl;
		}
	}

	for(i=0; i<m_nChannel; i++) (*m_images)[i]->SetROI(roi);
	m_Segment.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	return true;
}

void ZMultiMRISegmentation::Lowpass(ZMatrix<float>& meanres, vector< ZMatrix<float> >& meaninvcov)
{
	ZProgressIndicator progress("Lowpass Filtering (XY)", 40, m_nbLowpass, m_nVerbose==PROGRESS);

	vector< vector<float> > vec(m_nChannel + m_nChannel*m_nChannel);

	int i=0;
	for(; i<m_nChannel; i++)
	{
		vec[i].resize(m_nSize);

		for(int k=0; k<m_nSize; k++) vec[i][k] = meanres[k][i];
	}

	for(; i<m_nChannel + m_nChannel*m_nChannel; i++)
	{
		vec[i].resize(m_nSize);

		int m = (i - m_nChannel) / m_nChannel;
		int n = (i - m_nChannel) % m_nChannel;

		for(int k=0; k<m_nSize; k++) vec[i][k] = meaninvcov[k](m, n);
	}

	for(int niter = 0 ; niter < m_nbLowpass; niter++)
	{
		for(int h=0; h<m_nChannel + m_nChannel*m_nChannel; h++)
		{
			vector<float>::iterator p_img = vec[h].begin();
			vector<BYTE>::iterator p_mask = m_mask.begin();

			for(int k=0; k<m_nDepth; k++)
			{
				for(int j=0; j<m_nHeight; j++)
				{
					int i;
					for(i=0; i<m_nWidth && !*p_mask; i++, p_mask++, p_img++);
					if(i==m_nWidth) continue;
					
					BYTE m1=0, m2=*p_mask, m3=m2;
					float v1=0, v2=*p_img, v3=v2;
					for(;i<m_nWidth-1; i++, p_mask++, p_img++)
					{
						if(!*p_mask) {m1=false; continue;}
						
						m2 = *p_mask;		v2 = *p_img;
						m3 = *(p_mask+1);	v3 = *(p_img+1);

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

					if(*p_mask && m1) *p_img = (v1 + *p_img + *p_img) / 3;
					p_mask ++; p_img ++;
				}

				p_img = vec[h].begin() + k * m_nSliceSize;
				p_mask = m_mask.begin() + k * m_nSliceSize;

				for(int i=0; i<m_nWidth; i++, p_mask++, p_img++)
				{
					int j;

					vector<BYTE>::iterator pp_mask = p_mask;
					vector<float>::iterator pp_img = p_img;
					
					for(j=0; j<m_nHeight && !*pp_mask; j++, pp_mask+=m_nWidth, pp_img+=m_nWidth);
					if(j==m_nHeight) continue;

					BYTE	m1 = false, m2 = *pp_mask, m3=m2;
					float	v1 = 0,	 v2 = *pp_img,    v3=v2;
					for(; j<m_nHeight-1; j++, pp_mask+=m_nWidth, pp_img+=m_nWidth)
					{
						if(!*pp_mask) {m1=false; continue;}

						m2 = *pp_mask;				v2 = *pp_img;
						m3 = *(pp_mask+m_nWidth);	v3 = *(pp_img+m_nWidth);
						
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
					if(*pp_mask && m1) *pp_img = (v1 + *pp_img + *pp_img) / 3;
				}
			}
		}
		
		progress.Step(niter+1);
	}

	if(m_nVerbose==PROGRESS) cerr << endl;
	if(m_b3D)
	{
		ZProgressIndicator progress("Lowpass Filtering (Z)", 40, m_nbLowpass/12, m_nVerbose==PROGRESS);

		for(int niter = 0 ; niter < m_nbLowpass/12; niter++)
		{
			for(int h=0; h<m_nChannel + m_nChannel*m_nChannel; h++)
			{
				vector<BYTE>::iterator p_mask = m_mask.begin();
				vector<float>::iterator p_img = vec[h].begin();

 				for(int i=0; i<m_nSliceSize; i++, p_img++, p_mask++)
				{
					vector<BYTE>::iterator pp_mask = p_mask;
					vector<float>::iterator pp_img = p_img;
					
					int k;
					for(k=0; k<m_nDepth && !*pp_mask; k++, pp_mask+=m_nSliceSize, pp_img+=m_nSliceSize);
					if(k==m_nDepth) continue;

					BYTE m1=0, m2=*pp_mask, m3=m2;
					float v1=0, v2=*pp_img, v3=v2;
					for(; k<m_nDepth-1; k++, pp_mask+=m_nSliceSize, pp_img+=m_nSliceSize)
					{
						if(!*pp_mask) continue;

						m2 = *pp_mask; 					v2 = *pp_img;
						m3 = *(pp_mask+m_nSliceSize);	v3 = *(pp_img+m_nSliceSize);
						
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
				}
			}
			progress.Step(niter+1);
		}
		if(m_nVerbose==PROGRESS) cerr << endl;
	}

	for(i=0; i<m_nChannel; i++)
	{
		for(int k=0; k<m_nSize; k++) meanres[k][i] = vec[i][k];
	}

	for(; i<m_nChannel + m_nChannel*m_nChannel; i++)
	{
		int m = (i - m_nChannel) / m_nChannel;
		int n = (i - m_nChannel) % m_nChannel;

		for(int k=0; k<m_nSize; k++) meaninvcov[k](m, n) = vec[i][k];
	}
}

int ZMultiMRISegmentation::ICM(float BB)
{
	int change = 0;
	ZImage<BYTE>::iterator	p_seg(m_Segment);
	vector<BYTE>::iterator p_mask = m_mask.begin();
	ZMatrix<float>::row_iterator p_post = m_post.pbegin(), p_prob = m_prob.pbegin();

	for (; p_seg.more(); p_seg++, p_mask++, p_post++, p_prob++)
	{
		if(!*p_mask) continue;

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
			(*p_post)[c] = (*p_prob)[c] + v;
			(*p_prob)[c] = v;
		
			if(min > (*p_post)[c])
			{
				min = (*p_post)[c];
				clas = c;
			}
		}

		if(*p_seg != clas) change++, *p_seg = BYTE(clas);
	}
	return change;
}

bool ZMultiMRISegmentation::DetermineValidSlices()
{
	m_nBegin = 10000;
	bool v1, v2=false;

	int pos = 0;
	for(int k=0; k<m_nDepth; k++)
	{
		int nn = 0;
		for(int i=0; i<m_nSliceSize; i++, pos++)
		{
			bool flag = true;
			for(int c=0; c<m_nChannel; c++) 
			  if(m_Mri(pos, c) < 10) flag = false;  // MJ changed from 10 to 1 - test : ideally should be a robust value (e.g. 0.1 * (98 - 2 percentile) )

			if(flag) nn ++;
		}
		if(float(nn) / m_nSliceSize > 0.25) v1 = true;
		else v1 = false;

		if(v1 && !v2) m_nBegin = k;
		if(v2 && !v1) {	m_nEnd = k; break; }
		v2 = v1;
	}

	if(m_nBegin == 10000)
	{
#ifdef _WINDOWS
		AfxMessageBox("Empty input image. Nothing to do!");
#else
		cerr << "Empty input image. Nothing to do!\n";
#endif
		return false;
	}
	if(m_nEnd == 0) m_nEnd = m_nDepth-1;
	else m_nEnd --;

	if(m_nVerbose>=SLICE) cerr << "Parameter estimation will be performed from slice " 
						<< m_nBegin << " to " << m_nEnd << endl;

	return true;
}

ZGrayFloatImage* ZMultiMRISegmentation::GetBiasField(int channel)
{
	if(channel >= m_nChannel || channel < 0) 
	{
		ZWarning("ZMultiMRISegmentation::GetBiasField", "channel out of range!");
		return 0;
	}

	if(!m_FloatImage) m_FloatImage = new ZGrayFloatImage((*m_images)[0]->Width(), (*m_images)[0]->Height(), (*m_images)[0]->Depth());
	else m_FloatImage->Clear();
	m_FloatImage->SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	m_FloatImage->FillPixels(1);
	m_FloatImage->SetROI(m_left, m_right, m_top, m_bottom, m_front, m_back);

	ZImage<float>::iterator p_bias(*m_FloatImage);
	ZMatrix<float>::row_iterator p_mbias = m_Bias.pbegin();

	for(; p_bias.more(); p_bias++, p_mbias++) *p_bias = float(exp(-(*p_mbias)[channel]));

	m_FloatImage->FullROI();
	return m_FloatImage;
}

ZGrayFloatImage* ZMultiMRISegmentation::GetProbability(int wclass)
{
	if(wclass >= m_nbClasses || wclass < 0) 
	{
		ZWarning("ZMultiMRISegmentation::GetProbility", "class out of range!");
		return 0;
	}

	if(!m_FloatImage) m_FloatImage = new ZGrayFloatImage((*m_images)[0]->Width(), (*m_images)[0]->Height(), (*m_images)[0]->Depth());
	else m_FloatImage->Clear();

	m_FloatImage->SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	m_FloatImage->SetROI(m_left, m_right, m_top, m_bottom, m_front, m_back);

	ZImage<float>::iterator p_prob(*m_FloatImage);
	ZMatrix<float>::row_iterator p_post = m_post.pbegin();

	for(; p_prob.more(); p_prob++, p_post++) *p_prob = (*p_post)[wclass];

	m_FloatImage->FullROI();
	return m_FloatImage;
}

ZGrayFloatImage* ZMultiMRISegmentation::GetRestored(int channel)
{
	if(channel >= m_nChannel || channel < 0) 
	{
		ZWarning("ZMultiMRISegmentation::GetRestored", "channel out of range!");
		return 0;
	}

	if(!m_FloatImage) m_FloatImage = new ZGrayFloatImage((*m_images)[0]->Width(), (*m_images)[0]->Height(), (*m_images)[0]->Depth());
	else m_FloatImage->Clear();
	m_FloatImage->SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	m_FloatImage->SetROI(m_left, m_right, m_top, m_bottom, m_front, m_back);

	ZImage<float>::iterator p_res(*m_FloatImage);
	ZMatrix<float>::row_iterator p_mri= m_Mri.pbegin();
	ZMatrix<float>::row_iterator p_bias = m_Bias.pbegin();
	vector<BYTE>::iterator p_mask = m_mask.begin();
	for(; p_res.more(); p_res++, p_mri++, p_bias++, p_mask++) 
		if(*p_mask)
			*p_res = float(exp((*p_mri)[channel]));

	m_FloatImage->FullROI();
	return m_FloatImage;
}

void ZMultiMRISegmentation::PVE()
{
	static long seed=-1;

	if(m_nVerbose > QUIET) cerr << "Estimating PVE...\n";
	
// free some memory
	m_prob.clear();	m_apriori.clear();

	InitializePVE();

	ZVector<float> npve(m_nbClasses);

	float total;
	ZProgressIndicator progress("HMRF-EM for PVE", 40, 50, m_nVerbose==PROGRESS);
	for(int mrfiter=0; mrfiter<50; mrfiter++)
	{
		int c;
		total = 0;
		ZMatrix<float>::iterator p_mri = m_Mri.begin();
		ZMatrix<float>::iterator p_pve = m_PVE.begin();
		ZVector<float>::iterator p_eng = m_energy.pbegin();
		vector<BYTE>::iterator p_mask = m_mask.begin();
		for(; p_mri!=m_Mri.end(); p_mri++, p_eng++, p_pve++, p_mask++)
		{
			if(!*p_mask) continue;

			float aa = 1;
			for(c=1; c<m_nbClasses; c++) 
			{
				npve[c] = rand1(seed) * aa;
				aa -= npve[c];
			}
			npve[0] = aa;

			float cond = ComputeConditional(*p_mri, npve);
			float post = ComputePosterior(&(*p_mri), &(*p_pve), p_mask, npve);
			float eng = cond + post;

			float change = *p_eng - eng;
			if(change>0) *p_eng = eng, *p_pve = npve;

			total += *p_eng;
		}

		if(m_bPVEPara)
		{
			m_mean = 0; tv = 0;
			fill(m_covariance.begin(), m_covariance.end(), 0);
			
			ZMatrix<float>::iterator p_post = m_PVE.begin();
			p_mri = m_Mri.begin();
			p_mask = m_mask.begin();

			for(; p_mri<m_Mri.end(); p_mri+=m_nSampleRatio, p_post+=m_nSampleRatio, p_mask+=m_nSampleRatio)
			{
				if(!*p_mask) continue;

				for(c=0; c<m_nbClasses; c++)
				{
					float ppp = (*p_post)[c];
					for(int i=0; i<m_nChannel; i++) m_mean[c][i] += (*p_mri)[i] * ppp;
					tv[c] += ppp;
				}
			}

			for(c=0; c<m_nbClasses; c++) m_mean[c] /= tv[c];

			p_post = m_PVE.begin();
			p_mri = m_Mri.begin();
			p_mask = m_mask.begin();
			for(; p_mri<m_Mri.end(); p_mri+=m_nSampleRatio, p_post+=m_nSampleRatio, p_mask+=m_nSampleRatio)
			{
				if(!*p_mask) continue;

				for(c=0; c<m_nbClasses; c++)
				{
					tvc = *p_mri; tvc -= m_mean[c];

					float ppp = (*p_post)[c];
					for(int i=0; i<m_nChannel; i++)
					for(int j=0; j<m_nChannel; j++)
						m_covariance[c](i,j) += tvc[i] * tvc[j] * ppp;
				}
			}
				
			for(c=0; c<m_nbClasses; c++)
			{
				m_covariance[c] /= tv[c];
				m_icovariance[c] = inv(m_covariance[c]);
				m_stddev[c] = sqrt(Abs(det(m_covariance[c])));
			}
		}
		progress.Step(mrfiter+1);
	}
	if(m_nVerbose >= PARA) cerr << "Final energy " << total << endl; 

	if(m_nVerbose >= PARA) 
	{
		cerr << "    estimated parameters:\n";
		for(int c=0; c<m_nbClasses; c++)
			cerr << "\t" << exp(m_mean[c]) << endl;
		cerr << endl;
	}
}

void ZMultiMRISegmentation::InitializePVE()
{
	m_PVE.resize(m_nSize, m_nbClasses);
	m_energy.resize(m_nSize);

	ZMatrix<float>::iterator p_mri = m_Mri.begin();
	ZMatrix<float>::iterator p_pve = m_PVE.begin();
	vector<BYTE>::iterator p_mask = m_mask.begin();
	/*	for(; p_mri!=m_Mri.end(); p_mri++, p_pve++, p_mask++)
	{
		if(!*p_mask) continue;

		double	norm = 0; 
		for (int c = 0; c < m_nbClasses; c++)
		{
			double v = 0; 
			tvc = *p_mri; tvc -= m_mean[c];

			for(int i=0; i<m_nChannel; i++)
			{
				float u = 0;
				for(int j=0; j<m_nChannel; j++) u += tvc[j] * m_icovariance[c](j, i);
				v += u * tvc[i];
			}

			(*p_pve)[c] = v = float(exp(-v/2) / m_stddev[c] * M_1_SQRT_2PI);
			norm += v;
		}

		if(norm < TINY) *p_mask = 0; else *p_pve /= norm;
		}*/

        ZMatrix<float>::row_iterator p_post=m_post.pbegin();
        for(; p_mri!=m_Mri.end(); p_mri++, p_pve++, p_mask++, p_post++)
	  {
	    if(!*p_mask) continue;

	    for (int c = 0; c < m_nbClasses; c++)
		(*p_pve)[c]=(*p_post)[c];
	  }


	float total=0; p_pve=m_PVE.begin();
	p_mask = m_mask.begin(); p_mri = m_Mri.begin();
	ZVector<float>::iterator p_eng = m_energy.pbegin();
	for(; p_mri!=m_Mri.end(); p_mri++, p_mask++, p_pve++, p_eng++)
	{
		if(!*p_mask) continue;

		float cond = ComputeConditional(*p_mri, *p_pve);
		float post = ComputePosterior(&(*p_mri), &(*p_pve), p_mask, *p_pve);
		total += *p_eng = cond + post;
	}
	if(m_nVerbose >= PARA) cerr << "initial energy " << total << endl; 
}

float ZMultiMRISegmentation::ComputeConditional(const ZVector<float>& mri, const ZVector<float>& pve)
{
	int c;
	for(c=0; c<m_nbClasses; c++) tv[c] = pve[c] * pve[c];

	// calculate pixel mean and covariance

	m1 = m_covariance[0]; m1 *= tv[0];
	for(c=1; c<m_nbClasses; c++) m1 += m_covariance[c] * tv[c];
	
	m2 = inv(m1);

	for(c=0; c<m_nChannel; c++)
	{
		float a = 0;
		for(int i=0; i<m_nbClasses; i++) a+= pve[i] * m_mean(i,c);
		tvc[c] = mri[c] - a;
	}

	float cond = 0;
	for(c=0; c<m_nChannel; c++)
	{
		float v = 0;
		for(int i=0; i<m_nChannel; i++) v += tvc[i] * m2(i, c);
		cond += v * tvc[c];
	}

	cond = (cond + log(Abs(det(m1)))) / 2;

	return cond;
}

float ZMultiMRISegmentation::ComputePosterior(const ZVector<float>* p_mri, const ZVector<float>* p_pve, 
											vector<BYTE>::const_iterator p_mask, const ZVector<float>& rpv)
{
#define _POST(offset) \
	if(p_mask[offset]) \
	{ const ZVector<float>& pve = *(p_pve+offset); \
	  for(int c=0; c<m_nbClasses; c++) post += (pve[c] - rpv[c]) * (pve[c] - rpv[c]); }

	const ZVector<float>* pmri;
	float post=0;

	if(m_b3D)
	{
		if((pmri = p_mri - m_nSliceSize) > m_Mri.pbegin() && *(p_mask-m_nSliceSize)) 
		{
			const ZVector<float>& pve = *(p_pve - m_nSliceSize);
			float dif = 0;
			for(int c=0; c<m_nbClasses; c++) dif += (pve[c] - rpv[c]) * (pve[c] - rpv[c]);
			post += dif * diaB;
		}
		if((pmri = p_mri + m_nSliceSize) < m_Mri.pend() && *(p_mask+m_nSliceSize)) 
		{
			const ZVector<float>& pve = *(p_pve + m_nSliceSize);
			float dif = 0;
			for(int c=0; c<m_nbClasses; c++) dif += (pve[c] - rpv[c]) * (pve[c] - rpv[c]);
			post += dif * diaB;
		}
	}

	if((p_mri - 1) > m_Mri.pbegin()) _POST(-1); 
	if((p_mri + 1) < m_Mri.pend()) _POST(1); 

	if((p_mri - m_nWidth) > m_Mri.pbegin()) _POST(-m_nWidth); 
	if((p_mri + m_nWidth) < m_Mri.pend()) _POST(m_nWidth); 

	if((p_mri-1-m_nWidth) > m_Mri.pbegin()) _POST(-1-m_nWidth); 
	if((p_mri+1-m_nWidth) > m_Mri.pend()) _POST(1-m_nWidth); 

	if((p_mri-1+m_nWidth) < m_Mri.pbegin()) _POST(-1+m_nWidth); 
	if((p_mri+1+m_nWidth) < m_Mri.pend()) _POST(1+m_nWidth); 

	return post * PVEB;
}

ZGrayFloatImage* ZMultiMRISegmentation::GetPVE(int wclass)
{
	if(wclass >= m_nbClasses || wclass < 0) 
	{
		ZWarning("ZMultiMRISegmentation::GetPVE", "class out of range!");
		return 0;
	}

	if(!m_FloatImage) m_FloatImage = new ZGrayFloatImage((*m_images)[0]->Width(), (*m_images)[0]->Height(), (*m_images)[0]->Depth());
	else m_FloatImage->Clear();

	m_FloatImage->SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	m_FloatImage->SetROI(m_left, m_right, m_top, m_bottom, m_front, m_back);

	ZImage<float>::iterator p_ipve(*m_FloatImage);
	ZMatrix<float>::row_iterator p_pve = m_PVE.pbegin();

	for(; p_ipve.more(); p_ipve++, p_pve++) *p_ipve = (*p_pve)[wclass];

	m_FloatImage->FullROI();
	return m_FloatImage;
}
