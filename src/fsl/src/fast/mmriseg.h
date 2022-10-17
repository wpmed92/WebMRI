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

#ifndef __MULTIMRISEGMENTATION__H_
#define __MULTIMRISEGMENTATION__H_

#include <vector>
#include "imagecore.h"

// multi-spectral version of ZMRISegmentation
class ZMultiMRISegmentation
{
	enum VERBOSE { QUIET=0, MAIN=1, SLICE=2, ITER=3, PARA=4, PROGRESS=5 };

	int		m_nChannel;

	int		m_nWidth;			// image width
	int		m_nHeight;			// image height
	int		m_nDepth;			// image depth (number of slices)
	int		m_nSliceSize;		// width x height
	int		m_nSize;			// m_nSliceSize x m_nDepth
	
	int		m_nBegin;			// the first slice to be used for parameter estimation
	int		m_nEnd;				// the last slice to be used for parameter estimation

	float	pixdim[3];			// pixel dimmension in mm.

	bool	m_b3DLowpass;		// whether to do 3D lowpass filtering

	VERBOSE	m_nVerbose;			// verbose mode 

	BYTE	neighbour[27];		// for MRF neighbours of a pixel (second-order)
	float	B;					// MRF beta parameter
	float	diaB;				// MRF beta parameter used for diagonal neighbours
	float	ddiaB;				// MRF beta parameter used for neighbours in different slices
	float	TB;					// maximum prior of a pixel
	int		offset3d;			// offset of the first pixel in a 3x3x3 cube with respect to the center pixel
	int		offset2d;			// offset of the first pixel in a 3x3 cube with respect to the center pixel

	int		m_left, m_right, m_top, m_bottom, m_front, m_back;	// ROI of the images

	const std::vector<ZImageBase*>* m_images;

	std::vector<BYTE>	m_mask;

	ZGrayFloatImage		*m_FloatImage;

	void			EMloop(void);
	int ICM(float BB);
	void Lowpass(ZMatrix<float>& meanres, std::vector< ZMatrix<float> >& meaninvcov);
	void Initialize();
	// determine the range of slices for parameter estimation
	bool	DetermineValidSlices();

//////////////////////////////////////////////////////
// PVE members
//////////////////////////////////////////////////////
	void	PVE();			// PVE classification
	void	InitializePVE();
	float	ComputeConditional(const ZVector<float>& mri, const ZVector<float>& pve);
	float	ComputePosterior(const ZVector<float>* p_mri, const ZVector<float>* p_pve, 
				std::vector<BYTE>::const_iterator p_mask, const ZVector<float>& rpv);
public:

	ZMatrix<float>			m_Mri, m_Bias;		//m_nSize x m_nChannel
	ZGrayByteImage			m_Segment;
	
	ZMatrix<float> m_prob;			// likelihood probability map
	ZMatrix<float> m_post;			// posterior probability map
	ZVector<float> m_prior;
	ZMatrix<BYTE> m_apriori;		// apriori probability map

	int m_bapriori;

	int		m_nbClasses;			// number of classes
	int		m_nbIter;				// number of iterations
	int		m_nbLowpass;			// number of lowpass filtering for bias field

	UINT	m_nSampleRatio;			// the sample ratio of when doing parameter updating

	bool	m_b3D;					// 3D segmentation or 2D segmentation
	bool	m_b2DImage;				// the original image is a 2D image
	bool	m_bBias;				// whether to correct for the bias field

	ZMatrix<float>  m_mean;			// m_nbClasses x m_nChannel
	// m_nChannel x m_nChannel
	std::vector< ZMatrix<float> >	m_covariance, m_icovariance;		
	ZVector<float>			m_stddev;	// 1/||m_covariance||^2
	ZVector<UINT>			m_number;

	ZMatrix<float>	m_inimean;
	std::vector< ZMatrix<float> >	m_inicovariance;

	ZMultiMRISegmentation();
	~ZMultiMRISegmentation() { delete m_FloatImage; }
	
	void Create(const std::vector<ZImageBase*>& images, int nclass, bool b2d,
				int nbiter,	int nblowpass, float neighbor, UINT sampleratio, bool bBias=true);

	bool Segment(bool pve, const ZGrayByteImage* gm=0, const ZGrayByteImage* wm=0, const ZGrayByteImage* csf=0, int bapriori=0);

	ZGrayFloatImage* GetRestored(int channel);
	ZGrayFloatImage* GetBiasField(int channel);
	ZGrayFloatImage* GetProbability(int wclass);

	void SetVerboseMode(int verbose=1)  
	{ 
		if(verbose<0 || verbose>5) verbose = 1;
		m_nVerbose = VERBOSE(verbose); 
	}

//////////////////////////////////////////////////////
// PVE members
//////////////////////////////////////////////////////
	ZVector<float>	tv;				// temporary vector (m_nbClasses)
	ZMatrix<float>	m_PVE;			// m_nSize x m_nbClasses
	ZVector<float>	m_energy;		// energy
	ZGrayFloatImage m_PVEimage;		// PVE image
	bool			m_bPVEPara;		// whether updating parameters during PVE
	float			PVEB;			// beta for PVE MRF

	void SetPVEPara(bool bPara, float beta) { m_bPVEPara = bPara; PVEB = beta; }
	ZGrayFloatImage* GetPVE(int wclass);
};

#endif
