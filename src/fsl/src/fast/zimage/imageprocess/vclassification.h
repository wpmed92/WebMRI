/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for multispectral classification.

#ifndef __VCLASSIFICATION__H_
#define __VCLASSIFICATION__H_

void Reform(const std::vector<ZImageBase*>& images, ZMatrix<float>& mat);

void Classification(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, 
					ZMatrix<float>& mean, std::vector< ZMatrix<float> >* covariance=0, float* prior=0);
void GetStatistics(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, 
				   ZMatrix<float>& mean, std::vector< ZMatrix<float> >& covariance, float* prior);
int EM(const ZMatrix<float>& images, int nclass, int NITER, bool bNoZero, ZMatrix<float>& mean, 
	   std::vector< ZMatrix<float> >& covariance, float* prior);
int EM(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, 
	   ZMatrix<float>& mean, std::vector< ZMatrix<float> >& covariance, float* prior);
int TreeEM(const ZMatrix<float>& images, int nclass, int NITER, bool bNoZero, ZMatrix<float>& mean, 
		   std::vector< ZMatrix<float> >& covariance, float* prior, std::ostream* report=0);
int TreeEM(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, 
		   ZMatrix<float>& mean, std::vector< ZMatrix<float> >& covariance, float* prior, std::ostream* report=0);

int KMean(const ZMatrix<float>& images, int nclass, int NITER, bool bNoZero, ZMatrix<float>& mean, std::ostream* report=0);
int KMean(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, 
		  ZMatrix<float>& mean, std::vector< ZMatrix<float> >& covariance, float* prior, std::ostream* report=0);
int TreeKMean(const ZMatrix<float>& images, int nclass, bool bNoZero, ZMatrix<float>& mean, std::ostream* report=0);
int TreeKMean(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, std::ostream* report=0);
int TreeKMean(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, ZMatrix<float>& mean, 
			  std::vector< ZMatrix<float> >& covariance, float* prior, std::ostream* report=0);

inline int TreeKMean(const std::vector<ZImageBase*>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, std::ostream* report=0)
{
	ZMatrix<float> mat;  Reform(images, mat);
	seg.Create(images[0]->Width(), images[0]->Height(), images[0]->Depth());
	return TreeKMean(mat, seg, nclass, bNoZero, report);
}

inline int TreeKMean(const std::vector<ZImageBase*>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, ZMatrix<float>& mean, 
			  std::vector< ZMatrix<float> >& covariance, float* prior, std::ostream* report=0)
{
	ZMatrix<float> mat;  Reform(images, mat);
	seg.Create(images[0]->Width(), images[0]->Height(), images[0]->Depth());
	return TreeKMean(mat, seg, nclass, bNoZero, mean, covariance, prior, report);
}

#endif
