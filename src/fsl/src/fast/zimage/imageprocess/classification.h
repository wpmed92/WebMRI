/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for classification.

#ifndef __CLASSIFICATION__H_
#define __CLASSIFICATION__H_

void Classification(const ZImageBase& image, ZGrayByteImage& seg, int nclass, bool bNoZero, float* mean, float* stddev=0, float* prior=0);
void GetStatistics(const ZImageBase& image, const ZGrayByteImage& seg, int nclass, bool bNoZero, float* mean, float* stddev, float* prior);
void GetStatistics(const ZImageBase& image, int nclass, bool bNoZero, float* thres, float* mean, float* stddev, float* prior);
int EM(const ZImageBase& image, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior=0, std::ostream* report=0);
int EM(const ZImageBase& image, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior=0, std::ostream* report=0);
int TreeEM(const ZImageBase& image, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, std::ostream* report=0);
int TreeEM(const ZImageBase& image, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, std::ostream* report=0);

int KMean(const ZImageBase& image, int nclass, int NITER, bool bNoZero, float* mean, std::ostream* report=0);
int KMean(const ZImageBase& image, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, std::ostream* report=0);
int TreeKMean(const ZImageBase& image, int nclass, bool bNoZero, float* m, std::ostream* report=0);
int TreeKMean(const ZImageBase& image, ZGrayByteImage& seg, int nclass, bool bNoZero, std::ostream* report=0);
int TreeKMean(const ZImageBase& image, ZGrayByteImage& seg, int nclass, bool bNoZero, float* mean, float* stddev, float* prior, std::ostream* report=0);

int MRF(const ZImageBase& image, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* priorm, bool ICM, float B=-1);
int HMRF_EM(const ZImageBase& image, ZGrayByteImage& seg, int nclass, bool bNoZero, float B=-1, float* m=0, float* s=0, float* p=0, std::ostream* report=0);

void Thresholding(const ZImageBase& image, ZGrayByteImage& seg, int nclass, const float* thres, bool bNoZero=false);
inline void Thresholding(const ZImageBase& image, ZGrayByteImage& seg, float thres, bool bNoZero=false)
{
	Thresholding(image, seg, 2, &thres, bNoZero);
}
inline void Thresholding(const ZImageBase& image, ZGrayByteImage& seg, float thres1, float thres2, bool bNoZero=false)
{
	float thres[2]; thres[0] = thres1; thres[1] = thres2;
	Thresholding(image, seg, 3, thres, bNoZero);
}
inline void Thresholding(const ZImageBase& image, ZGrayByteImage& seg, float thres1, float thres2, float thres3, bool bNoZero=false)
{
	float thres[3]; thres[0] = thres1; thres[1] = thres2; thres[2] = thres3;
	Thresholding(image, seg, 4, thres, bNoZero);
}

#endif
