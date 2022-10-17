/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */


#include <cmath>
#include <vector>
#include <iostream>

#include "common.h"
#include "imagecore.h"
#include "zmath.h"
#include "classification.h"

using namespace std;

void Classification(const ZImageBase& image, ZGrayByteImage& seg, int nclass, bool bNoZero, float* mean, float* stddev, float* prior)
{
	if(nclass < 2)
	{
		ZError("Classification", "The number of classes must be more than 1!");
		return;
	}

	seg.Create(image.Width(), image.Height(), image.Depth());

	float	pixdim[3];
	image.GetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	seg.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);

	PBYTE p_seg = seg.GetBuffer();
	if(stddev==0)
	{
		for(UINT k=0; k<image.Depth(); k++)
		for(UINT j=0; j<image.Height(); j++)
		for(UINT i=0; i<image.Width(); i++, p_seg++)
		{
			float v = image.GetIntensity(i, j, k);
			if(!bNoZero || v)
			{
				int cc = 0; float max = 1e10;
				for(int c=0; c<nclass; c++)
				{
					float e = v - mean[c];
					if(max > Abs(e)) max = Abs(e), cc=c;
				}
				*p_seg = bNoZero ? cc+1 : cc;
			}
		}
		return;
	}

	vector<float> variance(nclass);
	for(int c=0; c<nclass; c++) variance[c] = stddev[c] * stddev[c];

	for(UINT k=0; k<image.Depth(); k++)
	for(UINT j=0; j<image.Height(); j++)
	for(UINT i=0; i<image.Width(); i++, p_seg++)
	{
		float v = image.GetIntensity(i, j, k);
		if(!bNoZero || v)
		{
			double prob, max=0; int cc = 0;
			for(int c=0; c<nclass; c++)
			{
				double e = v - mean[c];
				e = 0.5f * e * e / variance[c];
				prob = exp(-e) / stddev[c] * M_1_SQRT_2PI;
				if(prior) prob *= prior[c];
				if(max < prob) max = prob, cc = c;
			}
			*p_seg = bNoZero ? cc+1 : cc;
		}
	}
}

void Thresholding(const ZImageBase& image, ZGrayByteImage& seg, int nclass, const float* thres, bool bNoZero)
{
	if(nclass < 2)
	{
		ZError("Classification", "The number of classes must be more than 1!");
		return;
	}

	seg.Create(image.Width(), image.Height(), image.Depth());

	float	pixdim[3];
	image.GetPixelDim(pixdim[0], pixdim[1], pixdim[2]);
	seg.SetPixelDim(pixdim[0], pixdim[1], pixdim[2]);

	PBYTE p_seg = seg.GetBuffer();
	for(UINT k=0; k<image.Depth(); k++)
	for(UINT j=0; j<image.Height(); j++)
	for(UINT i=0; i<image.Width(); i++, p_seg++)
	{
		float v = image.GetIntensity(i, j, k);
		if(bNoZero && !v) continue;

		int c;
		for(c=0; c<nclass-1; c++) if(v < thres[c]) break;
		*p_seg = bNoZero?c+1:c;
	}
}

void GetStatistics(const ZImageBase& image, const ZGrayByteImage& seg, int nclass, bool bNoZero, float* mean, float* stddev, float* prior)
{
	if(seg.Size() != image.Size()) return;

	UINT depth=image.Depth(), height=image.Height(), width=image.Width();
	int c;

	for(c=0; c<nclass; c++) mean[c] = stddev[c] = prior[c] = 0;

	ZGrayByteImage::const_iterator p_seg(seg);
	UINT k;
	for(k=0; k<depth; k++)
	for(UINT j=0; j<height; j++)
	for(UINT i=0; i<width; i++, p_seg++)
	{
		int c = *p_seg;
		if(bNoZero) if(!c) continue; else c--;

		mean[c] += image.GetIntensity(i, j, k);
		prior[c] ++;
	}

	float norm = 0;
	for(c=0; c<nclass; c++) { norm += prior[c]; if(prior[c] > 0) mean[c] /= prior[c]; }

	p_seg.reset();
	for(k=0; k<depth; k++)
	for(UINT j=0; j<height; j++)
	for(UINT i=0; i<width; i++, p_seg++)
	{
		int c = *p_seg;
		if(bNoZero) if(!c) continue; else c--;

		float inten=image.GetIntensity(i, j, k);
		stddev[c] += (inten-mean[c]) * (inten-mean[c]);
	}

	for(c=0; c<nclass; c++) 
	{
		if(prior[c] > 0) stddev[c] = float(sqrt(stddev[c] / prior[c])); 
		prior[c] /= norm;
	}
}

void GetStatistics(const ZImageBase& image, int nclass, bool bNoZero, float* thres, float* mean, float* stddev, float* prior)
{
	ZGrayByteImage seg;
	Thresholding(image, seg, nclass, thres);
	GetStatistics(image, seg, nclass, bNoZero, mean, stddev, prior);
}

static int EMHist(const ZGrayByteImage& image, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, ostream* report)
{
	float *weight = new float[nclass];
	UINT  hist[256]; int c, size=0, iter;

	for(c=0; c<nclass; c++)	weight[c] = prior ? prior[c] : 1.0f / nclass;

	memset(hist, 0, 256*sizeof(UINT));
	ZGrayByteImage::const_iterator p_image(image);
	for(; p_image.more(); p_image++)
	{
		if(!bNoZero || *p_image) { size++; hist[*p_image]++; }
	}

	vector<double*> prob(nclass);
	vector<float> variance(nclass);
	for(c=0; c<nclass; c++)
	{
		prob[c] = new double[256];
		variance[c] = stddev[c] * stddev[c];
	}

	vector<int> label(256);
	for(iter=0; iter<NITER; iter++)
	{
		int i, c, change=0;

		for(i=0; i<256; i++)
		{
			if(hist[i]==0) continue;

			double norm = 0, max=0; int cc=0;
			for(c=0; c<nclass; c++)
			{
				double val = i - mean[c];
				val = 0.5f * val * val / variance[c];
				val = exp(-val) / stddev[c] * M_1_SQRT_2PI;
				prob[c][i] = val * weight[c];
				if(prob[c][i] > max) max = prob[c][i], cc=c;
				norm += prob[c][i]; 
			}
			if(label[i] != cc) change++, label[i] = cc;
			for (c = 0; c < nclass; c++) prob[c][i] /= norm;
		}

		if(report)
		{
			*report  << setprecision (4) << setiosflags (ios::fixed)
				<< "iter=" << iter << "\t" ;
			for(c=0; c<nclass; c++)
				*report	<< "U" << c << "=" << mean[c] << "\t" 
					<< "S" << c << "=" << stddev[c] << "\t"
					<< "P" << c << "=" << weight[c] << "\t";
					
			*report << "change=" << change << endl;
		}
		if(change == 0) break;

		vector<float> nor(nclass);

		for(c=0; c<nclass; c++) mean[c] = stddev[c] = variance[c] = nor[c] = 0;

		for(i=0; i<256; i++)
		{
			for (c = 0; c < nclass; c++)
			{
				double tmp = prob[c][i] * hist[i];
				mean[c] += tmp * i;
				nor[c] += tmp;
			}
		}

		for(c=0; c<nclass; c++) if(nor[c] > 1) mean[c] /= nor[c];

		for(i=0; i<256; i++)
		{
			if(hist[i] > 0)
			for (c = 0; c < nclass; c++)
			{
				float v = i - mean[c];
				variance[c] += (prob[c][i] * v * v) * hist[i];
			}
		}

		for(c=0; c<nclass; c++)
		{
			variance[c] /= nor[c];
			stddev[c] = float(sqrt(variance[c]));
			weight[c] = nor[c] / size;
		}
	}

	for(c=0; c<nclass; c++) { if(prior) prior[c] = weight[c]; delete []prob[c]; }
	delete []weight;

	return iter;
}

int EM(const ZImageBase& image, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, ostream* report)
{
	if(nclass < 2)
	{
		ZError("EM", "The number of classes must be more than 1!");
		return 0;
	}

	if(is8bit(image)) return EMHist((const ZGrayByteImage&)image, nclass, NITER, bNoZero, mean, stddev, prior, report);

	float *weight = new float[nclass];

	int c, iter, size = 0;
	for(c=0; c<nclass; c++)	weight[c] = prior ? prior[c] : 1.0f / nclass;

	ZGrayFloatImage* pimage;
	if(isFloat(image)) pimage = (ZGrayFloatImage*)(&image); 
	else { pimage = new ZGrayFloatImage; image.SendPixelsTo(*pimage); }
	
	ZGrayFloatImage& fimage = *pimage;

	ZGrayFloatImage::const_iterator p_image(fimage);
	for(; p_image.more(); p_image++)
	{
		if(!bNoZero || *p_image) size ++;
	}

	vector<float*> prob(nclass);
	vector<float> variance(nclass);
	for(c=0; c<nclass; c++)
	{
		prob[c] = new float[size];
		variance[c] = stddev[c] * stddev[c];
	}

	double error = 1e20;
//	ZProgressIndicator progress("progress", 70, NITER);
	for(iter=0; iter<NITER; iter++)
	{
		int i, c;

		double tmp=0;
		for(i=0, p_image.reset(); p_image.more(); p_image++)
		{
			if(bNoZero && !*p_image) continue;
			double norm = 0;
			for(c=0; c<nclass; c++)
			{
				double val = *p_image - mean[c];
				val = 0.5 * val * val / variance[c];
				val = exp(-val) / stddev[c] * M_1_SQRT_2PI;
				prob[c][i] = val * weight[c];
				norm += prob[c][i]; 
			}
			tmp -= log(norm);
			for (c = 0; c < nclass; c++) prob[c][i] /= norm;
			i++;
		}

		if(report)
		{
			*report  << setprecision (4) << setiosflags (ios::fixed)
				<< "iter=" << iter << "\t" ;
			for(c=0; c<nclass; c++)
				*report	<< "U" << c << "=" << mean[c] << "\t" 
					<< "S" << c << "=" << stddev[c] << "\t"
					<< "P" << c << "=" << weight[c] << "\t";
					
			*report << "error=" << tmp << endl;
		}
		tmp /= size;
		if(Abs(error-tmp) < 1e-10) break; error = tmp;

		vector<float> nor(nclass);

		for(c=0; c<nclass; c++) mean[c] = stddev[c] = variance[c] = nor[c] = 0;

		for(i=0, p_image.reset(); p_image.more(); p_image++)
		{
			if(bNoZero && !*p_image) continue;
			for (c = 0; c < nclass; c++)
			{
				mean[c] += prob[c][i] * *p_image;
				nor[c] += prob[c][i];
			}
			i++;
		}

		for(c=0; c<nclass; c++) if(nor[c] > 1) mean[c] /= nor[c];

		for(i=0, p_image.reset(); p_image.more(); p_image++)
		{
			if(bNoZero && !*p_image) continue;
			for (c = 0; c < nclass; c++)
			{
				float v = *p_image - mean[c];
				variance[c] += prob[c][i] * v * v;
			}
			i++;
		}

		for(c=0; c<nclass; c++)
		{
			variance[c] /= nor[c];
			stddev[c] = float(sqrt(variance[c]));
			weight[c] = nor[c] / size;
		}
//		progress.Step(iter);
	}

	for(c=0; c<nclass; c++) { if(prior) prior[c] = weight[c]; }
	delete []weight;

	if(!isFloat(image)) delete pimage;
	return iter;
}

int EM(const ZImageBase& image, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, ostream* report)
{
	int iter = EM(image, nclass, NITER, bNoZero, mean, stddev, prior, report);
	Classification(image, seg, nclass, bNoZero, mean, stddev, prior);
	return iter;
}

int TreeEM(const ZImageBase& image, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, ostream* report)
{
	if(nclass < 2)
	{
		ZError("TreeEM", "The number of classes must be more than 1!");
		return 0;
	}

	if(nclass > 10)
	{
		ZWarning("TreeEM", "The number of classes should be less than 10!");
		nclass = 10;
	}

	int i, c;
	vector<float> m(nclass), s(nclass), p(nclass);

	float min = 0, max = 256; image.MinMax(min, max); 
	image.Statistics(m[0], s[0], bNoZero);
	
	if(report) *report << "starting with 1 class -- mean: " << m[0] << "  stddev: " << s[0] << endl;

	int totaliter=0;

	for(int nc=1; nc<nclass; nc++)
	{
		int splitnode = max_element(&(*s.begin()), &s[nc])-&(*s.begin());
		for(i=nc-1; i>splitnode; i--) 
		{
			m[i+1] = m[i];
			s[i+1] = s[i] * nc / (nc+1);
		}
		float mm = m[splitnode], ss = s[splitnode];
		m[splitnode] = mm - ss/2;
		s[splitnode] = ss * nc / (nc+1);
		if(splitnode > 0) 
		{ 
			if(m[splitnode] < m[splitnode-1]) 
				m[splitnode] = (mm + m[splitnode-1]) / 2;
		}
		else if(m[splitnode] < min) m[splitnode] = min + 5;

		m[splitnode+1] = mm + ss/2;
		s[splitnode+1] = ss / 2;
		if(splitnode < nc-1) 
		{ 
			if(m[splitnode+1] > m[splitnode+2]) 
				m[splitnode+1] = (mm + m[splitnode+2]) / 2;
		}
		if(m[splitnode+1] > max) m[splitnode+1] = max - 5;

		for(i=0; i<nc+1; i++) p[i] = 1.0f / (nc+1);
		if(nc+1 == nclass) NITER *= 2;
		totaliter += EM(image, nc+1, NITER, bNoZero, &(*m.begin()), &(*s.begin()), &(*p.begin()));

		if(report)
		{
			*report << "splite node " << splitnode << " total " << nc+1 << " classes" << endl;
			*report << "mean/stddev/prior for each classes:\t";
			for(i=0; i<nc; i++)
				*report << mean[i] << " / " << stddev[i] << " / " << prior[i] << endl << "\t\t\t\t\t";
			*report << mean[i] << " / " << stddev[i] << " / " << prior[i] << endl;
			*report << "Total number of iterations so far: " << totaliter << endl << endl;
		}
	}

	for(c=0; c<nclass; c++) mean[c] = m[c];
	for(c=0; c<nclass; c++) stddev[c] = s[c];
	for(c=0; c<nclass; c++) prior[c] = p[c];

	return totaliter;
}

int TreeEM(const ZImageBase& image, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, ostream* report)
{
	int iter = TreeEM(image, nclass, NITER, bNoZero, mean, stddev, prior, report);
	Classification(image, seg, nclass, bNoZero, mean, stddev, prior);
	return iter;
}

static int KMeanHist(const ZGrayByteImage& image, int nclass, int NITER, bool bNoZero, float* mean, ostream* report)
{
	UINT  hist[256];

	int iter, size=0;
	memset(hist, 0, 256*sizeof(UINT));
	ZGrayByteImage::const_iterator p_image(image);
	for(; p_image.more(); p_image++)
	{
		if(!bNoZero || *p_image) { size++; hist[*p_image]++; }
	}

	vector<int> num(nclass), label(256);
	vector<float> nmean(nclass);

	for(iter=0; iter<NITER; iter++)
	{
		int i, c, change = 0;

		for(c=0; c<nclass; c++) num[c] = 0, nmean[c] = 0;
		for(i=0; i<256; i++)
		{
			if(hist[i]==0) continue;

			float min = 1e10; int cc = 0;
			for(c=0; c<nclass; c++)
			{
				float val = i - mean[c];
				val = Abs(val);
				if(min > val) min = val, cc = c;
			}

			if(label[i] != cc) change ++, label[i] = cc;
			num[cc] += hist[i]; nmean[cc] += i * hist[i];
		}
		for(c=0; c<nclass; c++) mean[c] = nmean[c]/num[c];

		if(report)
		{
			*report  << setprecision (4) << setiosflags (ios::fixed)
				<< "iter=" << iter << "\t" ;
			for(c=0; c<nclass; c++)
				*report << "U" << c << "=" << mean[c] << "\t" 
					<< "P" << c << "=" << float(num[c])/size << "\t";
					
			*report << "change=" << change << endl;
		}
		if(change == 0) break;
	}

	return iter;
}

int KMean(const ZImageBase& image, int nclass, int NITER, bool bNoZero, float* mean, ostream* report)
{
	if(nclass < 2)
	{
		ZError("KMean", "The number of classes must be more than 1!");
		return 0;
	}

	if(is8bit(image)) return KMeanHist((const ZGrayByteImage&)image, nclass, NITER, bNoZero, mean, report);

	int iter, size=0;
	ZGrayFloatImage* pimage;
	if(isFloat(image)) pimage = (ZGrayFloatImage*)(&image); 
	else { pimage = new ZGrayFloatImage; image.SendPixelsTo(*pimage); }
	
	ZGrayFloatImage& fimage = *pimage;

	ZGrayFloatImage::const_iterator p_image(fimage);
	for(; p_image.more(); p_image++)
	{
		if(!bNoZero || *p_image) size++;
	}

	if(report)
	{
		*report  << setprecision (4) << setiosflags (ios::fixed)
			<< "initial settings:\t";
		for(int c=0; c<nclass; c++)
			*report << "U" << c << "=" << mean[c] << "\t";
		*report << endl; 
	}

	vector<UINT> num(nclass);
	vector<float> nmean(nclass);

	for(iter=0; iter<NITER; iter++)
	{
		int c;

		for(c=0; c<nclass; c++) num[c] = 0, nmean[c] = 0;
		for(p_image.reset(); p_image.more(); p_image++)
		{
			if(bNoZero && !*p_image) continue;

			float min = 1e10; int cc = 0;
			for(c=0; c<nclass; c++)
			{
				float val = *p_image - mean[c];
				val = Abs(val);
				if(min > val) min = val, cc = c;
			}

			num[cc] ++; nmean[cc] += *p_image;
		}
		bool quit=true;
		for(c=0; c<nclass; c++) 
		{
			nmean[c]/=num[c];
			if(mean[c] == nmean[c]) continue;
			quit = false;
			mean[c] = nmean[c];
		}

		if(report)
		{
			*report  << setprecision (4) << setiosflags (ios::fixed)
				<< "iter=" << iter << "\t" ;
			for(c=0; c<nclass; c++)
				*report << "U" << c << "=" << mean[c] << "\t" 
					<< "P" << c << "=" << float(num[c])/size << "\t";
			*report << endl;
		}
		if(quit) break; 
	}

	if(!isFloat(image)) delete pimage;
	return iter;
}

int KMean(const ZImageBase& image, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, ostream* report)
{
	int iter = KMean(image, nclass, NITER, bNoZero, mean, report);
	Classification(image, seg, nclass, bNoZero, mean);
	GetStatistics(image, seg, nclass, bNoZero, mean, stddev, prior);
	return iter;
}

int TreeKMean(const ZImageBase& image, int nclass, bool bNoZero, float* m, ostream* report)
{
	if(nclass < 2)
	{
		ZError("TreeKMean", "The number of classes must be more than 1!");
		return 0;
	}

	if(nclass > 10)
	{
		ZWarning("TreeKMean", "The number of classes should be less than 10!");
		nclass = 10;
	}

	int i;
	vector<float> mean(nclass), stddev(nclass), prior(nclass);

	float min = 0, max = 256; image.MinMax(min, max); 
	image.Statistics(mean[0], stddev[0], bNoZero);

	if(report) *report << "starting with 1 class -- mean: " << mean[0] << "  stddev: " << stddev[0] << endl;

	int iter = 0;
	ZGrayByteImage seg;
	for(int nc=1; nc<nclass; nc++)
	{
		int splitnode = max_element(&(*stddev.begin()), &stddev[nc])-&(*stddev.begin());
		for(i=nc-1; i>splitnode; i--) mean[i+1] = mean[i];
		float m = mean[splitnode], s = stddev[splitnode];
		mean[splitnode] = m - s;
		if(splitnode > 0) 
		{ 
			if(mean[splitnode] < mean[splitnode-1]) 
				mean[splitnode] = (m + mean[splitnode-1]) / 2;
		}
		else if(mean[splitnode] < min) mean[splitnode] = min + 5;

		mean[splitnode+1] = m + s;
		if(splitnode < nc-1) 
		{ 
			if(mean[splitnode+1] > mean[splitnode+2]) 
				mean[splitnode+1] = (m + mean[splitnode+2]) / 2;
		}
		if(mean[splitnode+1] > max) mean[splitnode+1] = max - 5;

		iter += KMean(image, seg, nc+1, 500, bNoZero, &(*mean.begin()), &(*stddev.begin()), &(*prior.begin()));

		if(report)
		{
			*report << "splite node " << splitnode << " total " << nc+1 << " classes" << endl;
			*report << "mean/stddev/prior for each classes:\t";
			for(i=0; i<nc; i++)
				*report << mean[i] << " / " << stddev[i] << " / " << prior[i] << endl << "\t\t\t\t\t";
			*report << mean[i] << " / " << stddev[i] << " / " << prior[i] << endl;
			*report << "Total number of iterations so far: " << iter << endl << endl;
		}
	}

	if(m) for(int c=0; c<nclass; c++) m[c] = mean[c];

	return iter;
}

int TreeKMean(const ZImageBase& image, ZGrayByteImage& seg, int nclass, bool bNoZero, ostream* report)
{
	vector<float> mean(nclass);
	int iter = TreeKMean(image, nclass, bNoZero, &(*mean.begin()), report);
	Classification(image, seg, nclass, bNoZero, &(*mean.begin()));
	return iter;
}

int TreeKMean(const ZImageBase& image, ZGrayByteImage& seg, int nclass, bool bNoZero, float* mean, float* stddev, float* prior, ostream* report)
{
	int iter = TreeKMean(image, nclass, bNoZero, mean, report);
	Classification(image, seg, nclass, bNoZero, mean);
	GetStatistics(image, seg, nclass, bNoZero, mean, stddev, prior);
	return iter;
}

int MRF(const ZImageBase& image, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, float* mean, float* stddev, float* prior, bool ICM, float B)
{
	if(nclass < 2)
	{
		ZError("MRF", "The number of classes must be more than 1!");
		return 0;
	}

	ZGrayFloatImage* pimage;
	if(isFloat(image)) pimage = (ZGrayFloatImage*)(&image); 
	else { pimage = new ZGrayFloatImage; image.SendPixelsTo(*pimage); }
	
	ZGrayFloatImage& fimage = *pimage;

	int width = image.Width(), height = image.Height(), depth = image.Depth(), slicesize = width*height;

	ZGrayByteImage mask(width, height, depth);
	seg.Create(width, height, depth);

	ZGrayFloatImage::const_iterator p_image(fimage);
	if(bNoZero)
	{
		ZGrayByteImage::iterator p(mask);
		for(; p_image.more(); p_image++, p++) if(*p_image == 0) *p = 1;
	}

	int i, c;

	vector<float*> prob(nclass);
	vector<float> variance(nclass);
	for(c=0; c<nclass; c++)
	{
		prob[c] = new float[image.Size()];
		variance[c] = stddev[c] * stddev[c];
	}

	ZGrayByteImage::iterator p_mask(mask), p_seg(seg);
	for(i=0, p_image.reset(); p_image.more(); p_image++, p_seg++, p_mask++, i++)
	{
		if(*p_mask) continue; 

		double	max = 1e20;
		int		clas=0;
		for(c=0; c<nclass; c++)
		{
			double val = *p_image - mean[c];
			val = 0.5 * val * val / variance[c] + log(stddev[c]);
			prob[c][i] = float(val);
			if(val < max) max = val, clas = c;
		}
		
		*p_seg = clas;
	}

	bool b3D = image.Depth() > 10;

	vector<float> priorprob(nclass);
	vector<BYTE> neighbour(b3D?27:9);
	int	offset3d = - slicesize - width - 1;
	int offset2d = - width - 1;

	float diaB = 0.8f, ddiaB = 0.6f;
	float TB2D = 4 + diaB*4, TB3D = TB2D + 2 + ddiaB*8;

	int iter;
	float BB = B > 0? B : (ICM ? 0.2f : 0.5);

//	ZProgressIndicator progress("progress", 70, NITER);
	for(iter=0; iter<NITER; iter++)
	{
		ZGrayByteImage::const_iterator p_mask(mask);
		ZGrayByteImage::iterator p_seg(seg);
		
		int pos=0, change = 0;
		for (; p_seg.more(); p_mask++, p_seg++, pos++)
		{
			if(*p_mask) continue;

			PBYTE	neighbour1 = &(*neighbour.begin()), src;

			if(b3D)
			{
				src = p_seg + offset3d;
				for(int k=0; k<3; k++, src+=slicesize)
				{
					PBYTE ptr = src;
					for(int j=0; j<3; j++, ptr+=width)
					for(int i=0; i<3; i++, neighbour1++)
					{
						if(ptr+i < seg.GetBuffer() || ptr+i >= seg.GetLimit())
							*neighbour1 = 0;
						else 
							*neighbour1 = ptr[i];
					}
				}
			}
			else 
			{
				PBYTE ptr = p_seg + offset2d;
				for(int j=0; j<3; j++, ptr+=width)
				for(int i=0; i<3; i++, neighbour1++)
				{
					if(ptr+i < seg.GetBuffer() || ptr+i >= seg.GetLimit())
						*neighbour1 = 0;
					else *neighbour1 = ptr[i]; 
				}
			}

			fill(priorprob.begin(), priorprob.end(), 0);

			PBYTE center = b3D ? &neighbour[13] : &neighbour[4];
			priorprob[center[-1]] ++;
			priorprob[center[1]] ++;
			priorprob[center[-3]] ++;
			priorprob[center[3]] ++;

			priorprob[center[-4]] += diaB;
			priorprob[center[4]] += diaB;
			priorprob[center[-2]] += diaB;
			priorprob[center[2]] += diaB;

			if(b3D)
			{
				priorprob[center[-9]] ++;
				priorprob[center[9]] ++;

				priorprob[center[-12]] += ddiaB;
				priorprob[center[12]] += ddiaB;
				priorprob[center[-10]] += ddiaB;
				priorprob[center[10]] += ddiaB;
				priorprob[center[-8]] += ddiaB;
				priorprob[center[8]] += ddiaB;
				priorprob[center[-6]] += ddiaB;
				priorprob[center[6]] += ddiaB;
			}

			for(c=0; c<nclass; c++) priorprob[c] = (b3D?TB3D:TB2D) - priorprob[c];

			int	 clas = 0;
			if(ICM)
			{
				float min = 1e10;
				for(c=0; c<nclass; c++)
				{
					float postProb = prob[c][pos] + BB * priorprob[c];
				
					if(min > postProb)
					{
						min = postProb;
						clas = c;
					}
				}
			}
			else
			{
				float temprature = float(3.0 / (log(2.0+iter)));

				static long	seed = -10;
				vector<double>  postProb(nclass);
				for(c=0; c<nclass; c++)
				{
					double tmp = - (BB * priorprob[c] + prob[c][pos]);
					postProb[c] = exp(tmp) / temprature;
				}

				for(c=1; c<nclass; c++) postProb[c] += postProb[c-1];

				double randprob = postProb[nclass-1] * rand1(seed);

				for(c=1; c<nclass; c++)
				{
					if(randprob > postProb[c-1] && randprob < postProb[c]) break;
				}
				
				if(c == nclass) c = 0;
				clas = BYTE(c);
			}

			if(*p_seg != BYTE(clas)) change++, *p_seg = BYTE(clas);
		}

		if(ICM)
		{
			if(change == 0) break;
			if(B < 0 && BB < 1) BB += 0.1f;
		}
		else
		{
			if(change < int(image.Size()/1000)) break;
			if(B<0) BB += 0.5f / NITER;
		}
	
//		progress.Step(iter);
	}
	if(bNoZero)
	{
		ZGrayByteImage::iterator p_mask(mask), p_seg(seg);
		for(; p_seg.more(); p_seg++, p_mask++) if(!*p_mask) (*p_seg)++;
	}
	for(c=0; c<nclass; c++) delete []prob[c];

	if(!isFloat(image)) delete pimage;
	return iter;
}

int HMRF_EM(const ZImageBase& image, ZGrayByteImage& seg, int nclass, bool bNoZero, float B, float* m, float* s, float* p, ostream* report)
{
	if(nclass < 2)
	{
		ZError("MRF_EM", "The number of classes must be more than 1!");
		return 0;
	}

	ZVector<float> mean(nclass), stddev(nclass), weight(nclass), variance(nclass), pprior(nclass);
	TreeKMean(image, seg, nclass, bNoZero, &(*mean.begin()), &(*stddev.begin()), &(*weight.begin()));

	ZGrayFloatImage* pimage;
	if(isFloat(image)) pimage = (ZGrayFloatImage*)(&image); 
	else { pimage = new ZGrayFloatImage; image.SendPixelsTo(*pimage); }
	
	ZGrayFloatImage& fimage = *pimage;

	int isize=0, width = image.Width(), height = image.Height(), depth = image.Depth(), slicesize = width*height, size=slicesize*depth;

	ZGrayByteImage mask(width, height, depth);

	ZGrayFloatImage::const_iterator p_image(fimage);
	if(bNoZero)
	{
		ZGrayByteImage::iterator p(mask);
		for(; p_image.more(); p_image++, p++) if(*p_image == 0) *p = 1; else isize++;
		seg -= 1;
	}
	else isize = size;

	int i, c, iter;

	vector<float*> post(nclass), prob(nclass);
	for(c=0; c<nclass; c++)
	{
		post[c] = new float[size];
		prob[c] = new float[size];
		variance[c] = stddev[c] * stddev[c];
	}

	bool b3D = depth>10;
	vector<BYTE> neighbour(b3D?27:9);
	int	offset3d = - slicesize - width - 1;
	int offset2d = - width - 1;

	float diaB = 0.8f, ddiaB = 0.6f;
	float TB2D = 4 + diaB*4, TB3D = TB2D + 2 + ddiaB*8;

	if(B==0) B=-1;
	float BB = B > 0? B : 0.1f;
	for(iter=0; iter<100; iter++)
	{
		ZGrayByteImage::iterator p_mask(mask), p_seg(seg);
		for(i=0, p_image.reset(); p_image.more(); p_image++, p_mask++, i++)
		{
			if(*p_mask) continue; 

			for(c=0; c<nclass; c++)
			{
				double val = *p_image - mean[c];
				val = 0.5 * val * val / variance[c] + log(stddev[c]);
				prob[c][i] = float(val);
			}
		}

		weight = 0;
		int change=0;
		for (i=0, p_seg.reset(), p_mask.reset(); p_seg.more(); p_mask++, p_seg++, i++)
		{
			if(*p_mask) continue;

			PBYTE	neighbour1 = &(*neighbour.begin()), src;

			if(b3D)
			{
				src = p_seg + offset3d;
				for(int k=0; k<3; k++, src+=slicesize)
				{
					PBYTE ptr = src;
					for(int j=0; j<3; j++, ptr+=width)
					for(int i=0; i<3; i++, neighbour1++)
					{
						if(ptr+i < seg.GetBuffer() || ptr+i >= seg.GetLimit())
							*neighbour1 = 0;
						else 
							*neighbour1 = ptr[i];
					}
				}
			}
			else 
			{
				PBYTE ptr = p_seg + offset2d;
				for(int j=0; j<3; j++, ptr+=width)
				for(int i=0; i<3; i++, neighbour1++)
				{
					if(ptr+i < seg.GetBuffer() || ptr+i >= seg.GetLimit())
						*neighbour1 = 0;
					else *neighbour1 = ptr[i]; 
				}
			}

			pprior = 0;

			PBYTE center = b3D ? &neighbour[13] : &neighbour[4];
			pprior[center[-1]] ++;
			pprior[center[1]] ++;
			pprior[center[-3]] ++;
			pprior[center[3]] ++;

			pprior[center[-4]] += diaB;
			pprior[center[4]] += diaB;
			pprior[center[-2]] += diaB;
			pprior[center[2]] += diaB;

			if(b3D)
			{
				pprior[center[-9]] ++;
				pprior[center[9]] ++;

				pprior[center[-12]] += ddiaB;
				pprior[center[12]] += ddiaB;
				pprior[center[-10]] += ddiaB;
				pprior[center[10]] += ddiaB;
				pprior[center[-8]] += ddiaB;
				pprior[center[8]] += ddiaB;
				pprior[center[-6]] += ddiaB;
				pprior[center[6]] += ddiaB;
			}

			for(c=0; c<nclass; c++) 
			{
				if(!pprior[c]) pprior[c] = -6;
				pprior[c] = (b3D?TB3D:TB2D) - pprior[c];
			}

			int	 clas = 0;

			float min = 1e10;
			for(c=0; c<nclass; c++)
			{
				post[c][i] = prob[c][i] + BB * pprior[c];
				prob[c][i] = BB * pprior[c];
			
				if(min > post[c][i])
				{
					min = post[c][i];
					clas = c;
				}
			}

			weight[clas]++;
			if(*p_seg != clas) change++, *p_seg = BYTE(clas);
		}

		weight /= isize;

		if(report)
		{
			*report  << setprecision (4) << setiosflags (ios::fixed)
				<< "iter=" << iter << "\t" ;
			for(c=0; c<nclass; c++)
				*report	<< "U" << c << "=" << mean[c] << "\t" 
					<< "S" << c << "=" << stddev[c] << "\t"
					<< "P" << c << "=" << weight[c] << "\t";
					
			*report << "change=" << change << endl;
		}
		if(change==0) break;

		if(B < 0 && BB < -B) BB += 0.05f;

		mean = 0; stddev = 0; variance = 0; 
		vector<float> nor(nclass);

		for(i=0, p_image.reset(), p_mask.reset(); p_image.more(); p_image++, p_mask++, i++)
		{
			if(*p_mask) continue;

			double norm = 0;
			for(c=0; c<nclass; c++) norm += post[c][i] = float(exp(-post[c][i]));
			if(norm < 1e-10) { *p_mask = true;  continue; }
			for(c=0; c<nclass; c++) post[c][i] /= norm;

			for (c = 0; c < nclass; c++)
			{
				mean[c] += post[c][i] * *p_image;
				nor[c] += post[c][i];
			}
		}

		for(c=0; c<nclass; c++) if(nor[c] > 1) mean[c] /= nor[c];

		for(i=0, p_image.reset(), p_mask.reset(); p_image.more(); p_image++, p_mask++, i++)
		{
			if(*p_mask) continue;

			for (c = 0; c < nclass; c++)
			{
				float val = *p_image - mean[c];
				variance[c] += post[c][i] * val * val;
			}
		}

		for(c=0; c<nclass; c++)
		{
			if(nor[c] > 1) 
			{
				variance[c] /= nor[c];
				stddev[c] = float(sqrt(variance[c]));
			}
		}

		for(i=0, p_image.reset(), p_mask.reset(); p_image.more(); p_image++, p_mask++, i++)
		{
			if(*p_mask) continue; 

			for(c=0; c<nclass; c++)
			{
				double val = *p_image - mean[c];
				val = 0.5 * val * val / variance[c] + log(stddev[c]);
				post[c][i] = prob[c][i] + val;
			}
		}
	}
	if(m!=0) for(int c=0; c<nclass; c++) m[c] = mean[c];
	if(s!=0) for(int c=0; c<nclass; c++) s[c] = stddev[c];
	if(p!=0) for(int c=0; c<nclass; c++) p[c] = weight[c];

	for(c=0; c<nclass; c++) delete [](prob[c]), delete [](post[c]);

	if(bNoZero)
	{
		ZGrayByteImage::iterator p_mask(mask), p_seg(seg);
		for(; p_seg.more(); p_seg++, p_mask++) if(!*p_mask) (*p_seg)++;
	}
	if(!isFloat(image)) delete pimage;
	return iter;
}
