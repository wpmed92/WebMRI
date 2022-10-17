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
#include "vclassification.h"

using namespace std;

void Reform(const std::vector<ZImageBase*>& images, ZMatrix<float>& mat)
{
	if(images.size() == 0)
	{
		ZError("Reform", "Images is empty!");
		return;
	}

	UINT size = images[0]->Size(), spec=images.size();
	for(UINT i=1; i<spec; i++)
	{
		if(size != images[i]->Size())
		{
			ZError("Reform", "Images have different sizes!");
			return;
		}
	}

	mat.resize(size, spec, 0);

	UINT width=images[0]->Width(), height=images[0]->Height(), depth=images[0]->Depth();
	ZVector<float> v(spec);
	ZMatrix<float>::row_iterator p_mat = mat.pbegin();
	for(UINT k=0; k<depth; k++)
	for(UINT j=0; j<height; j++)
	for(UINT i=0; i<width; i++, p_mat++)
	{
		bool zero = false;
		for(UINT m=0; m<spec; m++) 
			if((v[m] = images[m]->GetIntensity(i, j, k)) == 0) { zero=true; break; }

		if(!zero) *p_mat = v;
	}
}

void Classification(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, ZMatrix<float>& mean, vector< ZMatrix<float> >* covariance, float* prior)
{
	if(nclass < 2)
	{
		ZError("Classification", "The number of classes must be more than 1!");
		return;
	}

	if(seg.Size() != images.NRows())
	{
		ZError("Classification", "The size of the seg is different from that of images!");
		return;
	}

	int spec = images.NCols();

	ZImage<BYTE>::iterator p_seg(seg);
	ZMatrix<float>::const_row_iterator p_image = images.pbegin();
	if(covariance==0)
	{
		for(; p_image!=images.pend(); p_image++, p_seg++)
		{
			if(!bNoZero || (*p_image)[0])
			{
				int cc = 0; float max = 1e10;
				for(int c=0; c<nclass; c++)
				{
					float e = sqrdis(*p_image, mean[c]);
					if(max > e) max = e, cc=c;
				}
				*p_seg = bNoZero ? cc+1 : cc;
			}
		}
		return;
	}

	vector< ZMatrix<float> > icovariance(nclass);
	ZVector<float> stddev(nclass);
	for(int c=0; c<nclass; c++) 
	{
		icovariance[c] = inv((*covariance)[c]);
		stddev[c] = sqrt(1.0f/ Abs(det((*covariance)[c])));
	}

	ZVector<float> dif(spec);
	for(; p_image!=images.pend(); p_image++, p_seg++)
	{
		if(!bNoZero || (*p_image)[0])
		{
			double max=0; int cc = 0;
			for(int c=0; c<nclass; c++)
			{
				dif = *p_image;  dif -= mean[c];
				double v = 0;
				for(int i=0; i<spec; i++)
				{
					double u = 0;
					for(int j=0; j<spec; j++) u += dif[j] * icovariance[c](j, i);
					v += u * dif[i];
				}
				double prob = exp(-v/2) * stddev[c] * M_1_SQRT_2PI;
				if(prior) prob *= prior[c];
				if(max < prob) max = prob, cc = c;
			}
			*p_seg = bNoZero ? cc+1 : cc;
		}
	}
}

void GetStatistics(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, 
					 ZMatrix<float>& mean, vector< ZMatrix<float> >& covariance, float* prior)
{
	if(seg.Size() != images.NRows()) return;

	int c, spec = images.NCols();
	for(c=0; c<nclass; c++) covariance[c] = 0, prior[c] = 0;
	mean = 0; 
	ZGrayByteImage::const_iterator p_seg(seg);
	ZMatrix<float>::const_row_iterator p_image = images.pbegin();
	for(; p_image!=images.pend(); p_image++, p_seg++)
	{
		int c = *p_seg;
		if(bNoZero) if(!c) continue; else c--;

		mean[c] += *p_image;
		prior[c] ++;
	}

	float norm = 0;
	for(c=0; c<nclass; c++) { norm += prior[c]; if(prior[c] > 0) mean[c] /= prior[c]; }

	ZVector<float> dif(spec);
	p_seg.reset(); p_image = images.pbegin();
	for(; p_image!=images.pend(); p_image++, p_seg++)
	{
		int c = *p_seg;
		if(bNoZero) if(!c) continue; else c--;

		dif = *p_image; dif -= mean[c];

		for(int i=0; i<spec; i++)
		for(int j=0; j<spec; j++)
			covariance[c](i,j) += dif[i] * dif[j];
	}

	for(c=0; c<nclass; c++) 
	{
		if(prior[c] > 0) covariance[c] *= 1.0f / prior[c];
		prior[c] /= norm;
	}
}

int EM(const ZMatrix<float>& images, int nclass, int NITER, bool bNoZero, ZMatrix<float>& mean, 
	   vector< ZMatrix<float> >& covariance, float* prior)
{
	if(nclass < 2)
	{
		ZError("EM", "The number of classes must be more than 1!");
		return 0;
	}

	int spec = images.NCols();
	float *weight = new float[nclass];

	int c, iter, size = 0;
	for(c=0; c<nclass; c++)	weight[c] = prior ? prior[c] : 1.0f / nclass;

	ZMatrix<float>::const_row_iterator p_image = images.pbegin();
	for(; p_image!=images.pend(); p_image++)
	{
		if(!bNoZero || (*p_image)[0]) size ++;
	}

	vector< ZMatrix<float> > icovariance(nclass);
	ZVector<float> stddev(nclass);
	for(c=0; c<nclass; c++) 
	{
		icovariance[c] = inv(covariance[c]);
		stddev[c] = sqrt(1.0f/ Abs(det(covariance[c])));
	}

	vector<float*> prob(nclass);
	for(c=0; c<nclass; c++) prob[c] = new float[size];

	double error = 1e20;
//	ZProgressIndicator progress("progress", 70, NITER);
	ZVector<float> dif(spec), nor(nclass);
	for(iter=0; iter<NITER; iter++)
	{
		double tmp=0; 

		int c, pos=0;
		for(p_image=images.pbegin(); p_image!=images.pend(); p_image++)
		{
			if(bNoZero && !(*p_image)[0]) continue;
			double norm = 0;
			for(c=0; c<nclass; c++)
			{
				double v = 0; 
				dif = *p_image; dif -= mean[c];
				for(int i=0; i<spec; i++)
				{
					double u = 0;
					for(int j=0; j<spec; j++) u += dif[j] * icovariance[c](j, i);
					v += u * dif[i];
				}

				prob[c][pos] = float(exp(-v/2) * stddev[c] * M_1_SQRT_2PI * weight[c]);
				norm += prob[c][pos]; 
			}
			tmp -= log(norm);
			for (c = 0; c < nclass; c++) prob[c][pos] /= norm;
			pos++;
		}

		tmp /= size;
		if(Abs(error-tmp) < 1e-10) break; error = tmp;

		nor = 0; mean = 0; stddev = 0; 
		for(c=0; c<nclass; c++) covariance[c] = 0;

		pos=0;
		for(p_image=images.pbegin(); p_image!=images.pend(); p_image++)
		{
			if(bNoZero && !(*p_image)[0]) continue;
			for (c = 0; c < nclass; c++)
			{
				float ppp = prob[c][pos];
				nor[c] += ppp;
				for(int i=0; i<spec; i++) mean(c, i) += (*p_image)[i] * ppp;
			}
			pos++;
		}

		for(c=0; c<nclass; c++) if(nor[c] > 1) mean[c] /= nor[c];

		pos=0;
		for(p_image=images.pbegin(); p_image!=images.pend(); p_image++)
		{
			if(bNoZero && !(*p_image)[0]) continue;
			for (c = 0; c < nclass; c++)
			{
				dif = *p_image; dif -= mean[c];
				float ppp = prob[c][pos];
				for(int i=0; i<spec; i++)
				for(int j=0; j<spec; j++)
					covariance[c](i,j) += dif[i] * dif[j] * ppp;
			}
			pos++;
		}

		for(c=0; c<nclass; c++)
		{
			if(nor[c] > 1) covariance[c] *= 1.0f / nor[c];
			icovariance[c] = inv(covariance[c]);
			stddev[c] = sqrt(1.0f / Abs(det(covariance[c])));
			weight[c] = nor[c] / size;
		}

//		progress.Step(iter);
	}

	for(c=0; c<nclass; c++) { if(prior) prior[c] = weight[c]; } 
	delete []weight;

	return iter;
}

int EM(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, 
	   ZMatrix<float>& mean,  vector< ZMatrix<float> >& covariance, float* prior)
{
	int iter = EM(images, nclass, NITER, bNoZero, mean, covariance, prior);
	Classification(images, seg, nclass, bNoZero, mean, &covariance, prior);
	return iter;
}

int TreeEM(const ZMatrix<float>& images, int nclass, int NITER, bool bNoZero, 
		   ZMatrix<float>& mean,  vector< ZMatrix<float> >& covariance, float* prior, ostream* report)
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

	int spec = images.NCols(), i, c;

	ZMatrix<float> m(nclass, spec);
	vector< ZMatrix<float> > v(nclass); for(c=0; c<nclass; c++) v[c].resize(spec, spec);
	ZVector<float> p(nclass), dif(spec);

	images.RowStatistics(m[0], v[0], bNoZero);

	ZVector<float> mm, eigenvalues, stddev;
	ZMatrix<float> eigenvectors, transm;

	if(report) *report << "starting with 1 class -- mean: " << m[0] << endl;

	int totaliter=0;

	int splitnode = 0;
	const ZMatrix<float> *pixels = &images;
	vector< ZMatrix<float> > simage;
	for(int nc=1; nc<nclass; nc++)
	{
		if(Eigensystem(*pixels, mm, eigenvectors, eigenvalues) == false)
		{
			ZError("TreeEM", "Fail to compute eigen system!");
			return -1;
		}
		transm = trans(eigenvectors);
		mm = transm * mm;

		stddev = sqrt(eigenvalues);

		for(i=nc-1; i>splitnode; i--) 
		{
			m[i+1] = m[i];
			v[i+1] = v[i] * nc / (nc+1);
		}
		m[splitnode] = eigenvectors * (mm - stddev); 
		v[splitnode] /= 3;
		m[splitnode+1] = eigenvectors * (mm + stddev);
		v[splitnode+1] /= 3;
		for(i=0; i<nc+1; i++) p[i] = 1.0f / (nc+1);
		if(nc+1 == nclass) NITER *= 2;
		totaliter += EM(images, nc+1, NITER, bNoZero, m, v, p.pbegin());

		if(report) 
		{
			*report << "splite node " << splitnode << " total " << nc+1 << " classes" << endl;
			*report << "mean of each classes:\t";
			for(i=0; i<nc; i++)
				*report << m[i] << endl;
		}
		
		ZVector<float> difs(nc+1);
		ZVector<int> num(nc+1);
		ZMatrix<float>::const_row_iterator p_image = images.pbegin();
		for(; p_image!=images.pend(); p_image++)
		{
			if(bNoZero && !(*p_image)[0]) continue;

			dif = *p_image; double max = 1e10;
			int cc = 0;
			for(c=0; c<=nc; c++) 
			{
				double v = sqrdis(dif, m[c]);
				if(v < max) max = v, cc = c;
			}
			difs[cc] += max; num[cc] ++; simage[cc].push_back(*p_image);
		}

		for(c=0; c<=nc; c++) if(num[c]) difs[c] /= num[c];
		splitnode = max_element(difs.begin(), difs.end())-difs.begin();
		pixels = &simage[splitnode];
	}

	mean = m;
	for(c=0; c<nclass; c++) covariance[c] = v[c];
	for(c=0; c<nclass; c++) prior[c] = p[c];

	return totaliter;
}

int TreeEM(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, 
		   ZMatrix<float>& mean, vector< ZMatrix<float> >& covariance, float* prior, ostream* report)
{
	int iter = TreeEM(images, nclass, NITER, bNoZero, mean, covariance, prior, report);
	Classification(images, seg, nclass, bNoZero, mean, &covariance, prior);
	return iter;
}

int KMean(const ZMatrix<float>& images, int nclass, int NITER, bool bNoZero, ZMatrix<float>& mean, ostream* report)
{
	if(nclass < 2)
	{
		ZError("KMean", "The number of classes must be more than 1!");
		return 0;
	}

	int iter, spec = images.NCols(), size=0;
	ZMatrix<float>::const_row_iterator p_image = images.pbegin();
	for(; p_image!=images.pend(); p_image++)
	{
		if(!bNoZero || (*p_image)[0]) size ++;
	}

	if(report)
	{
		*report  << setprecision (4) << setiosflags (ios::fixed)
			<< "initial settings:\tU" << 0 << "=" << mean[0] << "\n";
		for(int c=1; c<nclass; c++)
			*report << "\t\t\tU" << c << "=" << mean[c] << "\n";
		*report << endl; 
	}

	ZVector<UINT> num(nclass);
	ZMatrix<float> nmean(nclass, spec);
	ZVector<float> dif(spec);

	for(iter=0; iter<NITER; iter++)
	{
		int c;

		num = 0, nmean = 0;
		for(p_image = images.pbegin(); p_image!=images.pend(); p_image++)
		{
			if(bNoZero && !(*p_image)[0]) continue;

			float min = 1e10; int cc = 0;
			for(c=0; c<nclass; c++)
			{
				float val = sqrdis(*p_image, mean[c]);
				if(min > val) min = val, cc = c;
			}

			num[cc] ++; nmean[cc] += *p_image;
		}

		bool quit=true;
		for(c=0; c<nclass; c++) 
		{
			if(num[c]==0 || (mean[c] == (nmean[c]/=num[c]))) continue;
			quit = false;
			mean[c] = nmean[c];
		}

		if(report)
		{
			*report  << setprecision (4) << setiosflags (ios::fixed)
				<< "iter=" << iter << "\t" ;
			*report << "U" << 0 << "=" << mean[0] << "\t" 
				<< "P" << 0 << "=" << float(num[0])/size << "\n";
			for(c=1; c<nclass; c++)
				*report << "\tU" << c << "=" << mean[c] << "\t" 
					<< "P" << c << "=" << float(num[c])/size << "\n";
		}
		if(quit) break;
	}

	return iter;
}

int KMean(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, int NITER, bool bNoZero, 
		  ZMatrix<float>& mean, vector< ZMatrix<float> >& covariance, float* prior, ostream* report)
{
	int iter = KMean(images, nclass, NITER, bNoZero, mean, report);
	Classification(images, seg, nclass, bNoZero, mean);
	GetStatistics(images, seg, nclass, bNoZero, mean, covariance, prior);
	return iter;
}

int TreeKMean(const ZMatrix<float>& images, int nclass, bool bNoZero, ZMatrix<float>& mean, ostream* report)
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

	int spec = images.NCols(), i, c, size=0;

	ZMatrix<float> m(nclass, spec);

	ZMatrix<float>::const_row_iterator p_image = images.pbegin();
	for(; p_image!=images.pend(); p_image++)
	{
		if(!bNoZero || (*p_image)[0]) size++, m[0] += *p_image;
	}
	m[0] /= size;

	ZVector<float> mm, eigenvalues, stddev, dif(spec);
	ZMatrix<float> eigenvectors, transm;

	if(report) *report << "starting with 1 class -- mean: " << m[0] << endl;

	int totaliter=0, splitnode = 0;
	vector< ZMatrix<float> > simage(nclass);

	for(p_image = images.pbegin(); p_image!=images.pend(); p_image++)
	{
		if(bNoZero && !(*p_image)[0]) continue;

		simage[0].push_back(*p_image);
	}

	const ZMatrix<float> *pixels = &simage[0];
	for(int nc=1; nc<nclass; nc++)
	{
		if(Eigensystem(*pixels, mm, eigenvectors, eigenvalues) == false)
		{
			ZError("TreeEM", "Fail to compute eigen system!");
			return -1;
		}
		transm = trans(eigenvectors);
		mm = transm * mm;

		stddev = sqrt(eigenvalues);

		for(i=nc-1; i>splitnode; i--) m[i+1] = m[i];
		m[splitnode] = eigenvectors * (mm - stddev); 
		m[splitnode+1] = eigenvectors * (mm + stddev);

		if(report) 
		{
			*report << "splite node " << splitnode << " total " << nc+1 << " classes" << endl;
			*report << "mean of each classes:\t"; *report << m[0] << endl;
			for(i=1; i<=nc; i++) *report << "\t\t\t" << m[i] << endl;
			*report << "Total number of iterations so far: " << totaliter << endl << endl;
		}
		
		totaliter += KMean(images, nc+1, 100, bNoZero, m);

		ZVector<int> num(nc+1);
		for(c=0; c<nclass; c++) simage[c].clear();
		for(p_image = images.pbegin(); p_image!=images.pend(); p_image++)
		{
			if(bNoZero && !(*p_image)[0]) continue;

			dif = *p_image; double max = 1e10;
			int cc = 0;
			for(c=0; c<=nc; c++) 
			{
				double v = sqrdis(dif, m[c]);
				if(v < max) max = v, cc = c;
			}
			num[cc] ++; simage[cc].push_back(*p_image);
		}

		splitnode = max_element(num.pbegin(), num.pend())-num.pbegin();
		pixels = &simage[splitnode];
	}

	mean = m;

	return totaliter;
}

int TreeKMean(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, ostream* report)
{
	ZMatrix<float> mean;
	int iter = TreeKMean(images, nclass, bNoZero, mean, report);
	Classification(images, seg, nclass, bNoZero, mean);
	return iter;
}

int TreeKMean(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, ZMatrix<float>& mean, vector< ZMatrix<float> >& covariance, float* prior, ostream* report)
{
	int iter = TreeKMean(images, nclass, bNoZero, mean, report);
	Classification(images, seg, nclass, bNoZero, mean);
	GetStatistics(images, seg, nclass, bNoZero, mean, covariance, prior);
	return iter;
}

int HMRF(const ZMatrix<float>& images, ZGrayByteImage& seg, int nclass, bool bNoZero, ZMatrix<float>* m, vector< ZMatrix<float> >* cov, float* p, ostream* report)
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
/*
	int spec = images.NCols(), i, size=0;

	int m_nWidth = images[0]->Width();		// width of the new image with zeor boundaries removed
	int m_nHeight = images[0]->Height();		// height of the new image with zeor boundaries removed
	int m_nSliceSize = m_nWidth * m_nHeight;// slicesize of the new image with zeor boundaries removed
	int m_nDepth = images[0]->Depth();		// depth of the new image with zeor boundaries removed
	int m_nSize = m_nDepth * m_nSliceSize;	// size of the new image with zeor boundaries removed
	bool b3D = m_nDepth > 10;
	int offset3d = - m_nSliceSize - m_nWidth - 1;
	int offset2d = - m_nWidth - 1;
	
	float diaB = 0.8f, ddiaB = 0.6f;
	float TB2D = 4 + diaB*4, TB3D = TB2D + 2 + ddiaB*8;

	if(B==0) B=-1;
	float BB = B > 0? B : 0.1f;

//	TreeKMean(images, nclass, bNoZero, mean, report);
	ZMatrix<float> prob(m_nSize, nclass), post(m_nSize, nclass);
*/

	return 0;
}
