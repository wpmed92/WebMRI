/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "filter.h"
using namespace std;

#define splitchannel \
	ZGrayFloatImage* res[3]; int nchannel=1; \
	if(image.isColor()) \
	{ nchannel=3; ZImageBase* c[3]; \
		for(int i=0; i<3; i++) c[i] = GrayTypeCopyFrom(&image); \
		SplitRGBChannel(image, *c[0], *c[1], *c[2]); \
		for(int j=0; j<3; j++) { res[j]=new ZGrayFloatImage; c[j]->CopyPixelsTo(*res[j]); delete c[j]; } \
	} else { res[0]=new ZGrayFloatImage; image.CopyPixelsTo(*res[0]); }

#define combinechannel \
	ZImageBase* ret=res[0]; \
	if(image.isColor()) { ZImageBase* r = new ZRGBFloatImage; \
		CombineRGBChannel(*r, *res[0], *res[1], *res[2]); \
		delete res[0]; delete res[1]; delete res[2]; ret = r; }

ZImageBase*	Gradient(const ZImageBase& image)
{
	splitchannel
	ZGrayFloatImage grad;
	for(int c=0; c<nchannel; c++)
	{
		grad = *res[c];  grad.Inflate(1, 0, 1, 0, false, 0);
		UINT width = grad.Width();
		grad.SetROI(1, image.Width(), 1, image.Height(), 0, grad.Depth()-1);
		ZImage<float>::iterator in(grad), out(*res[c]);
		for(;in.more(); in++, out++)
		{
			float  xgrad = *in - *(in-1), ygrad = *in - *(in-width);
			*out = float(sqrt(xgrad * xgrad + ygrad * ygrad));
		}
	}
	combinechannel
	return ret;
}

ZImageBase*	XGradient(const ZImageBase& image)
{
	splitchannel
	ZGrayFloatImage grad;
	for(int c=0; c<nchannel; c++)
	{
		grad = *res[c];  grad.Inflate(1, 0, 0, 0, false, 0);
		grad.SetROI(1, image.Width(), 0, image.Height()-1, 0, image.Depth()-1);
		ZImage<float>::iterator in(grad), out(*res[c]);
		for(;in.more(); in++, out++) *out = *in - *(in-1);
	}
	combinechannel
	return ret;
}

ZImageBase*	YGradient(const ZImageBase& image)
{
	splitchannel
	ZGrayFloatImage grad;  int width = image.Width();
	for(int c=0; c<nchannel; c++)
	{
		grad = *res[c];  grad.Inflate(0, 0, 1, 0, false, 0);
		grad.SetROI(0, width-1, 1, image.Height(), 0, image.Depth()-1);
		ZImage<float>::iterator in(grad), out(*res[c]);
		for(;in.more(); in++, out++) *out = *in - *(in-width);
	}
	combinechannel
	return ret;
}

ZImageBase*	Gradient3D(const ZImageBase& image)
{
	splitchannel
	ZGrayFloatImage grad;
	for(int c=0; c<nchannel; c++)
	{
		grad = *res[c];  grad.Inflate(1, 0, 1, 0, 1, 0, false, 0);
		UINT width = grad.Width(), height = grad.Height(), xysize = width*height;
		grad.SetROI(1, image.Width(), 1, image.Height(), 1, image.Depth());

		ZImage<float>::iterator in(grad), out(*res[c]);

		for(;in.more(); in++, out++)
		{
			float xgrad = *in - *(in-1);
			float ygrad = *in - *(in-width);
			float zgrad = *in - *(in-xysize);
			*out = float(sqrt(xgrad * xgrad + ygrad * ygrad + zgrad * zgrad));
		}
	}
	combinechannel
	return ret;
}

void GaussianBlur(ZImageBase& image, float sigma, bool b2d)
{
	if(sigma < 0) sigma = 1;

	splitchannel
	for(int c=0; c<nchannel; c++)
	{
		float x, y, z;
		res[c]->GetPixelDim(x, y, z);
		GaussianBlur(*res[c], sigma/x, sigma/y, sigma/z, b2d);
	}
	combinechannel
	ZImageBase* pImage = TypeCopyFrom(&image);
	ret->SendPixelsTo(*pImage);
	image.Replace(*pImage);
	delete pImage;
	delete ret;
}

void GaussianBlur(ZGrayFloatImage& image, float x_sigma, float y_sigma, float z_sigma, bool b2d)
{
	int x_size = image.Width();
	int y_size = image.Height();
	int z_size = image.Depth();
	int row_size = image.BytesPerROIRow();
	int slice_size = image.PixelsPerSlice();
	if (x_sigma>0)
	{
		vector<float> array(x_size);

		int lp_mask_size = int(0.5+3*x_sigma);
		vector<float> lp_exp_orig(2*lp_mask_size+1);
		float* lp_exp = &(*lp_exp_orig.begin()) + lp_mask_size;

		for(int i=-lp_mask_size; i<=lp_mask_size; i++)
			lp_exp[i] = float(exp(-float(i*i)/ 2 / (x_sigma*x_sigma) ));

		float* p1 = image.ImageOrigin();
		for(int z=0; z<z_size; z++, p1+=slice_size)
		{
			float* p2 = p1;
			for(int y=0; y<y_size; y++, p2+=image.ImageWidth())
			{
				memcpy(&(*array.begin()), p2, row_size);
				float* p3 = p2;

				for(int x=0; x<x_size; x++)
				{
					float sum=0, total=0;
					int start=Max(x-lp_mask_size,0), stop=Min(x+lp_mask_size,x_size-1);

					for(int ii=start; ii<=stop; ii++)
					{
						sum += lp_exp[ii-x] * array[ii];
						total+=lp_exp[ii-x];
					}

					*p3++ = sum/total;
				}
			}
		}
	}

	if (y_sigma>0)
	{
		vector<float> array(y_size);

		int lp_mask_size = int(0.5+3*y_sigma);
		vector<float> lp_exp_orig(2*lp_mask_size+1);
		float* lp_exp = &(*lp_exp_orig.begin()) + lp_mask_size;

		for(int i=-lp_mask_size; i<=lp_mask_size; i++)
			lp_exp[i] = float(exp(-float(i*i)/ 2 / (y_sigma*y_sigma) ));

		float* p1 = image.ImageOrigin();
		for(int z=0; z<z_size; z++, p1+=slice_size)
		{
			float* p2 = p1;
			for(int x=0; x<x_size; x++, p2++)
			{
				float* p3 = p2;
				int y;
				for(y=0; y<y_size; y++, p3+=image.ImageWidth()) array[y] = *p3;
				p3 = p2;

				for(y=0; y<y_size; y++, p3+=image.ImageWidth())
				{
					float sum=0, total=0;
					int start = Max(y-lp_mask_size,0), stop=Min(y+lp_mask_size,y_size-1);

					for(int ii=start; ii<=stop; ii++)
					{
						sum += lp_exp[ii-y] * array[ii];
						total += lp_exp[ii-y];
					}

					*p3 = sum/total;
				}

			}
		}
	}

	if (!b2d && (z_sigma>0) && (z_size>1) )
	{
		vector<float> array(z_size);

		int lp_mask_size = int(0.5+3*z_sigma);
		vector<float> lp_exp_orig(2*lp_mask_size+1);
		float* lp_exp = &(*lp_exp_orig.begin()) + lp_mask_size;

		for(int i=-lp_mask_size; i<=lp_mask_size; i++)
			lp_exp[i] = float(exp(-float(i*i)/ 2 / (z_sigma*z_sigma) ));

		float* p1 = image.ImageOrigin();
		for(int y=0; y<y_size; y++, p1+=image.ImageWidth())
		{
			float* p2 = p1;
			for(int x=0; x<x_size; x++, p2++)
			{
				float* p3 = p2;
				int z;
				for(z=0; z<z_size; z++, p3+=slice_size) array[z] = *p3;
				p3 = p2;

				for(z=0; z<z_size; z++, p3+=slice_size)
				{
					float sum=0, total=0;
					int start = Max(z-lp_mask_size,0), stop = Min(z+lp_mask_size,z_size-1);

					for(int ii=start; ii<=stop; ii++)
					{
						sum+=lp_exp[ii-z] * array[ii];
						total+=lp_exp[ii-z];
					}

					*p3 = sum/total;
				}
			}
		}
	}
}
