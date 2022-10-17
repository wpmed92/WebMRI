/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declare the IO class for the AVW format.

#ifndef __AVWIMAGE_H
#define __AVWIMAGE_H

//AVW (A Visualization Workshop) is a comprehensive library of imaging functions 
//which permits full exploration and analysis of multidimensional and multimodal 
//biomedical image data sets. 

#include "format.h"

#define AVW_READ  0
#define AVW_WRITE 1

struct header_key          		//      header_key
{                               // off + size
    int sizeof_hdr;				//Must indicate the byte size of the header file
    char data_type[10];         
    char db_name[18];           // 14 + 18
    int extents;				//Should be 16384, the image file is created as contiguous with a minimum extent size.
    short int session_error;    // 36 + 2
    char regular;				//Must be `r' to indicate that all images and volumes are the same size. 
    char hkey_un0;
};           					// total=40

struct image_dimension       
{
    short int dim[8];		//array of the image dimensions
			//dim[0]      Number of dimensions in database; usually 4 
			//dim[1]      Image X dimension; number of pixels in an image row 
			//dim[2]      Image Y dimension; number of pixel rows in slice 
			//dim[3]      Volume Z dimension; number of slices in a volume 
			//cdim[4]     Time points, number of volumes in database.
    char vox_units[4];		//specifies the spatial units of measure for a voxel
    char cal_units[8];		//specifies the name of the calibration unit 
    short int unused1;			/* 24 + 2    */
    short int datatype;		//datatype for this image set 
			/*Acceptable values for datatype are
			#define DT_NONE				0 
			#define DT_UNKNOWN			0      // Unknown data type 
			#define DT_BINARY			1      //Binary (1 bit per voxel) 
			#define DT_UNSIGNED_CHAR	2      //Unsigned character (8 bits per voxel) 
			#define DT_SIGNED_SHORT		4      //Signed short (16 bits per voxel) 
			#define DT_SIGNED_INT		8      //Signed integer (32 bits per voxel) 
			#define DT_FLOAT			16     //Floating point (32 bits per voxel) 
			#define DT_COMPLEX			32     //Complex (64 bits per voxel; 2 floating point numbers) 
			#define DT_DOUBLE			64     //Double precision (64 bits per voxel) 
			#define DT_RGB				128
			#define DT_ALL				255
			*/

    short int bitpix;		//number of bits per pixel; 1, 8, 16, 32, or 64
    short int dim_un0;		//unused 
    float pixdim[8];		//Parallel array to dim[], giving real world measurements in mm. and ms. 
        	/* 
        		pixdim[] specifies the voxel dimensions:
        		pixdim[1] - voxel width in mm
        		pixdim[2] - voxel height in mm
        		pixdim[3] - slice thickness in mm
        	*/
    float vox_offset;		//byte offset in the .img file at which voxels start. This value can be 
                            //negative to specify that the absolute value is applied for every image
    float funused1;                       	/* 72 + 4    */
    float funused2;                      	/* 76 + 4    */
    float funused3;                      	/* 80 + 4    */
    float cal_max;			//specify the range of calibration values 
    float cal_min;
    int compressed;                     	/* 92 + 4    */
    int verified;                     	/* 96 + 4    */
    int glmax, glmin;		//The maximum and minimum pixel values for the entire database
};          				/* total=108 */


///////////////////////////////////////////////////////////
//The data_history substructure is not required, but the orient field is used to 
//indicate individual slice orientation and determines whether the Movie program 
//will attempt to flip the images before displaying a movie sequence.
//
//orient:       slice orientation for this dataset. 
//0      transverse unflipped 
//1      coronal unflipped 
//2      sagittal unflipped 
//3      transverse flipped 
//4      coronal flipped 
//5      sagittal flipped
///////////////////////////////////////////////////////////
struct data_history          		/*      data_history     */
{                                    	/* off + size*/
    char descrip[80];                	/* 0 + 80    */
    char aux_file[24];               	/* 80 + 24   */
    char orient;                     	/* 104 + 1   */
    char originator[10];             	/* 105 + 10  */
    char generated[10];              	/* 115 + 10  */
    char scannum[10];                	/* 125 + 10  */
    char patient_id[10];             	/* 135 + 10  */
    char exp_date[10];               	/* 145 + 10  */
    char exp_time[10];               	/* 155 + 10  */
    char hist_un0[3];                	/* 165 + 3   */
    int views;                       	/* 168 + 4   */
    int vols_added;                  	/* 172 + 4   */
    int start_field;                 	/* 176 + 4   */
    int field_skip;                  	/* 180 + 4   */
    int omax,omin;                   	/* 184 + 8   */
    int smax,smin;                   	/* 192 + 8   */
};                     			/* total=200 */

struct dsr                 		/*      dsr              */
{                                  		/* off + size*/
    struct header_key hk;          		/* 0 + 40    */
    struct image_dimension dime;   		/* 40 + 108  */
    struct data_history hist;      		/* 148 + 200 */
};                     			/* total=348 */
        
/* Acceptable values for hdr.dime.datatype */

#define DT_NONE				0
#define DT_UNKNOWN			0
#define DT_BINARY			1
#define DT_UNSIGNED_CHAR		2
#define DT_SIGNED_SHORT			4
#define DT_SIGNED_INT			8
#define DT_FLOAT			16
#define DT_COMPLEX			32
#define DT_DOUBLE			64
#define DT_RGB				128
#define DT_ALL				255

class ZAVWImage  : public ZImageFile
{
protected:
	UINT		m_nVolsize, m_nVoxelbytes;
	UINT		m_nVolume;

	dsr			m_Header;

	bool		m_bHeader;		//header already read

	std::fstream	m_HeaderStream;

	bool ReadHeader();
	void InitHeader();
	void ReverseDsr();
	void ReverseBuffer(void *buffer, int no_elements);
	
	bool ReadVolumes(PBYTE buffer);
	bool WriteVolumes(const ZImageBase&, bool intensityscale);
	void SeekVolume(UINT vols);

	void SetDim(int x, int y, int z);
	void SetVoxDim(float x, float y, float z, float tr);
	void SetMinMax(int min, int max);
	void SetVoxUnits(const char *units);
	void SetOrientation(const char orient);
	void SetOriginator(const short orig[5]);

public:
	ZAVWImage() : ZImageFile(), m_bHeader(false) {}
	~ZAVWImage() {Close();}
	bool Open(const char *fname, bool read);
	void Close();
	bool IsOpen() {return m_ImageStream.is_open()!=0 && m_HeaderStream.is_open()!=0; }
	void SetVolume(int nVol) { m_nVolume = nVol;}

	ZImageBase* ReadFile();
	bool WriteFile(const ZImageBase& image, bool intensityscale=false);
};

#endif
