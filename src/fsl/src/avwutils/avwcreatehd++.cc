//     avwcpgeom++.cc - Copy certain parts of an AVW header
//     Mark Jenkinson, Steve Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2001-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include <fstream>
#include <iostream>

using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: avwcreatehd <xsize> <ysize> <zsize> <tsize> <xvoxsize> <yvoxsize> <zvoxsize> <tr> <xorigin> <yorigin> <zorigin> <datatype> <headername>" << endl;
  cout << "       avwcreatehd <nifti_xml_file> <headername>" << endl;
  cout << "  In the second form, an XML-ish form of nifti header is read (as output by avwhd -x)" << endl;
  cout << "  Note that stdin is used if '-' is used in place of a filename" << endl;

}


int avwcreatehd_main(int argc, char *argv[])
{
  FSLIO* fslio;
  void *buffer=NULL;
  char *hdrxml, *filename;
  int fileread=1, filetype=-1, existingimage=0;
  size_t bufsize=0;
  short x,y,z,v, dt=-1;
  
  if (argc==3) {
      /* use the XML form of header specification */
      filename = argv[2];
  } else {
      filename = argv[13];
  }

  /* check if file already exists and if so, read the image contents */
  /* also store the size of this buffer for later in case it is wrong */ 
  if (FslFileExists(filename)) {
    /* buffer = FslReadAllVolumes(fslio,filename); */
    existingimage = 1;
    fslio = FslOpen(filename,"rb");
    FslGetDim(fslio,&x,&y,&z,&v);
    filetype = FslGetFileType(fslio);
    bufsize = x * y * z * v * (FslGetDataType(fslio,&dt) / 8);
    buffer = (void *) calloc(bufsize,1);
    FslReadVolumes(fslio,buffer,v);
    FslClose(fslio);
  }
  
  if (existingimage) {
    fslio = FslXOpen(filename,"wb",filetype);
  } else {
    fslio = FslOpen(filename,"wb");
    filetype = FslGetFileType(fslio);
  }

  if (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_MINC) {
    cerr << "Minc file type is not supported yet" << endl;
    exit(EXIT_FAILURE);
  }

  if (argc>3) {
    /* set uninteresting defaults */
    if (existingimage) {
      FslSetDataType(fslio,dt);
    } else {
      FslSetDataType(fslio,atoi(argv[12]));
    }
    FslSetDim(fslio,atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4])); 
    FslSetVoxDim(fslio,atof(argv[5]),atof(argv[6]),atof(argv[7]),atof(argv[8])); 
    
    {
      short short_array[100];
      short_array[0]=atoi(argv[9]);
      short_array[1]=atoi(argv[10]);
      short_array[2]=atoi(argv[11]);
      if ( (short_array[0]!=0) || (short_array[1]!=0) || (short_array[2]!=0) )
	{
           FslSetAnalyzeSform(fslio,short_array,
			 atof(argv[5]),atof(argv[6]),atof(argv[7]));
 	}
    }
    
  } else {
      /* read XML form */
      char *newstr, *oldfname, *oldiname;
      ifstream inputfile;
      
  
      if (strcmp(argv[1],"-")==0) {fileread=0;}
      newstr = (char *)calloc(10000,1);
      oldfname = (char *)calloc(10000,1);
      oldiname = (char *)calloc(10000,1);
      //hdrxml = (char *)calloc(65534,1);  /* too long, to be safe */
      hdrxml = new char[65534];
      if (fileread) 
      {
	inputfile.open (argv[1], ifstream::in | ifstream::binary);
        if (!inputfile.is_open())
        {      
	      cerr << "Cannot open file " << argv[1] << endl;
	      return EXIT_FAILURE;
	}
      }
       
      do 
      {
	  if (fileread) 
          {
	    inputfile.getline(newstr,9999);  // maybe use > for delimiting character remove while increase size
	  } 
          else 
          {
	      if (fgets(newstr,9999,stdin)==NULL) break;
	  }
	  strcat(hdrxml,newstr);
      }  while (strcmp(newstr + strlen(newstr) - 2,"/>")!=0);

      strcpy(oldfname,fslio->niftiptr->fname);
      strcpy(oldiname,fslio->niftiptr->iname);
      fslio->niftiptr = nifti_image_from_ascii(hdrxml);

      if (fslio->niftiptr == NULL) 
      {
	  cerr << "Incomplete or incorrect text: could not form header info" << endl;
	  return EXIT_FAILURE;
      }

      fslio->niftiptr->fname = oldfname;
      fslio->niftiptr->iname = oldiname;

      //free(hdrxml);
      delete [] hdrxml;
      if (fileread) inputfile.close();
  }

  /* reset filetype and datatype in case it has been overwritten */
  
  FslSetFileType(fslio,filetype);
  if (existingimage) {
    /* dt is only set if an image was previously read */ 
    FslSetDataType(fslio,dt);
  }

  fslio->niftiptr->byteorder = short_order();
//   if (strcmp(argv[argc-1],"-r")==0) { 
//       /* swap */
//       if (short_order()==MSB_FIRST) fslio->niftiptr->byteorder = LSB_FIRST;
//       else fslio->niftiptr->byteorder = MSB_FIRST;
//   }
  
  /* write header */
  
  FslWriteHeader(fslio);

  /* if previously read buffer is wrong size then make a zero image here */
  FslGetDim(fslio,&x,&y,&z,&v);
  if ( bufsize != ( x * y * z * v * (FslGetDataType(fslio,&dt)/8)) ) {
    if (bufsize>0) free(buffer);  /* only if previously read */
      buffer = (void *) calloc(x * y * z * v,FslGetDataType(fslio,&dt)/8);
  }

  /* write the data out - either from previous read or zeros */
  FslWriteVolumes(fslio,buffer,fslio->niftiptr->dim[4]);
  
  FslClose(fslio);
  return 0;
}


int main(int argc,char *argv[])
{
  if (argc != 14 && argc != 3) 
  {
    print_usage(string(argv[0]));
    return 1; 
  }
  return avwcreatehd_main(argc,argv); 
}


