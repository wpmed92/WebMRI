//     avwcpgeom++.cc - Copy certain parts of an AVW header
//     Mark Jenkinson, Steve Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2001-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"

using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: avwcpgeom <source> <destination> [-d]" << endl;
  cout << "-d : don't copy image dimensions" << endl;
}


int avwcpgeom_main(int argc, char *argv[])
{
  FSLIO *src = NULL, *dest = NULL, *destcopy = NULL;
  short x, y, z, v, copydim=1, t, scode, qcode, dt=-1;
  float vx, vy, vz, tr;
  int filetype;
  mat44 smat, qmat;
  void *buffer = NULL;
  char desthdrname[10000];
  size_t nbytes, nsrcbytes, ndestbytes;

  if (argc>3)
    copydim=0;
    
  src = FslOpen(argv[1], "r");

  dest = FslOpen(argv[2], "r");
  destcopy = (FSLIO *)calloc(sizeof(FSLIO),1);
  FslCloneHeader(destcopy,dest);

  if ((src == NULL) || (dest == NULL)) {
    perror("Error opening files");
    return EXIT_FAILURE;
  }

  FslGetDim(src, &x, &y, &z, &v);
  nsrcbytes = x * y * z * v * (FslGetDataType(dest, &t) / 8);
  FslGetDim(dest, &x, &y, &z, &v); 
  ndestbytes = x * y * z * v * (FslGetDataType(dest, &t) / 8);
  if (nsrcbytes > ndestbytes) nbytes=nsrcbytes; else nbytes=ndestbytes;
  if( (buffer = calloc(nbytes,1)) == NULL ) {
    perror("Unable to allocate memory for copy");
    return EXIT_FAILURE;
  }
  FslReadVolumes(dest, buffer, v);


  strcpy(desthdrname,dest->niftiptr->fname);
  filetype = FslGetFileType(dest);
  FslGetDataType(dest,&dt);
  FslClose(dest);
  dest = FslXOpen(desthdrname, "wb", filetype);
  FslCloneHeader(dest,destcopy); 

  scode = FslGetStdXform(src,&smat);
  FslSetStdXform(dest,scode,smat);
  qcode = FslGetRigidXform(src,&qmat);
  FslSetRigidXform(dest,qcode,qmat);
  
  FslGetVoxDim(src, &vx, &vy, &vz, &tr);
  FslSetVoxDim(dest, vx, vy, vz, tr);
  
  if (copydim) {
    FslGetDim(src, &x, &y, &z, &v);
    FslSetDim(dest, x, y, z, v);
  }


  /* Preserve the datatype - probably unneccesary now, but left for safety */ 
  FslSetDataType(dest,dt);

 FslWriteHeader(dest);
 FslWriteVolumes(dest, buffer, v);
 FslClose(dest);

 FslClose(src);

  return 0;
}


int main(int argc,char *argv[])
{
  if (argc < 3) 
  {
    print_usage(string(argv[0]));
    return 1; 
  }
  return avwcpgeom_main(argc,argv); 
}


