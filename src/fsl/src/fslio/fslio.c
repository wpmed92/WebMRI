/*
    fslio.c  (Input and output routines for images in FSL)

    Mark Jenkinson
    FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

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


#include "fslio.h"
#include "assert.h"

static int FslIgnoreMFQ=0;
static int FslOverrideOutputType=-1;

#define FSLIOERR(x) { fprintf(stderr,"Error:: %s\n",(x)); fflush(stderr); exit(EXIT_FAILURE); }


char* FslFileTypeString(int filetype)
{
  if (filetype==FSL_TYPE_ANALYZE)          return "ANALYZE-7.5";
  if (filetype==FSL_TYPE_NIFTI)            return "NIFTI-1+";
  if (filetype==FSL_TYPE_NIFTI_PAIR)       return "NIFTI-1";
  if (filetype==FSL_TYPE_ANALYZE_GZ)       return "ANALYZE-7.5";
  if (filetype==FSL_TYPE_NIFTI_GZ)         return "NIFTI-1+";
  if (filetype==FSL_TYPE_NIFTI_PAIR_GZ)    return "NIFTI-1";
  return "UNKNOWN";
}


int FslIsValidFileType(int filetype)
{
  if ( (filetype!=FSL_TYPE_ANALYZE)    && (filetype!=FSL_TYPE_ANALYZE_GZ) && 
       (filetype!=FSL_TYPE_NIFTI)      && (filetype!=FSL_TYPE_NIFTI_GZ) && 
       (filetype!=FSL_TYPE_NIFTI_PAIR) && (filetype!=FSL_TYPE_NIFTI_PAIR_GZ) &&
       (filetype!=FSL_TYPE_MINC)       && (filetype!=FSL_TYPE_MINC_GZ) ) {
    fprintf(stderr,"Error: unrecognised file type: %d\n",filetype);
    return 0;
  }
  return 1;
}


int FslBaseFileType(int filetype)
{
  /* returns -1 to indicate error - unrecognised filetype */
  if ( (filetype==FSL_TYPE_ANALYZE_GZ)    || (filetype==FSL_TYPE_ANALYZE) )
    return FSL_TYPE_ANALYZE;
  if ( (filetype==FSL_TYPE_NIFTI_GZ)      || (filetype==FSL_TYPE_NIFTI) )
    return FSL_TYPE_NIFTI;
  if ( (filetype==FSL_TYPE_NIFTI_PAIR_GZ) || (filetype==FSL_TYPE_NIFTI_PAIR) )
    return FSL_TYPE_NIFTI_PAIR;
  if ( (filetype==FSL_TYPE_MINC_GZ)       || (filetype==FSL_TYPE_MINC) )  
    return FSL_TYPE_MINC;
  fprintf(stderr,"Error: unrecognised file type (%d)\n",filetype);
  return -1;
}


int FslGetFileType2(const FSLIO *fslio, int quiet)
{
  FSLIO *mutablefslio;
  if (fslio==NULL)  FSLIOERR("FslGetFileType: Null pointer passed for FSLIO");
  if ( (fslio->file_mode==FSL_TYPE_MINC) || (fslio->file_mode==FSL_TYPE_MINC_GZ) ) {
    return fslio->file_mode;
  }
  if ( !FslIsValidFileType(fslio->file_mode) )    return -1;

  if (fslio->niftiptr!=NULL) {   /* check that it is nifti_type and filetype are consistent */
    if (fslio->niftiptr->nifti_type != FslBaseFileType(fslio->file_mode)) {
      if (!quiet) {
	fprintf(stderr,"Warning: nifti structure and fsl structure disagree on file type\n");
	fprintf(stderr,"nifti = %d and fslio = %d\n",fslio->niftiptr->nifti_type,fslio->file_mode);
      }
      mutablefslio = (FSLIO *) fslio;  /* dodgy and will generate warnings */
      mutablefslio->niftiptr->nifti_type = FslBaseFileType(fslio->file_mode);
      return fslio->file_mode;
    }
 }
  return fslio->file_mode;
}

int FslGetFileType(const FSLIO *fslio)
{ 
  return FslGetFileType2(fslio,0);
}



void FslSetFileType(FSLIO *fslio, int filetype)
{
  if (fslio==NULL)  FSLIOERR("FslSetFileType: Null pointer passed for FSLIO");
  if ( (filetype==FSL_TYPE_MINC) || (filetype==FSL_TYPE_MINC_GZ) ) {
    fslio->file_mode = filetype;
    return;
  }
  if (! FslIsValidFileType(filetype)) { return; } 
  fslio->file_mode = filetype;  /* indicates general nifti - details in niftiptr */
  if (fslio->niftiptr!=NULL) { 
    fslio->niftiptr->nifti_type = FslBaseFileType(filetype); 
    nifti_set_iname_offset(fslio->niftiptr);
  }
}



int FslIsSingleFileType(int filetype)
{
  if ( (filetype==FSL_TYPE_NIFTI) || (filetype==FSL_TYPE_NIFTI_GZ) || 
       (filetype==FSL_TYPE_MINC)  || (filetype==FSL_TYPE_MINC_GZ) )
    return 1;
  return 0;
}


int FslIsCompressedFileType(int filetype)
{
  if ( filetype >=100 ) return 1;
  return 0;
}


int FslGetWriteMode(const FSLIO *fslio)
{
  if (fslio==NULL)  FSLIOERR("FslGetWriteMode: Null pointer passed for FSLIO");
  return fslio->write_mode;
} 


void FslSetWriteMode(FSLIO *fslio, int mode)
{
  if (fslio==NULL)  FSLIOERR("FslSetWriteMode: Null pointer passed for FSLIO");
  fslio->write_mode = mode;
}


int FslGetEnvOutputType()
{
  /* return type is one of FSL_TYPE_* or -1 to indicate error */
  char *otype;
  if (FslOverrideOutputType>=0)  return FslOverrideOutputType;
  otype = getenv("FSLOUTPUTTYPE");
  if (otype == NULL) {
    fprintf(stderr,"ERROR:: Environment variable FSLOUTPUTTYPE is not set!\n");
    fprintf(stderr,"Please make sure that the appropriate configuration file is sourced by your shell (e.g. by putting it in .profile).\n");
    fprintf(stderr,"e.g. bash or sh users add the line \". ${FSLDIR}/etc/fslconf/fsl.sh\"\n");
    fprintf(stderr,"e.g. tcsh or csh users add the line \"source ${FSLDIR}/etc/fslconf/fsl.csh\"\n");
    exit(EXIT_FAILURE);
  }
  if (strcmp(otype,"ANALYZE")==0) { return FSL_TYPE_ANALYZE; }
  if (strcmp(otype,"ANALYZE_GZ")==0) { return FSL_TYPE_ANALYZE_GZ; }
  if (strcmp(otype,"NIFTI")==0) { return FSL_TYPE_NIFTI; }
  if (strcmp(otype,"NIFTI_GZ")==0) { return FSL_TYPE_NIFTI_GZ; }
  if (strcmp(otype,"NIFTI_PAIR")==0) { return FSL_TYPE_NIFTI_PAIR; }
  if (strcmp(otype,"NIFTI_PAIR_GZ")==0) { return FSL_TYPE_NIFTI_PAIR_GZ; }
  if (strcmp(otype,"MINC")==0) { return FSL_TYPE_MINC; }
  if (strcmp(otype,"MINC_GZ")==0) { return FSL_TYPE_MINC_GZ; }
  fprintf(stderr,"ERROR:: Unrecognised value (%s) of environment variable FSLOUTPUT\n",otype);
  fprintf(stderr,"Legal values are: ANALYZE, NIFTI, NIFTI_PAIR, MINC, ANALYZE_GZ, NIFTI_GZ, NIFTI_PAIR_GZ, MINC_GZ\n");
  exit(EXIT_FAILURE);
  return -1;
}
    

int FslFileType(const char* fname) 
{
  /* return type is FSL_TYPE_* or -1 to indicate undetermined */
  /* use name as first priority but if that is ambiguous then resolve using environment */
  int flen;
  int retval=-1;
  if (fname==NULL) return retval;
  flen = strlen(fname);
  if (flen<5) return retval;  /* smallest name + extension is a.nii */
  if (strcmp(fname + flen - 4,".nii")==0)  retval=FSL_TYPE_NIFTI;
  if (strcmp(fname + flen - 7,".nii.gz")==0)  retval=FSL_TYPE_NIFTI_GZ;
  if (strcmp(fname + flen - 4,".mnc")==0)  retval=FSL_TYPE_MINC;
  if (strcmp(fname + flen - 7,".mnc.gz")==0)  retval=FSL_TYPE_MINC;
  if (strcmp(fname + flen - 4,".hdr")==0)  retval=FSL_TYPE_NIFTI_PAIR;
  if (strcmp(fname + flen - 4,".img")==0)  retval=FSL_TYPE_NIFTI_PAIR;
  if (strcmp(fname + flen - 7,".hdr.gz")==0)  retval=FSL_TYPE_NIFTI_PAIR_GZ;
  if (strcmp(fname + flen - 7,".img.gz")==0)  retval=FSL_TYPE_NIFTI_PAIR_GZ;
  if ( (retval==FSL_TYPE_NIFTI_PAIR) || (retval==FSL_TYPE_NIFTI_PAIR_GZ) ) {
    /* If it was hdr or img, check if Analyze was requested by environment */
    if ( (FslGetEnvOutputType() == FSL_TYPE_ANALYZE) && (retval == FSL_TYPE_NIFTI_PAIR) ) 
      retval=FSL_TYPE_ANALYZE;
    if ( (FslGetEnvOutputType() == FSL_TYPE_ANALYZE_GZ) && (retval == FSL_TYPE_NIFTI_PAIR_GZ) )
      retval=FSL_TYPE_ANALYZE_GZ;
  }
  return retval;
}


int FslGetReadFileType(const FSLIO *fslio)
{
  /* This function is used to return the best estimate of the true file type once
   a simple open has occurred - for now it is used after a nifti open call is made */
  int filetype=FSL_TYPE_ANALYZE;  /* unused default */
  if (fslio==NULL)  FSLIOERR("FslReadGetFileType: Null pointer passed for FSLIO");
  /* Don't use fslio->file_mode as it hasn't been set yet */
  if (fslio->niftiptr!=NULL) {   
    /* use the nifti_type and hdr or img name to determine the actual type */
    if (fslio->niftiptr->nifti_type == FSL_TYPE_ANALYZE) {
      if (FslIsCompressedFileType(FslFileType(fslio->niftiptr->iname))) {
	filetype = FSL_TYPE_ANALYZE_GZ;
      } else {
	filetype = FSL_TYPE_ANALYZE;
      }
    }
    if (fslio->niftiptr->nifti_type == FSL_TYPE_NIFTI_PAIR) {
      if (FslIsCompressedFileType(FslFileType(fslio->niftiptr->iname))) {
	filetype = FSL_TYPE_NIFTI_PAIR_GZ;
      } else {
	filetype = FSL_TYPE_NIFTI_PAIR;
      }
    }
    if (fslio->niftiptr->nifti_type == FSL_TYPE_NIFTI) {
      if (FslIsCompressedFileType(FslFileType(fslio->niftiptr->fname))) {
	filetype = FSL_TYPE_NIFTI_GZ;
      } else {
	filetype = FSL_TYPE_NIFTI;
      }
    }
    
  }
  if (fslio->mincptr!=NULL) { 
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    filetype = FSL_TYPE_MINC;
  }
  return filetype;
}


int FslFileExists(const char *filename)
{ 
  /* return 1 if file(s) exists, otherwise return 0 */
  char *hdrname = nifti_findhdrname(filename);
  char *imgname = NULL;
  if ((hdrname!=NULL) && 
      ((imgname = nifti_findimgname(filename, FslBaseFileType(FslFileType(hdrname))))!=NULL))
    { 
      free(hdrname); 
      free(imgname); 
      return 1; 
    }
  return 0;
}

char *FslMakeBaseName(const char *fname)
{
  char *basename;
  int blen;
  basename = nifti_makebasename(fname);
  blen = strlen(basename);
#ifdef HAVE_ZLIB
  if (strcmp(basename + blen-7,".mnc.gz") == 0) 
    { basename[blen-7]='\0'; return basename; }
#endif
  if (strcmp(basename + blen-4,".mnc") == 0) 
    { basename[blen-4]='\0'; return basename; }
  return basename;
}


void FslGetHdrImgNames(const char* filename, const FSLIO* fslio, 
		       char** hdrname, char** imgname)
{
  char *basename;
  int filetype;
  basename = FslMakeBaseName(filename);
  *hdrname = (char *)calloc(sizeof(char),strlen(basename)+8);
  *imgname = (char *)calloc(sizeof(char),strlen(basename)+8);
  strcpy(*hdrname,basename);
  strcpy(*imgname,basename);
  filetype = FslGetFileType(fslio);
  if (filetype==FSL_TYPE_NIFTI_GZ) {
    strcat(*hdrname,".nii.gz");
    strcat(*imgname,".nii.gz");
    free(basename);
    return;
  }
  if (filetype==FSL_TYPE_NIFTI) {
    strcat(*hdrname,".nii");
    strcat(*imgname,".nii");
    free(basename);
    return;
  }
  if (filetype==FSL_TYPE_MINC_GZ) {
    strcat(*hdrname,".mnc.gz");
    strcat(*imgname,".mnc.gz");
    free(basename);
    return;
  }
  if (filetype==FSL_TYPE_MINC) {
    strcat(*hdrname,".mnc");
    strcat(*imgname,".mnc");
    free(basename);
    return;
  }
  if ( (filetype==FSL_TYPE_NIFTI_PAIR_GZ) || (filetype==FSL_TYPE_ANALYZE_GZ) ) {
    strcat(*hdrname,".hdr.gz");
    strcat(*imgname,".img.gz");
    free(basename);
    return;
  }
  if ( (filetype==FSL_TYPE_NIFTI_PAIR) || (filetype==FSL_TYPE_ANALYZE) ) {
    strcat(*hdrname,".hdr");
    strcat(*imgname,".img");
    free(basename);
    return;
  }

  fprintf(stderr,"Error: Unrecognised filetype (%d)\n",FslGetFileType(fslio));
  free(basename);
  /* Failure */
  *hdrname = NULL;
  *imgname = NULL;
}


void FslSetInit(FSLIO* fslio)
{
  /* set some sensible defaults */
  fslio->niftiptr = NULL;
  fslio->mincptr  = NULL;
  FslSetFileType(fslio,FslGetEnvOutputType());
  FslSetWriteMode(fslio,0);
  fslio->written_hdr = 0;
}


FSLIO *FslInit()
{
  /* set up space for struct and some sensible defaults */
  FSLIO *fslio;
  fslio = (FSLIO *) calloc(1,sizeof(FSLIO));
  FslSetInit(fslio);
  return fslio;
}



void FslInit4Write(FSLIO* fslio, const char* filename, int ft)
{
  /* ft determines filetype if ft>=0*/ 
  int imgtype;

  FslSetWriteMode(fslio,1);

  /* Determine file type from image name (first priority) or environment (default) */
  imgtype = FslFileType(filename);
  if (imgtype<0)  imgtype = FslGetEnvOutputType();

  if (ft >= 0) imgtype = ft;

  if (!FslIsValidFileType(imgtype)) {
    fprintf(stderr,"Error: Failed to determine file type for writing in FslOpen()\n");
    exit(EXIT_FAILURE);
  }
  
  if ( (FslBaseFileType(imgtype)!=FSL_TYPE_MINC) ) {
    FslInitHeader(fslio, NIFTI_TYPE_FLOAT32,
		  1, 1, 1, 3,  0.0, 0.0, 0.0, 0.0,  4, "mm");
    
    FslSetFileType(fslio,imgtype);  /* this is after InitHeader as niftiptr set there */
    
    /* determine the header and image filename */
    FslGetHdrImgNames(filename,fslio,&(fslio->niftiptr->fname),&(fslio->niftiptr->iname));
    if ( (fslio->niftiptr->fname == NULL) || (fslio->niftiptr->iname == NULL) ) { 
      fprintf(stderr,"Error: cannot find filenames for %s\n",filename); 
    }

  } else if (FslBaseFileType(imgtype)==FSL_TYPE_MINC) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    return;
  } else {
    fprintf(stderr,"Error:: unrecognised image type requested\n");
    return;
  }

  return;
}



void FslInitHeader(FSLIO *fslio, short t, 
		   size_t x, size_t y, size_t z, size_t v,
		   float vx, float vy, float vz, float tr,
		   size_t dim,
		   const char* units)
{
  /* NB: This function does not set the file type or write mode*/

  if (fslio==NULL)  FSLIOERR("FslInitHeader: Null pointer passed for FSLIO");
  
  fslio->niftiptr = nifti_simple_init_nim();
  /* make nifti type consistent with fslio */
  fslio->niftiptr->nifti_type = FslBaseFileType(fslio->file_mode);

  fslio->mincptr = NULL;

  FslSetDataType(fslio,t);
  FslSetDim(fslio,x,y,z,v);
  FslSetVoxDim(fslio,vx,vy,vz,tr);
  FslSetTimeUnits(fslio,"s");
  FslSetDimensionality(fslio,dim);
}


void FslCloneHeader(FSLIO *dest, const FSLIO *src)
{
  /* only clone the information that is stored in the disk version of the header */
  /*  - therefore _not_ the filenames, output type, write mode, etc */ 

  char *fname=NULL, *iname=NULL;
  void *data=NULL;
  int filetype, writemode;
  int preserve_nifti_values = 0;
  if (dest==NULL)  FSLIOERR("FslCloneHeader: Null pointer passed for FSLIO");
  if (src==NULL)   FSLIOERR("FslCloneHeader: Null pointer passed for FSLIO");

  if (src->niftiptr!=NULL) {
    /* preserve the filenames, output type and write mode */
    if (dest->niftiptr != NULL) {
      fname = dest->niftiptr->fname;
      iname = dest->niftiptr->iname;
      data = dest->niftiptr->data;
      preserve_nifti_values = 1;
    }
    filetype = FslGetFileType2(dest,1);
    writemode = FslGetWriteMode(dest);

    /* copy _all_ info across */
    dest->niftiptr = nifti_copy_nim_info(src->niftiptr);

    /* restore old values */
    if (preserve_nifti_values) {
      dest->niftiptr->fname = fname;
      dest->niftiptr->iname = iname; 
      dest->niftiptr->data = data; 
    } else { 
	/* destroy the values that the nifti copy creates */
      dest->niftiptr->fname = NULL;
      dest->niftiptr->iname = NULL; 
      dest->niftiptr->data = NULL; 
    }
    FslSetFileType(dest,filetype);
    FslSetWriteMode(dest,writemode);
  }

  if (src->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }

}


int  fsl_fileexists(const char* fname)
{
  znzFile fp;
  int retval = 0;
  fp = znzopen( fname , "rb" , 1 );
  retval = !znz_isnull(fp);
  znzclose(fp);

  return retval;
}


int FslCheckForMultipleFileNames(const char* filename)
{
  char *basename, *tmpname;
  int singlecount=0, hdrcount=0, imgcount=0, ambiguous=0;
  basename =  nifti_makebasename(filename);
  tmpname = (char *)calloc(strlen(basename) + 10,sizeof(char));

  strcpy(tmpname,basename);
  strcat(tmpname,".nii"); 
  if (fsl_fileexists(tmpname)) { singlecount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".nii.gz"); 
  if (fsl_fileexists(tmpname)) { singlecount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".mnc"); 
  if (fsl_fileexists(tmpname)) { singlecount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".mnc.gz"); 
  if (fsl_fileexists(tmpname)) { singlecount++; }

  strcpy(tmpname,basename);
  strcat(tmpname,".img"); 
  if (fsl_fileexists(tmpname)) { imgcount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".img.gz"); 
  if (fsl_fileexists(tmpname)) { imgcount++; }

  strcpy(tmpname,basename);
  strcat(tmpname,".hdr"); 
  if (fsl_fileexists(tmpname)) { hdrcount++; }
  strcpy(tmpname,basename);
  strcat(tmpname,".hdr.gz"); 
  if (fsl_fileexists(tmpname)) { hdrcount++; }
 
  ambiguous = 1;
  if ( (hdrcount==1) && (imgcount==1) && (singlecount==0) )  { ambiguous=0; }
  if ( (hdrcount==0) && (imgcount==0) && (singlecount==1) )  { ambiguous=0; }

  /* treat no image found as not ambiguous - want opening errors instead */
  if ( (hdrcount==0) && (imgcount==0) && (singlecount==0) )  { ambiguous=0; }

  free(tmpname);
  free(basename);
  return ambiguous;
}



int check_for_multiple_filenames(const char* filename)
{
  char *basename, *tmpname;
  char *otype;
  if (FslCheckForMultipleFileNames(filename))
    {  /* take action */
      basename =  nifti_makebasename(filename);
      tmpname = (char *)calloc(strlen(basename) + 10,sizeof(char));
      fprintf(stderr,"\n\n\nWARNING!!!! Multiple image files detected:\n");
      /* list the offending files */
      strcpy(tmpname,basename);
      strcat(tmpname,".nii"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".nii.gz"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".mnc"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".mnc.gz"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".img"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".img.gz"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".hdr"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      strcpy(tmpname,basename);
      strcat(tmpname,".hdr.gz"); 
      if (fsl_fileexists(tmpname)) { fprintf(stderr,"%s ",tmpname); }
      fprintf(stderr,"\n\n");

      if (!FslIgnoreMFQ) {
	otype = getenv("FSLMULTIFILEQUIT");
	if (otype!=NULL) {
	  fprintf(stderr,"STOPPING PROGRAM\n");
	  exit(EXIT_FAILURE);
	}
      }
      return 1;
    }
  return 0;
}



FSLIO *FslOpen(const char *filename, const char *opts)
{
  /* Note: -1 for filetype indicates that FslXOpen should determine filetype for itself */
  return FslXOpen(filename,opts,-1);
}


FSLIO *FslXOpen(const char *filename, const char *opts, int filetype)
{
  /* Opens a file for either reading or writing */
  /* The filetype specifies the type of file to be written */
  /*    - legals values are as defined by FSL_TYPE_* */
  /*    - if value is less than zero, then it is ignored and the type is determined */
  /*          by the filename extension or, failing that, the environment default */
  /*    - files to be read are automatically read whether compressed or not */
  /* Also, reading uses the file extension and will fail if that file does not exist */
  /*    - for a more robust read, pass the basename in as then all types will be tried */

  FSLIO *fslio;
  char bopts[1024];
  size_t i, bi;
  int imgtype;

  fslio = FslInit();

  bi=0;
  for(i=0;i<strlen(opts);i++) {
    if (opts[i]=='w') { FslSetWriteMode(fslio,1); }
    if (opts[i]!='b' && opts[i]!='t') { bopts[bi++]=opts[i]; }
  }
  /* add in 'b' (at the end) for windows compatibility */
  bopts[bi++]='b';
  bopts[bi]='\0';
  

  if (FslGetWriteMode(fslio)==1) {
    
    /** ====================== Open file for writing ====================== **/
   
    FslInit4Write(fslio,filename,filetype);
    imgtype = FslGetFileType(fslio);
    fslio->written_hdr = 0;

    /* open the image file - not the header */
    fslio->fileptr = znzopen(fslio->niftiptr->iname,bopts,FslIsCompressedFileType(imgtype));
    if (znz_isnull(fslio->fileptr)) { 
      fprintf(stderr,"Error: failed to open file %s\n",fslio->niftiptr->iname); 
      return NULL;
    }

    if (!FslIsSingleFileType(imgtype)) {
      /* set up pointer at end of iname_offset for dual file formats (not singles) */
      FslSeekVolume(fslio,0);
    }
    return fslio;

  }



  /** ======================== Open file for reading ====================== **/

  check_for_multiple_filenames(filename);

  /* see if the extension indicates a minc file */
  imgtype = FslFileType(filename);
  if ((imgtype>=0) && (FslBaseFileType(imgtype)==FSL_TYPE_MINC)) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    return NULL;
  }

  /* otherwise open nifti file: read header and open img file (may be same file) */
  fslio->fileptr = nifti_image_open(filename,bopts,&(fslio->niftiptr));
  if (znz_isnull(fslio->fileptr)) { 
    fprintf(stderr,"Error: failed to open file %s\n",filename); 
    return NULL;
  }

  /* set the file type given what has been read - it uses nifti_type and filenames */
  imgtype = FslGetReadFileType(fslio);
  FslSetFileType(fslio,imgtype);
  FslSetWriteMode(fslio,0);

  /* if it is a nifti file but has inconsistent left-right ordering in the sform and qform then complain/crash */
  if (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_NIFTI) {
    if (FslGetLeftRightOrder(fslio) == FSL_INCONSISTENT)
      {
	fprintf(stderr,"ERROR: inconsistent left-right order stored in sform and qform in file %s\n",filename); 
	fprintf(stderr,"       Using sform instead of qform values\n\n");
	/* return NULL; */
      }
  }


  if (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE) {
    /* For the ANALYZE case in FSL, must cheat and grab the originator field! */
    /* Note that the header file is always separate here and closed by now */
    struct dsr ahdr;
    short orig[5];
    FslReadRawHeader(&ahdr,fslio->niftiptr->fname);
    if (fslio->niftiptr->byteorder != short_order()) {
      AvwSwapHeader(&ahdr);
    }
    /* Read the origin and set the sform up (if origin is non-zero) */
    /* Note that signed pixdims are passed in to set the LR orientation */
    memcpy(orig,&(ahdr.hist.originator),10);
    FslSetAnalyzeSform(fslio, orig, fslio->niftiptr->pixdim[1],
		       fslio->niftiptr->pixdim[2], fslio->niftiptr->pixdim[3]);
  }

  /* from now on force all vox dims to be positive - LR info is in sform */
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->dx = fabs(fslio->niftiptr->dx);
    fslio->niftiptr->dy = fabs(fslio->niftiptr->dy);
    fslio->niftiptr->dz = fabs(fslio->niftiptr->dz);
    fslio->niftiptr->pixdim[1] = fabs(fslio->niftiptr->pixdim[1]);
    fslio->niftiptr->pixdim[2] = fabs(fslio->niftiptr->pixdim[2]);
    fslio->niftiptr->pixdim[3] = fabs(fslio->niftiptr->pixdim[3]);
  }
  /* set up pointer at end of iname_offset , ready for reading */
  FslSeekVolume(fslio,0);  

  return fslio;

}



void* FslReadAllVolumes(FSLIO* fslio, char* filename)
{
  /* this reads both header and image - no need for Open or Close calls */ 
  /*  - it returns the location of data block (allocated by this function) */
  /*  - the best call to make before this is FslInit() or a calloc() for fslio */

  int imgtype;
  if (fslio==NULL)  FSLIOERR("FslReadAllVolumes: Null pointer passed for FSLIO");

  /* see if the extension indicates a minc file */
  imgtype = FslFileType(filename);
  if ((imgtype>=0) && (FslBaseFileType(imgtype)==FSL_TYPE_MINC)) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    return NULL;
  }

  /** otherwise it is a nifti file - so read it! **/
  fslio->mincptr = NULL;
  /* make sure an FslOpen hasn't locked the file */
  if (!znz_isnull(fslio->fileptr)) FslClose(fslio);  
  
  fslio->niftiptr = nifti_image_read(filename,1);
  FslSetFileType(fslio,fslio->niftiptr->nifti_type);
  FslSetWriteMode(fslio,0);
  return fslio->niftiptr->data;
}



size_t FslReadVolumes(FSLIO *fslio, void *buffer, size_t nvols)
{
  /* Returns number of volumes read */
  int volbytes;
  size_t retval=0;
  if (fslio==NULL)  FSLIOERR("FslReadVolumes: Null pointer passed for FSLIO");
  if (znz_isnull(fslio->fileptr))  FSLIOERR("FslReadVolumes: Null file pointer");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->data = buffer;
    volbytes = FslGetVolSize(fslio)  * fslio->niftiptr->nbyper;
    retval = nifti_read_buffer(fslio->fileptr,fslio->niftiptr->data,nvols*volbytes,fslio->niftiptr);
    retval /= volbytes;
  }

  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return retval;
}



void FslWriteAllVolumes(FSLIO *fslio, const void *buffer)
{
  short x,y,z,t;

  /* writes all data from buffer (using size info from fslio) to file */
  if (fslio==NULL)  FSLIOERR("FslWriteAllVolumes: Null pointer passed for FSLIO");

  FslGetDim(fslio,&x,&y,&z,&t);
  FslWriteHeader(fslio);
  FslWriteVolumes(fslio,buffer,t);
  return;
}



size_t FslWriteVolumes(FSLIO *fslio, const void *buffer, size_t nvols)
{
  /* The dimensions and datatype must be set before calling this function */
  int retval=0;
  if (fslio==NULL)  FSLIOERR("FslWriteVolumes: Null pointer passed for FSLIO");
  if ( (!fslio->written_hdr) && (FslIsSingleFileType(FslGetFileType(fslio))) &&
       (FslIsCompressedFileType(FslGetFileType(fslio))) )
    { FSLIOERR("FslWriteVolumes: header must be written before data for single compressed file types"); }
  
  if (fslio->niftiptr!=NULL) {
    long int nbytes, bpv;
    bpv = fslio->niftiptr->nbyper;  /* bytes per voxel */
    nbytes = nvols * FslGetVolSize(fslio) * bpv;

    if ( (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE)
	 && (FslGetLeftRightOrder(fslio)==FSL_NEUROLOGICAL) ) {
      /* If it is Analyze and Neurological order then SWAP DATA into Radiological order */
      /* This is nasty - but what else can be done?!? */
      char *tmpbuf, *inbuf;
      long int x, b, n, nrows;
      short nx, ny, nz, nv;
      inbuf = (char *) buffer;
      tmpbuf = (char *)calloc(nbytes,1);
      FslGetDim(fslio,&nx,&ny,&nz,&nv);
      nrows = nbytes / (nx * bpv);
      for (n=0; n<nrows; n++) {
	for (x=0; x<nx; x++) {
	  for (b=0; b<bpv; b++) {
	    tmpbuf[b +  bpv * (n*nx + nx - 1 - x)] = inbuf[b + bpv * (n*nx + x)];
	  }
	}
      }
      retval = nifti_write_buffer(fslio->fileptr, tmpbuf, nbytes);
      free(tmpbuf);
    } else {
      retval = nifti_write_buffer(fslio->fileptr, buffer, nbytes);
    }
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return retval;
}


void FslWriteHeader(FSLIO *fslio)
{
  short sform_code, qform_code;
  mat44 smat, qmat;

  /* writes header and opens img file ready for writing */
  if (fslio==NULL)  FSLIOERR("FslWriteHeader: Null pointer passed for FSLIO");

  if (fslio->niftiptr!=NULL) {
    fslio->written_hdr = 1;
    if (znz_isnull(fslio->fileptr)) FSLIOERR("FslWriteHeader: no file opened!");
    
    /* modify niftiptr for FSL-specific purposes */
    strcpy(fslio->niftiptr->descrip,"FSL3.3");
    /* set qform to equal sform if currently unset (or vice versa) */
    qform_code = FslGetRigidXform(fslio,&qmat);
    sform_code = FslGetStdXform(fslio,&smat);
    if ( (sform_code != NIFTI_XFORM_UNKNOWN) && 
	 (qform_code == NIFTI_XFORM_UNKNOWN) ) {
      FslSetRigidXform(fslio,sform_code,smat);
    }
    if ( (qform_code != NIFTI_XFORM_UNKNOWN) && 
	 (sform_code == NIFTI_XFORM_UNKNOWN) ) {
      FslSetStdXform(fslio,qform_code,qmat);
    }

    if (FslIsSingleFileType(FslGetFileType(fslio))) {
      /* write header info but don't close the file */
      nifti_image_write_hdr_img2(fslio->niftiptr,2,"wb",&(fslio->fileptr));
      /* set up pointer at end of iname_offset for single files only */
      FslSeekVolume(fslio,0);
    } else {
      /* open a new hdr file, write it and close it */
      nifti_image_write_hdr_img(fslio->niftiptr,0,"wb");
    }
  }

  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return;
}


size_t FslReadSliceSeries(FSLIO *fslio, void *buffer, short slice, size_t nvols)
{
  size_t slbytes,volbytes;
  size_t n, orig_offset;
  short x,y,z,v,type;

  if (fslio==NULL)  FSLIOERR("FslReadSliceSeries: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    
    FslGetDim(fslio,&x,&y,&z,&v);
    
    if ((slice<0) || (slice>=z)) FSLIOERR("FslReadSliceSeries: slice outside valid range");
    
    slbytes = x * y * (FslGetDataType(fslio, &type) / 8);
    volbytes = slbytes * z;
    
    orig_offset = znztell(fslio->fileptr);
    znzseek(fslio->fileptr, slbytes*slice, SEEK_CUR);
    
    for (n=0; n<nvols; n++) {
      if (n>0) znzseek(fslio->fileptr, volbytes - slbytes, SEEK_CUR);
      if (znzread((char *)buffer+n*slbytes, 1, slbytes, fslio->fileptr) != slbytes)
	FSLIOERR("FslReadSliceSeries: failed to read values");
     if (fslio->niftiptr->byteorder != short_order())
	swap_Nbytes(slbytes / fslio->niftiptr->swapsize, fslio->niftiptr->swapsize,
		    (char *)buffer+n*slbytes);
     }
    
    
    /* restore file pointer to original position */
    znzseek(fslio->fileptr,orig_offset,SEEK_SET);
    return n;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


size_t FslReadRowSeries(FSLIO *fslio, void *buffer, short row, short slice, size_t nvols)
{
  size_t rowbytes,slbytes,volbytes;
  size_t n, orig_offset;
  short x,y,z,v,type;
  
  if (fslio==NULL)  FSLIOERR("FslReadRowSeries: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    
    FslGetDim(fslio,&x,&y,&z,&v);
    
    if ((slice<0) || (slice>=z)) FSLIOERR("FslReadRowSeries: slice outside valid range");
    if ((row<0) || (row>=y)) FSLIOERR("FslReadRowSeries: row outside valid range");
    
    rowbytes = x * (FslGetDataType(fslio, &type)) / 8;
    slbytes = rowbytes * y;
    volbytes = slbytes * z;
    
    orig_offset = znztell(fslio->fileptr);
    znzseek(fslio->fileptr, rowbytes*row + slbytes*slice, SEEK_CUR);
    
    for (n=0; n<nvols; n++){
      if (n>0) znzseek(fslio->fileptr, volbytes - rowbytes, SEEK_CUR);
      if (znzread((char *)buffer+n*rowbytes, 1, rowbytes, fslio->fileptr) != rowbytes)
	FSLIOERR("FslReadRowSeries: failed to read values");
      if (fslio->niftiptr->byteorder != short_order())
	swap_Nbytes(rowbytes / fslio->niftiptr->swapsize, fslio->niftiptr->swapsize,
		    (char *)buffer+n*rowbytes);
    }
    
    /* restore file pointer to original position */
    znzseek(fslio->fileptr,orig_offset,SEEK_SET);
    return n;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


size_t FslReadTimeSeries(FSLIO *fslio, void *buffer, short xVox, short yVox, short zVox, 
			 size_t nvols)
{
  size_t volbytes, offset, orig_offset;
  size_t n;
  short xdim,ydim,zdim,v,wordsize;

  if (fslio==NULL)  FSLIOERR("FslReadTimeSeries: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {

    FslGetDim(fslio,&xdim,&ydim,&zdim,&v);
    
    if ((xVox<0) || (xVox >=xdim)) FSLIOERR("FslReadTimeSeries: voxel outside valid range");
    if ((yVox<0) || (yVox >=ydim)) FSLIOERR("FslReadTimeSeries: voxel outside valid range");
    if ((zVox<0) || (zVox >=zdim)) FSLIOERR("FslReadTimeSeries: voxel outside valid range");
    
    wordsize = fslio->niftiptr->nbyper;
    volbytes = xdim * ydim * zdim * wordsize;
    
    orig_offset = znztell(fslio->fileptr);
    offset = ((ydim * zVox + yVox) * xdim + xVox) * wordsize;
    znzseek(fslio->fileptr,offset,SEEK_CUR);
    
    for (n=0; n<nvols; n++) {
      if (n>0) znzseek(fslio->fileptr, volbytes - wordsize, SEEK_CUR);
      if (znzread((char *)buffer+(n*wordsize), 1, wordsize,fslio->fileptr) != wordsize)
	FSLIOERR("FslReadTimeSeries: failed to read values"); 
      if (fslio->niftiptr->byteorder != short_order())
	swap_Nbytes(1,fslio->niftiptr->swapsize,(char *)buffer+(n*wordsize));
    }
    
    /* restore file pointer to original position */
    znzseek(fslio->fileptr,orig_offset,SEEK_SET);
    return n;

  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


size_t FslReadCplxVolumes(FSLIO *fslio, void *buffer, size_t nvols, char mode)
{
  if (fslio==NULL)  FSLIOERR("FslReadCplxVolumes: Null pointer passed for FSLIO");
  fprintf(stderr,"Warning:: FslReadCplxVolumes is not yet supported\n");
  return 0;
}

size_t FslWriteCplxVolumes(FSLIO *fslio, void *buffer, size_t nvols, char mode)
{
  if (fslio==NULL)  FSLIOERR("FslWriteCplxVolumes: Null pointer passed for FSLIO");
  fprintf(stderr,"Warning:: FslWriteCplxVolumes is not yet supported\n");
  return 0;
}

int FslSeekVolume(FSLIO *fslio, size_t vols)
{
  int offset;
  if (fslio==NULL)  FSLIOERR("FslSeekVolume: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    offset = fslio->niftiptr->iname_offset + 
      vols * FslGetVolSize(fslio) * fslio->niftiptr->nbyper;
    if (znz_isnull(fslio->fileptr)) FSLIOERR("FslSeekVolume: Null file pointer");
    return znzseek(fslio->fileptr,offset,SEEK_SET);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


size_t FslGetVolSize(FSLIO *fslio)
{
  /* returns number of voxels per 3D volume */
  if (fslio==NULL)  FSLIOERR("FslGetVolSize: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    return (fslio->niftiptr->nx * fslio->niftiptr->ny * fslio->niftiptr->nz);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
}


void FslSetDim(FSLIO *fslio, short x, short y, short z, short v)
{
  int ndim;
  if (fslio==NULL)  FSLIOERR("FslSetDim: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {

    ndim=4;
    if (v<=1) {ndim--; if (z<=1) {ndim--; if (y<=1) {ndim--; if (x<=1) {ndim--;}}}}

    fslio->niftiptr->ndim = ndim;

    if (x>=1) fslio->niftiptr->nx = x; else fslio->niftiptr->nx=1;
    if (y>=1) fslio->niftiptr->ny = y; else fslio->niftiptr->ny=1;
    if (z>=1) fslio->niftiptr->nz = z; else fslio->niftiptr->nz=1;
    if (v>=1) fslio->niftiptr->nt = v; else fslio->niftiptr->nt=1;
    fslio->niftiptr->nu = 1;
    fslio->niftiptr->nv = 1;
    fslio->niftiptr->nw = 1;

    /* deal with stupid redundancies */
    fslio->niftiptr->dim[0] = fslio->niftiptr->ndim ;
    fslio->niftiptr->dim[1] = fslio->niftiptr->nx;
    fslio->niftiptr->dim[2] = fslio->niftiptr->ny;
    fslio->niftiptr->dim[3] = fslio->niftiptr->nz;
    fslio->niftiptr->dim[4] = fslio->niftiptr->nt;
    fslio->niftiptr->dim[5] = fslio->niftiptr->nu;
    fslio->niftiptr->dim[6] = fslio->niftiptr->nv;
    fslio->niftiptr->dim[7] = fslio->niftiptr->nw;

    fslio->niftiptr->nvox =  fslio->niftiptr->nx * fslio->niftiptr->ny * fslio->niftiptr->nz
      * fslio->niftiptr->nt * fslio->niftiptr->nu * fslio->niftiptr->nv * fslio->niftiptr->nw ;

  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetDim(FSLIO *fslio, short *x, short *y, short *z, short *v)
{
  if (fslio==NULL)  FSLIOERR("FslGetDim: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *x = fslio->niftiptr->nx;
    *y = fslio->niftiptr->ny;
    *z = fslio->niftiptr->nz;
    *v = fslio->niftiptr->nt;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetDimensionality(FSLIO *fslio, size_t dim)
{
  if (fslio==NULL)  FSLIOERR("FslSetDimensionality: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->ndim = dim;
    fslio->niftiptr->dim[0] = dim;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetDimensionality(FSLIO *fslio, size_t *dim)
{
  if (fslio==NULL)  FSLIOERR("FslGetDimensionality: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *dim = fslio->niftiptr->ndim;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetVoxDim(FSLIO *fslio, float x, float y, float z, float tr)
{
  if (fslio==NULL)  FSLIOERR("FslSetVoxDim: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->dx = fabs(x);
    fslio->niftiptr->dy = fabs(y);
    fslio->niftiptr->dz = fabs(z);
    fslio->niftiptr->dt = fabs(tr);
    fslio->niftiptr->pixdim[1] = fabs(x);
    fslio->niftiptr->pixdim[2] = fabs(y);
    fslio->niftiptr->pixdim[3] = fabs(z);
    fslio->niftiptr->pixdim[4] = fabs(tr);
    /* set the units to mm and seconds */
    fslio->niftiptr->xyz_units  = NIFTI_UNITS_MM;
    fslio->niftiptr->time_units = NIFTI_UNITS_SEC;
 }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetVoxDim(FSLIO *fslio, float *x, float *y, float *z, float *tr)
{
  if (fslio==NULL)  FSLIOERR("FslGetVoxDim: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *x = fabs(fslio->niftiptr->dx);
    *y = fabs(fslio->niftiptr->dy);
    *z = fabs(fslio->niftiptr->dz);
    *tr = fabs(fslio->niftiptr->dt);
    /* now check the units and convert to mm and sec */
    if (fslio->niftiptr->xyz_units == NIFTI_UNITS_METER) 
    { *x *= 1000.0;   *y *= 1000.0;   *z *= 1000.0; }
    if (fslio->niftiptr->xyz_units == NIFTI_UNITS_MICRON) 
    { *x /= 1000.0;   *y /= 1000.0;   *z /= 1000.0; }
    if (fslio->niftiptr->xyz_units == NIFTI_UNITS_MSEC) 
    { *tr /= 1000.0; }
    if (fslio->niftiptr->xyz_units == NIFTI_UNITS_USEC) 
    { *tr /= 1000000.0; }
    /* if it is Hz or other frequency then leave it */
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetCalMinMax(FSLIO *fslio, float *min, float *max)
{
  if (fslio==NULL)  FSLIOERR("FslGetCalMinMax: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *min = fslio->niftiptr->cal_min;
    *max = fslio->niftiptr->cal_max;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetCalMinMax(FSLIO *fslio, float  min, float  max)
{
  if (fslio==NULL)  FSLIOERR("FslSetCalMinMax: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->cal_min = min;
    fslio->niftiptr->cal_max = max;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetAuxFile(FSLIO *fslio,char *aux_file)
{
  if (fslio==NULL)  FSLIOERR("FslGetAuxFile: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    strncpy(aux_file,fslio->niftiptr->aux_file, 24);
    aux_file[23] = '\0';
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetAuxFile(FSLIO *fslio,const char *aux_file)
{
  if (fslio==NULL)  FSLIOERR("FslSetAuxFile: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    strncpy(fslio->niftiptr->aux_file, aux_file, 24);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetVoxUnits(FSLIO *fslio, const char *units)
{
  int unitcode=0;
  if (fslio==NULL)  FSLIOERR("FslSetVoxUnits: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    if (strcmp(units,nifti_units_string(NIFTI_UNITS_METER))==0) {
      unitcode = NIFTI_UNITS_METER;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_MM))==0) {
      unitcode = NIFTI_UNITS_MM;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_MICRON))==0) {
      unitcode = NIFTI_UNITS_MICRON;
    }
    fslio->niftiptr->xyz_units  = unitcode;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetVoxUnits(FSLIO *fslio, char *units)
{
  if (fslio==NULL)  FSLIOERR("FslGetVoxUnits: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    strcpy(units,nifti_units_string(fslio->niftiptr->xyz_units));
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}

void FslSetTimeUnits(FSLIO *fslio, const char *units)
{
  int unitcode=0;
  if (fslio==NULL)  FSLIOERR("FslSetTimeUnits: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    if (strcmp(units,nifti_units_string(NIFTI_UNITS_HZ))==0) {
      unitcode = NIFTI_UNITS_HZ;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_PPM))==0) {
      unitcode = NIFTI_UNITS_PPM;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_RADS))==0) {
      unitcode = NIFTI_UNITS_RADS;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_SEC))==0) {
      unitcode = NIFTI_UNITS_SEC;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_MSEC))==0) {
	fprintf(stderr,"Warning::Setting time units to msec is not fully recommended in fslio\n");
	unitcode = NIFTI_UNITS_MSEC;
    } else if (strcmp(units,nifti_units_string(NIFTI_UNITS_USEC))==0) {
	fprintf(stderr,"Warning::Setting time units to msec is not fully recommended in fslio\n");
	unitcode = NIFTI_UNITS_USEC;
    }
    fslio->niftiptr->time_units = unitcode;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetTimeUnits(FSLIO *fslio, char *units)
{
  if (fslio==NULL)  FSLIOERR("FslGetTimeUnits: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    strcpy(units,nifti_units_string(fslio->niftiptr->time_units));
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslSetDataType(FSLIO *fslio, short t)
{
  int nbytepix=0, ss=0;
  if (fslio==NULL)  FSLIOERR("FslSetDataType: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    fslio->niftiptr->datatype = t;
    nifti_datatype_sizes(t,&nbytepix,&ss);
    fslio->niftiptr->nbyper = nbytepix;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}

size_t FslGetDataType(FSLIO *fslio, short *t)
{
    /* returns bits per pixel */
  int nbytepix=32, ss=0;
  if (fslio==NULL)  FSLIOERR("FslGetDataType: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *t = fslio->niftiptr->datatype;
    nifti_datatype_sizes(*t,&nbytepix,&ss);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return (size_t) 8 * nbytepix;
}


void FslGetMMCoord(mat44 stdmat, float voxx, float voxy, float voxz, 
		   float *mmx, float *mmy, float *mmz) 
{
    *mmx = stdmat.m[0][0] * voxx + stdmat.m[0][1] * voxy + stdmat.m[0][2] * voxz 
	+ stdmat.m[0][3];
    *mmy = stdmat.m[1][0] * voxx + stdmat.m[1][1] * voxy + stdmat.m[1][2] * voxz 
	+ stdmat.m[1][3];
    *mmz = stdmat.m[2][0] * voxx + stdmat.m[2][1] * voxy + stdmat.m[2][2] * voxz 
	+ stdmat.m[2][3];
}


void FslGetVoxCoord(mat44 stdmat, float mmx, float mmy, float mmz, 
		   float *voxx, float *voxy, float *voxz) 
{
  mat44 mm2vox;

  mm2vox = mat44_inverse(stdmat);
    *voxx = mm2vox.m[0][0] * mmx + mm2vox.m[0][1] * mmy + mm2vox.m[0][2] * mmz 
	+ mm2vox.m[0][3];
    *voxy = mm2vox.m[1][0] * mmx + mm2vox.m[1][1] * mmy + mm2vox.m[1][2] * mmz 
	+ mm2vox.m[1][3];
    *voxz = mm2vox.m[2][0] * mmx + mm2vox.m[2][1] * mmy + mm2vox.m[2][2] * mmz 
	+ mm2vox.m[2][3];
}


void FslSetStdXform(FSLIO *fslio, short sform_code, mat44 stdmat)
{
    /* NB: stdmat must point to a 4x4 array */
  if (fslio==NULL)  FSLIOERR("FslSetStdXform: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      fslio->niftiptr->sform_code = sform_code;
      fslio->niftiptr->sto_xyz.m[0][0] = stdmat.m[0][0];
      fslio->niftiptr->sto_xyz.m[0][1] = stdmat.m[0][1];
      fslio->niftiptr->sto_xyz.m[0][2] = stdmat.m[0][2];
      fslio->niftiptr->sto_xyz.m[0][3] = stdmat.m[0][3];
      fslio->niftiptr->sto_xyz.m[1][0] = stdmat.m[1][0];
      fslio->niftiptr->sto_xyz.m[1][1] = stdmat.m[1][1];
      fslio->niftiptr->sto_xyz.m[1][2] = stdmat.m[1][2];
      fslio->niftiptr->sto_xyz.m[1][3] = stdmat.m[1][3];
      fslio->niftiptr->sto_xyz.m[2][0] = stdmat.m[2][0];
      fslio->niftiptr->sto_xyz.m[2][1] = stdmat.m[2][1];
      fslio->niftiptr->sto_xyz.m[2][2] = stdmat.m[2][2];
      fslio->niftiptr->sto_xyz.m[2][3] = stdmat.m[2][3];
      fslio->niftiptr->sto_xyz.m[3][0] = 0;
      fslio->niftiptr->sto_xyz.m[3][1] = 0;
      fslio->niftiptr->sto_xyz.m[3][2] = 0;
      fslio->niftiptr->sto_xyz.m[3][3] = 1;
      fslio->niftiptr->sto_ijk = mat44_inverse(fslio->niftiptr->sto_xyz);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


short FslGetStdXform(FSLIO *fslio, mat44 *stdmat)
{
    /* returns sform code  (NB: stdmat must point to a 4x4 array) */
    float dx,dy,dz,tr;
    if (fslio==NULL)  FSLIOERR("FslGetStdXform: Null pointer passed for FSLIO");
    if (fslio->niftiptr!=NULL) {
	stdmat->m[0][0] = fslio->niftiptr->sto_xyz.m[0][0];
	stdmat->m[0][1] = fslio->niftiptr->sto_xyz.m[0][1];
	stdmat->m[0][2] = fslio->niftiptr->sto_xyz.m[0][2];
	stdmat->m[0][3] = fslio->niftiptr->sto_xyz.m[0][3];
	stdmat->m[1][0] = fslio->niftiptr->sto_xyz.m[1][0];
	stdmat->m[1][1] = fslio->niftiptr->sto_xyz.m[1][1];
	stdmat->m[1][2] = fslio->niftiptr->sto_xyz.m[1][2];
	stdmat->m[1][3] = fslio->niftiptr->sto_xyz.m[1][3];
	stdmat->m[2][0] = fslio->niftiptr->sto_xyz.m[2][0];
	stdmat->m[2][1] = fslio->niftiptr->sto_xyz.m[2][1];
	stdmat->m[2][2] = fslio->niftiptr->sto_xyz.m[2][2];
	stdmat->m[2][3] = fslio->niftiptr->sto_xyz.m[2][3];
	stdmat->m[3][0] = 0.0;
	stdmat->m[3][1] = 0.0;
	stdmat->m[3][2] = 0.0;
	stdmat->m[3][3] = 1.0;
	
	/* the code below gives a default but it really should never be used */
        if (fslio->niftiptr->sform_code == NIFTI_XFORM_UNKNOWN) {
	    FslGetVoxDim(fslio,&dx,&dy,&dz,&tr);
	    stdmat->m[0][0] = -dx;  /* default Radiological convention */
	    stdmat->m[0][1] = 0;
	    stdmat->m[0][2] = 0;
	    stdmat->m[0][3] = 0;
	    stdmat->m[1][0] = 0;
	    stdmat->m[1][1] = dy;
	    stdmat->m[1][2] = 0;
	    stdmat->m[1][3] = 0;
	    stdmat->m[2][0] = 0;
	    stdmat->m[2][1] = 0;
	    stdmat->m[2][2] = dz;
	    stdmat->m[2][3] = 0;
	    stdmat->m[3][0] = 0.0;
	    stdmat->m[3][1] = 0.0;
	    stdmat->m[3][2] = 0.0;
	    stdmat->m[3][3] = 1.0;
	}
	return fslio->niftiptr->sform_code;
    }
    if (fslio->mincptr!=NULL) {
	fprintf(stderr,"Warning:: Minc is not yet supported\n");
    }
    return NIFTI_XFORM_UNKNOWN;
}


void FslSetRigidXform(FSLIO *fslio, short qform_code, mat44 rigidmat)
{
    /* NB: rigidmat must point to an allocated mat44 */
  float dx, dy, dz;
  if (fslio==NULL)  FSLIOERR("FslSetRigidXform: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      fslio->niftiptr->qform_code = qform_code;
      fslio->niftiptr->qto_xyz.m[0][0] = rigidmat.m[0][0];
      fslio->niftiptr->qto_xyz.m[0][1] = rigidmat.m[0][1];
      fslio->niftiptr->qto_xyz.m[0][2] = rigidmat.m[0][2];
      fslio->niftiptr->qto_xyz.m[0][3] = rigidmat.m[0][3];
      fslio->niftiptr->qto_xyz.m[1][0] = rigidmat.m[1][0];
      fslio->niftiptr->qto_xyz.m[1][1] = rigidmat.m[1][1];
      fslio->niftiptr->qto_xyz.m[1][2] = rigidmat.m[1][2];
      fslio->niftiptr->qto_xyz.m[1][3] = rigidmat.m[1][3];
      fslio->niftiptr->qto_xyz.m[2][0] = rigidmat.m[2][0];
      fslio->niftiptr->qto_xyz.m[2][1] = rigidmat.m[2][1];
      fslio->niftiptr->qto_xyz.m[2][2] = rigidmat.m[2][2];
      fslio->niftiptr->qto_xyz.m[2][3] = rigidmat.m[2][3];
      fslio->niftiptr->qto_xyz.m[3][0] = 0;
      fslio->niftiptr->qto_xyz.m[3][1] = 0;
      fslio->niftiptr->qto_xyz.m[3][2] = 0;
      fslio->niftiptr->qto_xyz.m[3][3] = 1;
      mat44_to_quatern(fslio->niftiptr->qto_xyz,&(fslio->niftiptr->quatern_b),
		       &(fslio->niftiptr->quatern_c),&(fslio->niftiptr->quatern_d),
		       &(fslio->niftiptr->qoffset_x),&(fslio->niftiptr->qoffset_y),
		       &(fslio->niftiptr->qoffset_z),&dx,&dy,&dz,&(fslio->niftiptr->qfac));
      fslio->niftiptr->qto_ijk = mat44_inverse(fslio->niftiptr->qto_xyz);

  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


short FslGetRigidXform(FSLIO *fslio, mat44 *rigidmat)
{
    /* returns qform code  (NB: rigidmat must point to an allocated mat44) */
    float dx,dy,dz,tr;
    if (fslio==NULL)  FSLIOERR("FslGetRigidXform: Null pointer passed for FSLIO");
    if (fslio->niftiptr!=NULL) {
	rigidmat->m[0][0] = fslio->niftiptr->qto_xyz.m[0][0];
	rigidmat->m[0][1] = fslio->niftiptr->qto_xyz.m[0][1];
	rigidmat->m[0][2] = fslio->niftiptr->qto_xyz.m[0][2];
	rigidmat->m[0][3] = fslio->niftiptr->qto_xyz.m[0][3];
	rigidmat->m[1][0] = fslio->niftiptr->qto_xyz.m[1][0];
	rigidmat->m[1][1] = fslio->niftiptr->qto_xyz.m[1][1];
	rigidmat->m[1][2] = fslio->niftiptr->qto_xyz.m[1][2];
	rigidmat->m[1][3] = fslio->niftiptr->qto_xyz.m[1][3];
	rigidmat->m[2][0] = fslio->niftiptr->qto_xyz.m[2][0];
	rigidmat->m[2][1] = fslio->niftiptr->qto_xyz.m[2][1];
	rigidmat->m[2][2] = fslio->niftiptr->qto_xyz.m[2][2];
	rigidmat->m[2][3] = fslio->niftiptr->qto_xyz.m[2][3];
	rigidmat->m[3][0] = 0.0;
	rigidmat->m[3][1] = 0.0;
	rigidmat->m[3][2] = 0.0;
	rigidmat->m[3][3] = 1.0;
	
	/* the code gives a default but it should never really be used */
        if (fslio->niftiptr->qform_code == NIFTI_XFORM_UNKNOWN) {
	  FslGetVoxDim(fslio,&dx,&dy,&dz,&tr);
	  rigidmat->m[0][0] = -dx;  /* default Radiological convention */
	  rigidmat->m[0][1] = 0;
	  rigidmat->m[0][2] = 0;
	  rigidmat->m[0][3] = 0;
	  rigidmat->m[1][0] = 0;
	  rigidmat->m[1][1] = dy;
	  rigidmat->m[1][2] = 0;
	  rigidmat->m[1][3] = 0;
	  rigidmat->m[2][0] = 0;
	  rigidmat->m[2][1] = 0;
	  rigidmat->m[2][2] = dz;
	  rigidmat->m[3][0] = 0.0;
	  rigidmat->m[3][1] = 0.0;
	  rigidmat->m[3][2] = 0.0;
	  rigidmat->m[3][3] = 1.0;
	}
	return fslio->niftiptr->qform_code;
    }
    if (fslio->mincptr!=NULL) {
	fprintf(stderr,"Warning:: Minc is not yet supported\n");
    }
    return NIFTI_XFORM_UNKNOWN;
}


void FslSetIntent(FSLIO *fslio, short intent_code, float p1, float p2, float p3)
{
  if (fslio==NULL)  FSLIOERR("FslSetIntent: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      fslio->niftiptr->intent_code = intent_code;
      fslio->niftiptr->intent_p1 = p1;
      fslio->niftiptr->intent_p2 = p2;
      fslio->niftiptr->intent_p3 = p3;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


short FslGetIntent(FSLIO *fslio, short *intent_code, float *p1, float *p2,
		   float *p3)
{
  /* also returns intent code */
  if (fslio==NULL)  FSLIOERR("FslGetIntent: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      *intent_code = fslio->niftiptr->intent_code;
      *p1 = fslio->niftiptr->intent_p1;
      *p2 = fslio->niftiptr->intent_p2;
      *p3 = fslio->niftiptr->intent_p3;
      return fslio->niftiptr->intent_code;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return NIFTI_INTENT_NONE;
}




void FslSetIntensityScaling(FSLIO *fslio, float slope, float intercept)
{
  if (fslio==NULL)  FSLIOERR("FslSetIntensityScaling: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
      fslio->niftiptr->scl_slope = slope;
      fslio->niftiptr->scl_inter = intercept;
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


int FslGetIntensityScaling(FSLIO *fslio, float *slope, float *intercept)
{
  /* returns 1 if scaling required or 0 otherwise */
  if (fslio==NULL)  FSLIOERR("FslGetIntensityScaling: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    *slope = fslio->niftiptr->scl_slope;
    *intercept = fslio->niftiptr->scl_inter;
    if (fabs(*slope)<1e-30) {
      *slope = 1.0;
      *intercept = 0.0;
      return 0;
    }
    if ( (fabs(*slope - 1.0)>1e-30) || (fabs(*intercept)>1e-30) ) {
      return 1;
    } else {
      return 0;
    }
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return 0;
 
}


mat33 mat44_to_mat33(mat44 x)
{
  mat33 y;
  int i,j;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      y.m[i][j] = x.m[i][j];
    }
  }
  return y;
}


int FslGetLeftRightOrder2(int sform_code, mat44 sform44, 
			  int qform_code, mat44 qform44)
{
  /* Determines if the image is stored in neurological or radiological convention */
  int order=FSL_RADIOLOGICAL;
  float dets=-1.0, detq=-1.0, det=-1.0;
  mat33 sform33, qform33;
  if (qform_code!=NIFTI_XFORM_UNKNOWN) { 
    qform33 = mat44_to_mat33(qform44);
    detq = mat33_determ(qform33);
    det = detq;
  }
  if (sform_code!=NIFTI_XFORM_UNKNOWN) { 
    sform33 = mat44_to_mat33(sform44);
    dets = mat33_determ(sform33);
    det = dets;
  }
  
  if (det<0.0) order=FSL_RADIOLOGICAL;
  else order=FSL_NEUROLOGICAL;
  /* check for inconsistency if both are set */
  if ( (sform_code!=NIFTI_XFORM_UNKNOWN) && 
       (qform_code!=NIFTI_XFORM_UNKNOWN) ) { 
    if (dets * detq < 0.0) order=FSL_INCONSISTENT;
  }
  return order;
}


int FslGetLeftRightOrder(FSLIO *fslio)
{
  int order=FSL_RADIOLOGICAL, sform_code, qform_code;
  mat44 sform44, qform44;
  if (fslio==NULL)  FSLIOERR("FslGetLeftRightOrder: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    sform_code = FslGetStdXform(fslio,&sform44);
    qform_code = FslGetRigidXform(fslio,&qform44);
    return FslGetLeftRightOrder2(sform_code,sform44,qform_code,qform44);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return order;
}

short FslGetVox2mmMatrix2(mat44 *vox2mm, int sform_code, mat44 sform44, 
			  int qform_code, mat44 qform44, 
			  float dx, float dy, float dz)
{
  short retcode=NIFTI_XFORM_UNKNOWN;
  int ii,jj;
  if (sform_code!=NIFTI_XFORM_UNKNOWN) { 
    for (ii=0; ii<4; ii++) { 
      for (jj=0; jj<4; jj++) {
	vox2mm->m[ii][jj] = sform44.m[ii][jj];
      } 
    }
    retcode=sform_code;
  } else if (qform_code!=NIFTI_XFORM_UNKNOWN) { 
    for (ii=0; ii<4; ii++) { 
      for (jj=0; jj<4; jj++) {
	vox2mm->m[ii][jj] = qform44.m[ii][jj];
      } 
    }
    retcode=qform_code;
  } else {
    /* default case - for FSLView is positive voxel to mm scalings */
    vox2mm->m[0][0] = dx;
    vox2mm->m[0][1] = 0.0;
    vox2mm->m[0][2] = 0.0;
    vox2mm->m[0][3] = 0.0;
    vox2mm->m[1][0] = 0.0;
    vox2mm->m[1][1] = dy;
    vox2mm->m[1][2] = 0.0;
    vox2mm->m[1][3] = 0.0;
    vox2mm->m[2][0] = 0.0;
    vox2mm->m[2][1] = 0.0;
    vox2mm->m[2][2] = dz;
    vox2mm->m[3][0] = 0.0;
    vox2mm->m[3][1] = 0.0;
    vox2mm->m[3][2] = 0.0;
    vox2mm->m[3][3] = 1.0;
    retcode=NIFTI_XFORM_UNKNOWN;
  }
  return retcode;
}


short FslGetVox2mmMatrix(FSLIO *fslio, mat44 *vox2mm)
{
  int sform_code, qform_code;
  float dx, dy, dz, tr;
  mat44 sform44, qform44;
  if (fslio==NULL)  FSLIOERR("FslGetVox2mmMatrix: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    sform_code = FslGetStdXform(fslio,&sform44);
    qform_code = FslGetRigidXform(fslio,&qform44);
    FslGetVoxDim(fslio,&dx,&dy,&dz,&tr);
    return FslGetVox2mmMatrix2(vox2mm,sform_code,sform44,
			       qform_code,qform44,dx,dy,dz);
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
  return NIFTI_XFORM_UNKNOWN;
}

void FslSetAnalyzeSform(FSLIO *fslio, const short *orig,
			float dx, float dy, float dz)
{
  /* Creates an sform matrix for an Analyze file */
  /* THIS ALWAYS CREATES A RADIOLOGICAL ORDERED SFORM */
  /* NB: the origin passed in here is in Analyze convention - starting at 1, not 0 */
  float x, y, z;
  if (fslio==NULL)  FSLIOERR("FslSetAnalyzeSform: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    if (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE) {
      /* default case */
      fslio->niftiptr->sform_code = NIFTI_XFORM_UNKNOWN;
    }
    /* ignore all zero origins - really all serious coord stuff should
       be done via the FslSetStdCoord call */
    if ((orig[0]!=0) || (orig[1]!=0) || (orig[2]!=0))
      {
	short origx=0, origy=0, origz=0;
	if ((orig[0]!=0) || (orig[1]!=0) || (orig[2]!=0)) {
	  /* convert to nifti conventions (start at 0 not 1) */
	  origx = orig[0] - 1;
	  origy = orig[1] - 1;
	  origz = orig[2] - 1;
	}
	if ( dx * dy * dz > 0 ) {
	  /* change neurological convention to radiological if necessary */
	  dx = -dx;
	}
	if ( (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE) 
	     || (fslio->niftiptr->sform_code == NIFTI_XFORM_UNKNOWN) ) {
	  /* make a default transform with the requested origin at xyz=000 */ 
	  fslio->niftiptr->sform_code = NIFTI_XFORM_ALIGNED_ANAT;
	  fslio->niftiptr->sto_xyz.m[0][0] = dx;
	  fslio->niftiptr->sto_xyz.m[0][1] = 0;
	  fslio->niftiptr->sto_xyz.m[0][2] = 0;
	  fslio->niftiptr->sto_xyz.m[0][3] = -(origx)*(dx);
	  fslio->niftiptr->sto_xyz.m[1][0] = 0;
	  fslio->niftiptr->sto_xyz.m[1][1] = dy;
	  fslio->niftiptr->sto_xyz.m[1][2] = 0;
	  fslio->niftiptr->sto_xyz.m[1][3] = -(origy)*(dy);
	  fslio->niftiptr->sto_xyz.m[2][0] = 0;
	  fslio->niftiptr->sto_xyz.m[2][1] = 0;
	  fslio->niftiptr->sto_xyz.m[2][2] = dz;
	  fslio->niftiptr->sto_xyz.m[2][3] = -(origz)*(dz);
	  fslio->niftiptr->sto_xyz.m[3][0] = 0;
	  fslio->niftiptr->sto_xyz.m[3][1] = 0;
	  fslio->niftiptr->sto_xyz.m[3][2] = 0;
	  fslio->niftiptr->sto_xyz.m[3][3] = 1;
	  fslio->niftiptr->sto_ijk = mat44_inverse(fslio->niftiptr->sto_xyz);
	} else {
	  /* update the existing origin */
	  /* find out what the existing xyz of the requested origin is */
	  x = fslio->niftiptr->sto_xyz.m[0][0] * origx
	    + fslio->niftiptr->sto_xyz.m[0][1] * origy
	    + fslio->niftiptr->sto_xyz.m[0][2] * origz
	    + fslio->niftiptr->sto_xyz.m[0][3];
	  y = fslio->niftiptr->sto_xyz.m[1][0] * origx
	    + fslio->niftiptr->sto_xyz.m[1][1] * origy
	    + fslio->niftiptr->sto_xyz.m[1][2] * origz
	    + fslio->niftiptr->sto_xyz.m[1][3];
	  z = fslio->niftiptr->sto_xyz.m[2][0] * origx
	    + fslio->niftiptr->sto_xyz.m[2][1] * origy
	    + fslio->niftiptr->sto_xyz.m[2][2] * origz
	    + fslio->niftiptr->sto_xyz.m[2][3];
	  /* subtract off whatever is currently the xyz of the origin */
	  fslio->niftiptr->sto_xyz.m[0][3] -= x;
	  fslio->niftiptr->sto_xyz.m[1][3] -= y;
	  fslio->niftiptr->sto_xyz.m[2][3] -= z;
	  fslio->niftiptr->sto_ijk = mat44_inverse(fslio->niftiptr->sto_xyz);
	}
	
      }
    
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}


void FslGetAnalyzeOrigin(FSLIO *fslio, short orig[5])
{
  /* NB: orig returned here is in Analyze convention - starting at 1, not 0 */
  if (fslio==NULL)  FSLIOERR("FslGetAnalyzeOrigin: Null pointer passed for FSLIO");
  if (fslio->niftiptr!=NULL) {
    /* Use sform or qform to determine the origin - default is zero */
      orig[0]=0; 
      orig[1]=0; 
      orig[2]=0; 
      orig[3]=0; 
      orig[4]=0;

      if (fslio->niftiptr->qform_code != NIFTI_XFORM_UNKNOWN) {
	orig[0]=(short) fslio->niftiptr->qto_ijk.m[0][3] + 1;
	orig[1]=(short) fslio->niftiptr->qto_ijk.m[1][3] + 1;
	orig[2]=(short) fslio->niftiptr->qto_ijk.m[2][3] + 1;
      } 

      if (fslio->niftiptr->sform_code != NIFTI_XFORM_UNKNOWN) {
	orig[0]=(short) fslio->niftiptr->sto_ijk.m[0][3] + 1;
	orig[1]=(short) fslio->niftiptr->sto_ijk.m[1][3] + 1;
	orig[2]=(short) fslio->niftiptr->sto_ijk.m[2][3] + 1;
      } 
  }
  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
  }
}



int FslClose(FSLIO *fslio)
{
  int retval=0, filetype;
  struct dsr *hdr;
  znzFile hptr=NULL;

  if (fslio==NULL)   return 0;

  /* close the (data) file */
  if (!znz_isnull(fslio->fileptr)) retval=znzclose(fslio->fileptr);

  /** ----- if writing the image, need to worry about the header bit ----- **/

  if ( (fslio->niftiptr!=NULL) && (FslGetWriteMode(fslio)==1) 
       && (fslio->written_hdr==0) ) {

    /* ensure that the type is set correctly */
    fslio->niftiptr->nifti_type = FslBaseFileType(FslGetFileType(fslio));

    /* must write the header now */
    filetype = FslGetFileType(fslio);
    strcpy(fslio->niftiptr->descrip,"FSL3.3");
    if (!FslIsSingleFileType(filetype)) {
      /* for file pairs - open new header file and write it */
      nifti_image_write_hdr_img(fslio->niftiptr,0,"wb");
    } else {
      /* for single files it is more complicated */
      if (!FslIsCompressedFileType(filetype)) {
	/* noncompressed -> reopen this file in r+ mode and write the header part again */
	nifti_image_write_hdr_img(fslio->niftiptr,0,"r+b");
      } else {
	/* compressed mode -> not possible! */
	fprintf(stderr,"Error:: header must be written before writing any other data.\n");
	return -1;
      }
    }
  }
    
  /* --- nasty hack to write the origin in Analyze files --- */

  if ( (FslGetWriteMode(fslio)==1) && (fslio->niftiptr!=NULL) && 
       (FslBaseFileType(FslGetFileType(fslio))==FSL_TYPE_ANALYZE) ) {
 
    /* read in the old header, change the origin and write it out again */
    hdr = (struct dsr *) calloc(1,sizeof(struct dsr));
    FslReadRawHeader(hdr,fslio->niftiptr->fname);
    if (fslio->niftiptr->byteorder != short_order()) { AvwSwapHeader(hdr); }
    
    /* calculate origin from sform (if set) */
    {
      short blah[5];
      FslGetAnalyzeOrigin(fslio,blah);
      memcpy(hdr->hist.originator,blah,5*sizeof(short));
    
      /* Write out in radiological order if origin is non-zero */
      /* set negative pixdim if needed to keep LR orientation consistent */
      if ( (blah[0]!=0) || (blah[1]!=0) || (blah[2]!=0) ) {
        if (hdr->dime.pixdim[1] * hdr->dime.pixdim[2] * hdr->dime.pixdim[3] > 0) {
          hdr->dime.pixdim[1] = - hdr->dime.pixdim[1]; 
        }
      }
    }

    /* swap back byte order and write out */
    if (fslio->niftiptr->byteorder != short_order()) { AvwSwapHeader(hdr); }
    hptr = znzopen(fslio->niftiptr->fname,"wb",FslIsCompressedFileType(FslGetFileType(fslio)));
    if (znz_isnull(hptr)) { 	
      fprintf(stderr,"Error:: Could not write origin data to header file %s.\n",
	      fslio->niftiptr->fname);
      return -1;
    };
    
    znzwrite(hdr,1,sizeof(struct dsr),hptr);
    znzclose(hptr);
    free(hdr);
  }

  if (fslio->mincptr!=NULL) {
    fprintf(stderr,"Warning:: Minc is not yet supported\n");
    return -1;
  }

  return retval;
}
  

void AvwSwapHeader(struct dsr *avw)
{
  char *ptr;

  ptr = (char *) &(avw->hk);
  swap_4bytes(1,ptr);		/* sizeof_hdr */
  ptr += 32;
  swap_4bytes(1,ptr);		/* extents */
  ptr += 4;
  swap_2bytes(1,ptr);		/* session_error */
  
  ptr = (char *) &(avw->dime);
  swap_2bytes(8,ptr);		/* dims */
  ptr += 28;
  swap_2bytes(4,ptr);		/* unused1, datatype, bitpix, dim_un0 */
  ptr += 8;
  swap_4bytes(18,ptr);		/* pixdim, vox_offset, ... */
				/* cal_min, compressed, ... glmin */

  ptr = (char *) &(avw->hist);
  ptr += 105;
  swap_2bytes(5,ptr);		/* originator (used to store origin) */
  ptr += 63;
  swap_4bytes(8,ptr);		/* views, ... smin */
}


int FslReadRawHeader(void *buffer, const char* filename)
{
  znzFile fp;
  int retval;
  fp = znzopen(filename,"rb",1);
  if (znz_isnull(fp)) {
    znzclose(fp);
    fprintf(stderr,"Could not open header %s\n",filename);
    return 0;
  }
  retval = znzread(buffer,1,348,fp);
  if (retval != 348) {
    znzclose(fp);
    fprintf(stderr,"Could not read header %s\n",filename);
    return retval;
  }
  znzclose(fp);
  return retval;
}

void FslSetOverrideOutputType(int type)
{
  if ( (type==-1) || (FslIsValidFileType(type)) ) {
    FslOverrideOutputType=type;
  } else {
    fprintf(stderr,"Invalid file type (%d) requested - ignoring this\n",type);
  }
}

int FslGetOverrideOutputType()
{
  return FslOverrideOutputType;
}

void FslSetIgnoreMFQ(int flag)
{
  assert((flag==0) || (flag==1));
  FslIgnoreMFQ=flag;
}


int FslGetIgnoreMFQ()
{
  return FslIgnoreMFQ;
}

