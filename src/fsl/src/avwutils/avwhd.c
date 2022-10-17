/* {{{ copyright */

/*  avwhd.c - show image header

    Steve Smith and Mark Jenkinson, FMRIB Image Analysis Group
    With memset(string) suggestion from Bettyann Chodkowski

    Copyright (C) 2000-2004 University of Oxford  */

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

/* }}} */
/* {{{ defines, includes and typedefs */

#include "fslio/fslio.h"
#include "fslio/dbh.h"

void usage(void);

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("\nUsage: avwhd [-x] <input>\n");
  printf("   -x : instead print an XML-style NIFTI header\n");
  exit(1);
}

/* }}} */

void ShowNifti(char *fileName, FSLIO* fslio)
{
    if (fslio == NULL) {
	fprintf(stderr,"ERROR: Could not open file\n");
	return;
    }
    if (fslio->niftiptr!=NULL) {
	printf("%s\n",nifti_image_to_ascii(fslio->niftiptr));
	return;
    }
    if (fslio->mincptr!=NULL) {
	fprintf(stderr,"ERROR: Minc is not currently supported\n");
	return;
    }
    return;
}

/* {{{ ShowHdr */

void ShowHdr(char *fileName, FSLIO* fslio)
{
  int i, ft, isanalyze=0;
  char string[128];
  struct dsr *hdr;
  mat44 mat;
  int icode, jcode, kcode;

  if (fslio == NULL) {
    fprintf(stderr,"ERROR: Could not open file\n");
    return;
  }

  ft = FslGetFileType(fslio);
  if (FslBaseFileType(ft)==FSL_TYPE_MINC) {
    fprintf(stderr,"ERROR: Minc is not currently supported\n");
    return;
  }

  if (fslio->niftiptr == NULL) {
    fprintf(stderr,"ERROR: Not an Analyze or Nifti file\n");
    return;
  }


  /* -------------------- ANALYZE CASE ----------------------- */


  if (FslBaseFileType(ft)==FSL_TYPE_ANALYZE) { 
    isanalyze=1; 
    /* load raw hdr structure */
    
    hdr = (struct dsr *)calloc(1,sizeof(struct dsr));
    FslReadRawHeader(hdr,fslio->niftiptr->fname);

    if (fslio->niftiptr->byteorder != short_order()) {
      printf("Byte swapping\n");
      AvwSwapHeader(hdr);
    }

    printf("filename       %s\n\n", fileName);
    
    /* Header Key */
    printf("sizeof_hdr     %d\n", hdr->hk.sizeof_hdr);
    printf("data_type      %s\n", hdr->hk.data_type);
    printf("db_name        %s\n", hdr->hk.db_name);
    printf("extents        %d\n", hdr->hk.extents);
    printf("session_error  %d\n", hdr->hk.session_error);
    printf("regular        %c\n", hdr->hk.regular);
    printf("hkey_un0       %c\n", hdr->hk.hkey_un0);
    
    /* Image Dimension */
    for(i=0;i<8;i++)
      printf("dim%d           %d\n", i, hdr->dime.dim[i]);
    
    memset(string,0,128);
    strncpy(string,hdr->dime.vox_units,4);
    printf("vox_units      %s\n", string);
    
    memset(string,0,128);
    strncpy(string,hdr->dime.cal_units,8);
    printf("cal_units      %s\n", string);
    printf("unused1        %d\n", hdr->dime.unused1);
    printf("datatype       %d\n", hdr->dime.datatype);
    printf("bitpix         %d\n", hdr->dime.bitpix);
    
    for(i=0;i<8;i++)
      printf("pixdim%d        %.10f\n",i, hdr->dime.pixdim[i]);
    
    printf("vox_offset     %6.4f\n",  hdr->dime.vox_offset);
    printf("funused1       %6.4f\n", hdr->dime.funused1);
    printf("funused2       %6.4f\n", hdr->dime.funused2);
    printf("funused3       %6.4f\n", hdr->dime.funused3);
    printf("cal_max        %6.4f\n", hdr->dime.cal_max);
    printf("cal_min        %6.4f\n", hdr->dime.cal_min);
    printf("compressed     %d\n", hdr->dime.compressed);
    printf("verified       %d\n", hdr->dime.verified);
    printf("glmax          %d\n", hdr->dime.glmax);
    printf("glmin          %d\n", hdr->dime.glmin);
    
    /* Data History */
    memset(string,0,128);
    strncpy(string,hdr->hist.descrip,80);
    printf("descrip        %s\n", string);
    memset(string,0,128);
    strncpy(string,hdr->hist.aux_file,24);
    printf("aux_file       %s\n", string);
    printf("orient         %d\n", hdr->hist.orient);
    
    memset(string,0,128);
    strncpy(string,hdr->hist.originator,10);
    printf("originator     %s\n", string);
    {
      short blah[5];
      memcpy(blah,hdr->hist.originator,5*sizeof(short));
      printf("origin1        %d\n",blah[0]);
      printf("origin2        %d\n",blah[1]);
      printf("origin3        %d\n",blah[2]);
    }
    
    memset(string,0,128);
    strncpy(string,hdr->hist.generated,10);
    printf("generated      %s\n", string);
    
    memset(string,0,128);
    strncpy(string,hdr->hist.scannum,10);
    printf("scannum        %s\n", string);
    
    memset(string,0,128);
    strncpy(string,hdr->hist.patient_id,10);
    printf("patient_id     %s\n", string);
    
    memset(string,0,128);
    strncpy(string,hdr->hist.exp_date,10);
    printf("exp_date       %s\n", string);
    
    memset(string,0,128);
    strncpy(string,hdr->hist.exp_time,10);
    printf("exp_time       %s\n", string);
    
    memset(string,0,128);
    strncpy(string,hdr->hist.hist_un0,10);
    printf("hist_un0       %s\n", string);
    
    printf("views          %d\n", hdr->hist.views);
    printf("vols_added     %d\n", hdr->hist.vols_added);
    printf("start_field    %d\n", hdr->hist.start_field);
    printf("field_skip     %d\n", hdr->hist.field_skip);
    printf("omax           %d\n", hdr->hist.omax);
    printf("omin           %d\n", hdr->hist.omin);
    printf("smin           %d\n", hdr->hist.smax);
    printf("smin           %d\n", hdr->hist.smin);

    printf("file_type      %s\n",FslFileTypeString(0));
    printf("file_code      0\n");

    return;
  }
  
  
  /* -------------------- NIFTI CASE ----------------------- */

  if (fslio->niftiptr->byteorder != short_order()) {
    printf("Byte swapping\n");
  }

  printf("filename       %s\n\n", fslio->niftiptr->fname);

  printf("sizeof_hdr     %d\n", 348);
  printf("data_type      %s\n", nifti_datatype_string(fslio->niftiptr->datatype));
  for(i=0;i<8;i++)
    printf("dim%d           %d\n", i, fslio->niftiptr->dim[i]);
  
  printf("vox_units      %s\n", nifti_units_string(fslio->niftiptr->xyz_units));
  printf("time_units     %s\n", nifti_units_string(fslio->niftiptr->time_units));
	
  printf("datatype       %d\n", fslio->niftiptr->datatype);
  printf("nbyper         %d\n", fslio->niftiptr->nbyper);
  printf("bitpix         %d\n", fslio->niftiptr->nbyper * 8);
  for(i=0;i<8;i++)
    printf("pixdim%d        %.10f\n",i, fslio->niftiptr->pixdim[i]);

  printf("vox_offset     %d\n",  fslio->niftiptr->iname_offset);

  printf("cal_max        %6.4f\n", fslio->niftiptr->cal_max);
  printf("cal_min        %6.4f\n", fslio->niftiptr->cal_min);

  printf("scl_slope      %f\n", fslio->niftiptr->scl_slope);
  printf("scl_inter      %f\n", fslio->niftiptr->scl_inter);

  printf("phase_dim      %d\n", fslio->niftiptr->phase_dim);
  printf("freq_dim       %d\n", fslio->niftiptr->freq_dim);
  printf("slice_dim      %d\n", fslio->niftiptr->slice_dim);

  printf("slice_name     %s\n", nifti_slice_string(fslio->niftiptr->slice_code));
  printf("slice_code     %d\n", fslio->niftiptr->slice_code);
  printf("slice_start    %d\n", fslio->niftiptr->slice_start);
  printf("slice_end      %d\n", fslio->niftiptr->slice_end);
  printf("slice_duration %f\n", fslio->niftiptr->slice_duration);

  printf("time_offset    %f\n", fslio->niftiptr->toffset);

  printf("intent         %s\n", nifti_intent_string(fslio->niftiptr->intent_code));
  printf("intent_code    %d\n", fslio->niftiptr->intent_code);
  printf("intent_name    %s\n", fslio->niftiptr->intent_name);
  printf("intent_p1      %f\n", fslio->niftiptr->intent_p1);
  printf("intent_p2      %f\n", fslio->niftiptr->intent_p2);
  printf("intent_p3      %f\n", fslio->niftiptr->intent_p3);

  printf("qform_name     %s\n", nifti_xform_string(fslio->niftiptr->qform_code));
  printf("qform_code     %d\n", fslio->niftiptr->qform_code);

  mat = fslio->niftiptr->qto_xyz;
  printf("qto_xyz:1      %f  %f  %f  %f\n",mat.m[0][0],mat.m[0][1],
	 mat.m[0][2],mat.m[0][3]);
  printf("qto_xyz:2      %f  %f  %f  %f\n",mat.m[1][0],mat.m[1][1],
	 mat.m[1][2],mat.m[1][3]);
  printf("qto_xyz:3      %f  %f  %f  %f\n",mat.m[2][0],mat.m[2][1],
	 mat.m[2][2],mat.m[2][3]);
  printf("qto_xyz:4      %f  %f  %f  %f\n",mat.m[3][0],mat.m[3][1],
	 mat.m[3][2],mat.m[3][3]);

  mat44_to_orientation(mat,&icode,&jcode,&kcode);
  printf("qform_xorient  %s\n",nifti_orientation_string(icode));
  printf("qform_yorient  %s\n",nifti_orientation_string(jcode));
  printf("qform_zorient  %s\n",nifti_orientation_string(kcode));


  printf("sform_name     %s\n", nifti_xform_string(fslio->niftiptr->sform_code));
  printf("sform_code     %d\n", fslio->niftiptr->sform_code);

  mat = fslio->niftiptr->sto_xyz;
  printf("sto_xyz:1      %f  %f  %f  %f\n",mat.m[0][0],mat.m[0][1],
	 mat.m[0][2],mat.m[0][3]);
  printf("sto_xyz:2      %f  %f  %f  %f\n",mat.m[1][0],mat.m[1][1],
	 mat.m[1][2],mat.m[1][3]);
  printf("sto_xyz:3      %f  %f  %f  %f\n",mat.m[2][0],mat.m[2][1],
	 mat.m[2][2],mat.m[2][3]);
  printf("sto_xyz:4      %f  %f  %f  %f\n",mat.m[3][0],mat.m[3][1],
	 mat.m[3][2],mat.m[3][3]);

  mat44_to_orientation(mat,&icode,&jcode,&kcode);
  printf("sform_xorient  %s\n",nifti_orientation_string(icode));
  printf("sform_yorient  %s\n",nifti_orientation_string(jcode));
  printf("sform_zorient  %s\n",nifti_orientation_string(kcode));


  printf("file_type      %s\n", FslFileTypeString(fslio->niftiptr->nifti_type));
  printf("file_code      %d\n", fslio->niftiptr->nifti_type);

  /* Data History */
  memset(string,0,128);
  strncpy(string,fslio->niftiptr->descrip,80);
  printf("descrip        %s\n", string);
  memset(string,0,128);
  strncpy(string,fslio->niftiptr->aux_file,24);
  printf("aux_file       %s\n", string);
  /* printf("orient         %d\n", hdr->hist.orient); */

  return;
}

/* }}} */
/* {{{ main */

int main(int argc, char **argv) 
{
  FSLIO* fslio=NULL;
  int argval=1, niftiform=0;

  if (argc<2)
    usage();
    
  if (strcmp(argv[1],"-x")==0) {
      niftiform=1;
      argval=2;
  }
  fslio = FslOpen(FslMakeBaseName(argv[argval]),"rb");
  FslClose(fslio);
  if (niftiform==0) { ShowHdr(argv[argval], fslio); }
  else { ShowNifti(argv[argval],fslio); }

  exit(0);
}

/* }}} */
