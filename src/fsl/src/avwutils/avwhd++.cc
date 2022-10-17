//     avwhd++.cc - show image header
//     Steve Smith, Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2005 University of Oxford  
//     COPYRIGHT  

#include "newimage/newimageall.h"
#include <iomanip>
#include <iostream>
using namespace NEWIMAGE;

void print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: avwhd [-x] <input>" << endl;
  cout << "       -x : instead print an XML-style NIFTI header" << endl;
}

void ShowNifti(FSLIO* fslio)
{
    if (fslio == NULL) {
        cerr << "ERROR: Could not open file" << endl;
	return;
    }
    if (fslio->niftiptr!=NULL) {
        cout << nifti_image_to_ascii(fslio->niftiptr) << endl;
	return;
    }
    if (fslio->mincptr!=NULL) {
        cerr << "ERROR: Minc is not currently supported" << endl;
	return;
    }
    return;
}

void ShowHdr(char *fileName, FSLIO* fslio)
{
  int i, ft, isanalyze=0;
  struct dsr *hdr;
  mat44 mat;
  int icode, jcode, kcode;

  if (fslio == NULL) 
  {
    cerr << "ERROR: Could not open file" << endl;;
    return;
  }

  ft = FslGetFileType(fslio);
  if (FslBaseFileType(ft)==FSL_TYPE_MINC) 
  {
    cerr << "ERROR: Minc is not currently supported" << endl;
    return;
  }

  if (fslio->niftiptr == NULL) 
  {
    cerr << "ERROR: Not an Analyze or Nifti file" << endl;
    return;
  }
  //-------------------- ANALYZE CASE ----------------------- 
  if (FslBaseFileType(ft)==FSL_TYPE_ANALYZE) 
  { 
    isanalyze=1; 
    //load raw hdr structure   
    hdr = (struct dsr *)calloc(1,sizeof(struct dsr));
    FslReadRawHeader(hdr,fslio->niftiptr->fname);
    if (fslio->niftiptr->byteorder != short_order()) 
    {
      cout << "Byte swapping" << endl;
      AvwSwapHeader(hdr);
    }
    cout << "filename       " << fileName << endl << endl;
    
    //Header Key
    cout << "sizeof_hdr     " << hdr->hk.sizeof_hdr<< endl;
    cout << "data_type      " << hdr->hk.data_type << endl;
    cout << "db_name        " << hdr->hk.db_name << endl;
    cout << "extents        " << hdr->hk.extents << endl;
    cout << "session_error  " << hdr->hk.session_error << endl;
    cout << "regular        " << hdr->hk.regular << endl;
    cout << "hkey_un0       " << hdr->hk.hkey_un0 << endl;
    //Image Dimension 
    for(i=0;i<8;i++) cout << "dim" << i << "           " << hdr->dime.dim[i] << endl;
    cout << "vox_units      " << hdr->dime.vox_units << endl;
    cout << "cal_units      " << hdr->dime.cal_units << endl;
    cout << "unused1        " << hdr->dime.unused1 << endl;
    cout << "datatype       " << hdr->dime.datatype << endl;
    cout << "bitpix         " << hdr->dime.bitpix << endl;
    cout.setf(ios::fixed);  //need cout.setf(ios::fixed) instead of << fixed in stream for tru64 comp
    for(i=0;i<8;i++) cout << "pixdim" << i << "        " << setprecision(10) << hdr->dime.pixdim[i] << endl;       
    cout.precision(4);                                   
    cout << "vox_offset     " << setw(6) << hdr->dime.vox_offset << endl;
    cout << "funused1       " << setw(6) << hdr->dime.funused1 << endl;
    cout << "funused2       " << setw(6) << hdr->dime.funused2 << endl;
    cout << "funused3       " << setw(6) << hdr->dime.funused3 << endl;
    cout << "cal_max        " << setw(6) << hdr->dime.cal_max << endl;
    cout << "cal_min        " << setw(6) << hdr->dime.cal_min << endl;
    cout << "compressed     " << hdr->dime.compressed << endl;
    cout << "verified       " << hdr->dime.verified << endl;
    cout << "glmax          " << hdr->dime.glmax << endl;
    cout << "glmin          " << hdr->dime.glmin << endl ;
    //Data History 
    cout << "descrip        " <<  hdr->hist.descrip << endl;
    cout << "aux_file       " <<  hdr->hist.aux_file << endl;
    cout << "orient         " << (int)hdr->hist.orient << endl; //need cast else blank
    cout << "originator     " << hdr->hist.originator << endl;
    /*{
      short blah[5];
      memcpy(blah,hdr->hist.originator,5*sizeof(short));
      cout << "origin1        " << blah[0] << endl;
      cout << "origin2        " <<blah[1] << endl;
      cout << "origin3        " << blah[2] << endl;
      }*/
    cout << "origin1        " << (short)hdr->hist.originator[0] << endl; //These lines don't work on
    cout << "origin2        " << (short)hdr->hist.originator[2] << endl; //all platforms... but WHICH 
    cout << "origin3        " << (short)hdr->hist.originator[4] << endl; //ones - that is the question
    cout << "generated      " << hdr->hist.generated << endl;
    cout << "scannum        " << hdr->hist.scannum << endl;
    cout << "patient_id     " <<  hdr->hist.patient_id << endl;
    cout << "exp_date       " << hdr->hist.exp_date << endl;
    cout << "exp_time       " << hdr->hist.exp_time << endl;
    cout << "hist_un0       " << hdr->hist.hist_un0 <<endl;
    cout << "views          " << hdr->hist.views << endl;
    cout << "vols_added     " << hdr->hist.vols_added << endl;
    cout << "start_field    " << hdr->hist.start_field << endl;
    cout << "field_skip     " << hdr->hist.field_skip << endl;
    cout << "omax           " << hdr->hist.omax << endl;
    cout << "omin           " << hdr->hist.omin << endl;
    cout << "smin           " << hdr->hist.smax << endl;
    cout << "smin           " << hdr->hist.smin << endl;
    cout << "file_type      " << FslFileTypeString(0) << endl;
    cout << "file_code      0" << endl;
    return;
  }
    /* -------------------- NIFTI CASE ----------------------- */
  if (fslio->niftiptr->byteorder != short_order()) 
  { 
    cout << "Byte swapping" << endl;
  }

  cout << "filename       " <<  fslio->niftiptr->fname << endl << endl;
  cout << "sizeof_hdr     " <<  "348" << endl;
  cout << "data_type      " <<  nifti_datatype_string(fslio->niftiptr->datatype) << endl;
  for(i=0;i<8;i++)  cout << "dim" << i << "           " <<  fslio->niftiptr->dim[i] << endl;
  cout << "vox_units      " <<  nifti_units_string(fslio->niftiptr->xyz_units) << endl;
  cout << "time_units     " <<  nifti_units_string(fslio->niftiptr->time_units) << endl;
  cout << "datatype       " <<  fslio->niftiptr->datatype << endl;
  cout << "nbyper         " <<  fslio->niftiptr->nbyper << endl;
  cout << "bitpix         " <<  fslio->niftiptr->nbyper * 8 << endl;
  cout.setf(ios::fixed);  //need cout.setf(ios::fixed) instead of << fixed in stream for tru64 comp
  for(i=0;i<8;i++) cout << "pixdim" << i << "        " << setprecision(10) <<  fslio->niftiptr->pixdim[i] << endl;
  cout << "vox_offset     " <<   fslio->niftiptr->iname_offset << endl;
  cout.precision(4);
  cout << "cal_max        " << setw(6) <<  fslio->niftiptr->cal_max << endl;
  cout << "cal_min        " << setw(6) <<  fslio->niftiptr->cal_min << endl;
  cout.precision(6);
  cout << "scl_slope      " <<  fslio->niftiptr->scl_slope << endl;
  cout << "scl_inter      " <<  fslio->niftiptr->scl_inter << endl;
  cout << "phase_dim      " <<  fslio->niftiptr->phase_dim << endl;
  cout << "freq_dim       " <<  fslio->niftiptr->freq_dim  << endl;
  cout << "slice_dim      " <<  fslio->niftiptr->slice_dim << endl;
  cout << "slice_name     " <<  nifti_slice_string(fslio->niftiptr->slice_code) << endl;
  cout << "slice_code     " <<  fslio->niftiptr->slice_code << endl;
  cout << "slice_start    " <<  fslio->niftiptr->slice_start << endl;
  cout << "slice_end      " <<  fslio->niftiptr->slice_end << endl;
  cout << "slice_duration " <<  fslio->niftiptr->slice_duration << endl;
  cout << "time_offset    " <<  fslio->niftiptr->toffset << endl;
  cout << "intent         " <<  nifti_intent_string(fslio->niftiptr->intent_code) << endl;
  cout << "intent_code    " <<  fslio->niftiptr->intent_code << endl;
  cout << "intent_name    " <<  fslio->niftiptr->intent_name << endl;
  cout << "intent_p1      " <<  fslio->niftiptr->intent_p1 << endl;
  cout << "intent_p2      " <<  fslio->niftiptr->intent_p2 << endl;
  cout << "intent_p3      " <<  fslio->niftiptr->intent_p3 << endl;
  cout << "qform_name     " <<  nifti_xform_string(fslio->niftiptr->qform_code) << endl;
  cout << "qform_code     " <<  fslio->niftiptr->qform_code << endl;
  mat = fslio->niftiptr->qto_xyz;
  for(i=1;i<=4;i++) cout << "qto_xyz:" << i << "      " << mat.m[i-1][0] << "  " << mat.m[i-1][1] << "  " << mat.m[i-1][2] << "  " << mat.m[i-1][3] << endl;
  mat44_to_orientation(mat,&icode,&jcode,&kcode);
  cout << "qform_xorient  " << nifti_orientation_string(icode) << endl;
  cout << "qform_yorient  " << nifti_orientation_string(jcode) << endl;
  cout << "qform_zorient  " << nifti_orientation_string(kcode) << endl;
  cout << "sform_name     " << nifti_xform_string(fslio->niftiptr->sform_code) << endl;
  cout << "sform_code     " << fslio->niftiptr->sform_code << endl;
  mat = fslio->niftiptr->sto_xyz;
  for(i=1;i<=4;i++) cout << "sto_xyz:" << i << "      " << mat.m[i-1][0] << "  " << mat.m[i-1][1] << "  " << mat.m[i-1][2] << "  " << mat.m[i-1][3] << endl;
  mat44_to_orientation(mat,&icode,&jcode,&kcode);
  cout << "sform_xorient  " << nifti_orientation_string(icode) << endl;
  cout << "sform_yorient  " << nifti_orientation_string(jcode) << endl;
  cout << "sform_zorient  " << nifti_orientation_string(kcode) << endl;
  cout << "file_type      " << FslFileTypeString(fslio->niftiptr->nifti_type) << endl;
  cout << "file_code      " << fslio->niftiptr->nifti_type << endl;
  //Data History 
  cout << "descrip        " << fslio->niftiptr->descrip << endl;
  cout << "aux_file       " << fslio->niftiptr->aux_file << endl;
  /*cout << "orient         %d\n", hdr->hist.orient); */
  return;
}


int main(int argc,char *argv[])
{
  if (argc < 2) 
  {
    print_usage(string(argv[0]));
    return 1; 
  }
  FSLIO* fslio=NULL;
  int argval=1, niftiform=0;
  if (strcmp(argv[1],"-x")==0) 
  {
      niftiform=1;
      argval=2;
  }
  fslio = FslOpen(FslMakeBaseName(argv[argval]),"rb");
  FslClose(fslio);
  if (niftiform==0) ShowHdr(argv[argval], fslio); 
  else ShowNifti(fslio); 
  return 0;
}


