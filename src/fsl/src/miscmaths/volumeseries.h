/*  volumeseries.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

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

#include <iostream>
#include <fstream>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "volumeinfo.h"

#include <string>

using namespace NEWMAT;

namespace MISCMATHS {
  
#if !defined(__volumeseries_h)
#define __volumeseries_h

  // Rows are volumes
  // Columns are (time) series
  // num Rows is size of (time) series
  // num Cols is size of volumes
  class VolumeSeries : public Matrix
    {
    public:

      VolumeSeries() : Matrix() {}
      
      VolumeSeries(const VolumeInfo& pvolinfo, const ColumnVector& in) : 
	Matrix(),
	volinfo(pvolinfo),
	preThresholdPositions(in)
	{}

      VolumeSeries(int pnumVols, int pnumSeries) : 
	Matrix(pnumVols, pnumSeries), 
	means(pnumSeries){}
      
      VolumeSeries(int pnumVols, int pnumSeries, const VolumeInfo& pvolinfo, const ColumnVector& in) : 
	Matrix(pnumVols, pnumSeries), 
	volinfo(pvolinfo),
	preThresholdPositions(in),
	means(pnumSeries){}

      VolumeSeries& operator=(const VolumeSeries& vol) { 
	Matrix::operator=(vol); 
	preThresholdPositions = vol.preThresholdPositions;
	volinfo = vol.volinfo;
	means = vol.means;
	return *this;
      }
      VolumeSeries& operator=(const Matrix& mat) {
	Matrix::operator=(mat);
	return *this;
      }

      VolumeSeries& operator=(float pin) { 
	Matrix::operator=(pin); 
	return *this;
      }
      VolumeSeries(const VolumeSeries& vol){operator=(vol);}
      VolumeSeries(const Matrix& pmat) : Matrix(pmat) {}
      
      void thresholdSeries(float thresh, bool removeMean);
      void thresholdSeries();
      void unthresholdSeries();
      void unthresholdSeries(const VolumeInfo& pvolinfo,const ColumnVector& in);
      void removeSeriesMeans();

      const ColumnVector& getPreThresholdPositions() const { return preThresholdPositions; } 
      void setPreThresholdPositions(const ColumnVector& in) { preThresholdPositions = in; }

      int getNumVolumes() const { return Nrows(); }
      int getNumSeries() const { return Ncols(); }
      int nvoxels() const { return Ncols(); }
      int tsize() const { return Nrows(); }

      const ReturnMatrix getSeries(int i) const { ColumnVector tmp = Column(i);tmp.Release();return tmp; }

      ReturnMatrix getSeries(int i) { ColumnVector tmp = Column(i);tmp.Release();return tmp; }

      void setSeries(const ColumnVector& in, int i) { Column(i) = in; }

      const ReturnMatrix getVolume(int i) const { RowVector tmp = Row(i);tmp.Release();return tmp; }
      ReturnMatrix getVolume(int i) { RowVector tmp = Row(i);tmp.Release();return tmp; }

      void setVolume(const ColumnVector& in, int i) { Row(i) = in.AsRow(); }
      void setVolume(const RowVector& in, int i) { Row(i) = in; }

	  void read(const std::string& fname);
	  void writeAsInt(const std::string& fname);
	  void writeAsFloat(const std::string& fname);
	  void writeThresholdedSeriesAsFloat(const VolumeInfo& pvolinfo,const ColumnVector& in,const std::string& fname);
      
      void replaceMeans();

      const VolumeInfo& getInfo() const { return volinfo; }
      void setInfo(const VolumeInfo& pvolinfo) { volinfo = pvolinfo; }

      int getUnthresholdNumSeries() const { return volinfo.x*volinfo.y*volinfo.z; }
      
    protected:
      VolumeInfo volinfo;
      ColumnVector preThresholdPositions; 
      ColumnVector means;
    };

#endif

}




