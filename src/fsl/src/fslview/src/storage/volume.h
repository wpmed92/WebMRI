/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

// volume.h: interface for the Volume class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(VOLUME_H)
#define VOLUME_H

#include "fslio/fslio.h"
#include <boost/weak_ptr.hpp>
#include <string>

//! @brief Abstract interface for VolumeStore data.
class Volume
{
public:
  typedef boost::shared_ptr< Volume > Handle;

  virtual short inqX() const = 0;
  virtual short inqY() const = 0;
  virtual short inqZ() const = 0;

  virtual float normalized(short x, short y, short z) = 0;
  virtual float value(unsigned int offset) = 0;
  virtual float value(short x, short y, short z) = 0;
  virtual void  setValue(short x, short y, short z, float value) = 0;

  virtual float inqMin() const = 0;
  virtual float inqMax() const = 0;
  virtual float minValue() const = 0;
  virtual float maxValue() const = 0;
  virtual void calculateMinMax() = 0;
  virtual void calculateRobustMinMax(const bool ignoreZeros = false) = 0;

  virtual bool saveVolume(FSLIO *avw) = 0;

  Volume() {}
  virtual ~Volume() {}
};

//! @brief Templated volume storage class.
template <typename VoxelType>
class VolumeStore: public Volume
{
  friend class TestVolumeStore;

public:
  typedef boost::shared_ptr< VolumeStore<VoxelType> > Handle;

  virtual ~VolumeStore();
	
  VoxelType& operator()(short x, short y, short z);
  VoxelType& operator()(unsigned int offset);

  virtual float normalized(short x, short y, short z);
  virtual float value(short x, short y, short z);
  virtual float value(unsigned int offset) { return m_buffer[offset]; }
  virtual void  setValue(short x, short y, short z, float value);
  virtual short inqX() const;
  virtual short inqY() const;
  virtual short inqZ() const;

  virtual float inqMin() const { return m_min; }
  virtual float inqMax() const { return m_max; }
  virtual float minValue() const;
  virtual float maxValue() const;
  virtual void calculateMinMax();
  virtual void calculateRobustMinMax(const bool ignoreZeros);

  void setMin(float min){m_min = min;}
  void setMax(float max){m_max = max;}

  //  static Handle load(const std::string& filename, size_t n = 0);
  static Handle getVolume(FSLIO* avw, float min, float max, size_t n);
  static Handle clone(const Handle& src);
  static Handle create(short x, short y, short z, VoxelType* buffer);
  static Handle createBlank(short x, short y, short z);

  bool saveVolume(FSLIO *avw);

private:
  VolumeStore(short x, short y, short z, VoxelType* buffer);

  short m_x, m_y, m_z;
  float m_min, m_max;
  bool m_doscaling;
  float m_slope, m_intercept;
  VoxelType* m_buffer;
};

typedef VolumeStore<unsigned char> VolumeB;
typedef VolumeStore<short>         VolumeS;
typedef VolumeStore<int>           VolumeI;
typedef VolumeStore<float>         VolumeF;
typedef VolumeStore<double>        VolumeD;

#include "volume.inc"

#endif
