
/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(LOOKUPTABLE_H)
#define LOOKUPTABLE_H

#include <boost/shared_ptr.hpp>
#include <string>

#include <vector>
#include "stdio.h"

typedef unsigned char ColorTriplet[3];
typedef unsigned char ColorRGBA[4];

/**
 * @author James Saunders <jim@fmrib.ox.ac.uk>
 * @date   Mon Dec 23 17:59:36 2002
 * 
 * @brief  Implementation of a color look up table (lut).
 * 
 * Provides a color lut suitable for GL rendering along with methods
 * to create standard luts and read custom ones from file.
 */
class LookUpTable
{
public:
  typedef boost::shared_ptr< LookUpTable > Handle;

  static LookUpTable::Handle load(const std::string& filename); 

  static LookUpTable::Handle greyScale();
  static LookUpTable::Handle red();
  static LookUpTable::Handle blue();
  static LookUpTable::Handle green();
  static LookUpTable::Handle yellow();
  static LookUpTable::Handle redYellow();
  static LookUpTable::Handle blueLightblue();
  static LookUpTable::Handle pink();
  static LookUpTable::Handle hot();
  static LookUpTable::Handle cool();
  static LookUpTable::Handle copper();
  static LookUpTable::Handle spectrum();
  static LookUpTable::Handle render1();
  static LookUpTable::Handle render1t();
  static LookUpTable::Handle render2();
  static LookUpTable::Handle render2t();
  static LookUpTable::Handle render3();
  static LookUpTable::Handle cortical();
  static LookUpTable::Handle subcortical();

  void pushValue(unsigned char red, unsigned char green, unsigned char blue, int index);
  const ColorTriplet& inqValue(float f);  
  const ColorTriplet& inqValueIndex(float f);

  void allocateMemory(int size);  
  std::string inqLutName() const;
  bool isIndexLut() const;
  bool isVisible() const;
  bool isAutoSelectable() const;
  
  virtual ~LookUpTable();

private:   
  LookUpTable(const std::string& filename);
  LookUpTable();
  LookUpTable::Handle LoadStdLut(const char ** lutData,
                                 LookUpTable::Handle lut,
                                 int elementCount,
                                 const std::string  name);
  std::string extractName(std::string filename);

  ColorTriplet*  m_lookUpData;
  ColorTriplet   m_blackTriplet;
  int            m_lookUpDataSize;
  std::string    m_fileName;
  std::string    m_lutName;
  bool           m_isIndexLut;
  bool           m_isAutoSelectable;
  bool           m_isVisible;
};


inline std::string LookUpTable::inqLutName() const { return m_lutName; }
inline bool LookUpTable::isIndexLut() const { return m_isIndexLut; }
inline bool LookUpTable::isVisible() const { return m_isVisible; }
inline bool LookUpTable::isAutoSelectable() const { return m_isAutoSelectable; }

#endif
