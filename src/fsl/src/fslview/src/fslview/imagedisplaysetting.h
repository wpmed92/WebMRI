/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(IMAGEDISPLAYSETTING_H)
#define IMAGEDISPLAYSETTING_H

#if defined(WIN32)
#pragma warning(disable:4786)
#endif

#include <boost/shared_ptr.hpp>
#include <memory>
#include "storage/image.h"
#include "lookuptable.h"
#include "bricon.h"
  
typedef enum {None,Lines,RGB} DtiDisplay;

class ImageDisplaySetting  
{
public:
  typedef boost::shared_ptr< ImageDisplaySetting > Handle;
  static Handle create(Image::Handle image, LookUpTable::Handle lut,float trans = 0.5, bool visible = true);
  
  void setTransparency(float trans);
  void setVisibility(bool visible); 
  void setLookUpTable(LookUpTable::Handle);
  void setSecondaryLookUpTable(LookUpTable::Handle);
  float inqTransparency();
  bool inqVisibility();
  void setDtiDisplay(DtiDisplay);
  int  inqDtiDisplay();
  bool inqTransMod();
  void setTransMod(bool);
  void setModTransparency(float);
  float inqModTransparency();
  void setModImage(Image::Handle);
  Handle clone();
  LookUpTable::Handle inqLookUpTable();
  LookUpTable::Handle inqSecondaryLookUpTable();
  BriCon::Handle inqBriCon();
  Image::Handle inqModImage();

  virtual ~ImageDisplaySetting();

private:
  ImageDisplaySetting(Image::Handle image, LookUpTable::Handle lut,
                      float trans, bool visible);
  ImageDisplaySetting(BriCon::Handle bricon, 
		      LookUpTable::Handle lut, LookUpTable::Handle slut,
                      float trans, bool visible,
                      int dti, Image::Handle mod,
                      bool transMod, float modTransVal);

  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;
};



#endif
