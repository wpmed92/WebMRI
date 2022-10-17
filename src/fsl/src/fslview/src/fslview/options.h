/*  FSLView - 2D/3D Interactive Image Viewer

    David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(FSLVIEWOPTIONS_H)
#define FSLVIEWOPTIONS_H

#include "utils/options.h"
#include "fslio/fslio.h"

#include <qfileinfo.h>

namespace Utilities {

  extern bool string_to_T(std::pair<float,float>&, const string&);

}

class OverlayOption
{
public:
  OverlayOption(const string& filename, 
		const Utilities::Option<string>& lutname, 
		const Utilities::Option< std::pair<float,float> >& bricon):
    m_filename(filename), m_lutname(lutname), m_bricon(bricon)
  {
  }
  
  QFileInfo fileInfo() const
  {
    return QFileInfo(FslMakeBaseName(m_filename.c_str()));
  }

  bool lutSpecified() const { return m_lutname.set(); }
  const string& lutname() const { return m_lutname.value(); }

  bool briconSpecified() const { return m_bricon.set(); }
  float min() const { return m_bricon.value().first; }
  float max() const { return m_bricon.value().second; }
  
  virtual ~OverlayOption() {}

private:
  //  OverlayOption() {}

  string m_filename;
  Utilities::Option<string> m_lutname;
  Utilities::Option< std::pair<float,float> > m_bricon;
};

typedef std::list<OverlayOption> OverlayOptionList;

#endif
