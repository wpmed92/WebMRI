/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include <string>
#include <boost/shared_ptr.hpp>

class Preferences
{
public:
  typedef boost::shared_ptr<Preferences> Handle;

  std::string FSLDir() const;
  std::string mni152() const;

  static Handle instance();

private:
  Preferences();

  static Handle m_instance;

  std::string m_fsldir;
  std::string m_mni;
};

