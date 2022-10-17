/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include "preferences.h"

using namespace std;

Preferences::Handle Preferences::m_instance;

Preferences::Preferences(): 
  m_fsldir("/usr/local/fsl"),
  m_mni("/usr/local/fsl/etc/standard/avg152T_brain.hdr") 
{}

//! @brief Access Preferences singleton
//!
//! @return a handle to be used as a pointer to a Preferences object
Preferences::Handle Preferences::instance()
{
  if(!m_instance)
    m_instance = Preferences::Handle(new Preferences());

  return m_instance;
}
   
//! @brief Prefered value of FSLDIR
//!
//! Returns the prefered value of FSLDIR
//!
//! @return The users prefered value of FSLDIR
string Preferences::FSLDir() const 
{
  return m_fsldir; 
}

//! @brief Prefered location of mni152T1_brain
//!
//! Return the prefered location where we can find mni152T1_brain image
//!
//! @return The users prefered location for the mni152T1_brain image
string Preferences::mni152() const
{
  return m_mni;
}

