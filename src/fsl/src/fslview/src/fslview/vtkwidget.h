/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(VTKWIDGET_H)
#define VTKWIDGET_H

#include "imagewidget.h"

#include <vtkProperty.h>

class QVTKWidget;

class VTKProperties
{
public:
  VTKProperties();
  VTKProperties(const VTKProperties& rhs);

  VTKProperties& operator=(const VTKProperties& rhs);

  void Swap(VTKProperties& other);

  int inqLowerThreshold() const { return m_lowerThreshold; }
  void setLowerThreshold(int t) { m_lowerThreshold = t; }

  int inqUpperThreshold() const { return m_upperThreshold; }
  void setUpperThreshold(int t) { m_upperThreshold = t; }

  float inqMcThreshold() const { return m_mcThreshold; }
  void setMcThreshold(float f) { m_mcThreshold = f; }

  int inqInterpMode() const { return m_interpMode; }
  void setInterpMode(int t) { m_interpMode = t; }

  void inqKernel(int& x, int&  y, int&  z) const 
  { x = m_kernelX; y = m_kernelY; z = m_kernelZ; }

  void setKernel(int x, int y, int z) 
  { m_kernelX = x; m_kernelY = y; m_kernelZ = z; }

  void inqStdDev(float& x, float&  y, float&  z) const 
  { x = m_stdDev[0]; y = m_stdDev[1]; z = m_stdDev[2]; }

  void setStdDev(float x, float y, float z) 
  { m_stdDev[0] = x; m_stdDev[1] = y; m_stdDev[2] = z; }

  void inqRadii(float& x, float&  y, float&  z) const 
  { x = m_radii[0]; y = m_radii[1]; z = m_radii[2]; }

  void setRadii(float x, float y, float z) 
  { m_radii[0] = x; m_radii[1] = y; m_radii[2] = z; }

  void inqColor(float& r, float& g, float& b)
  { r = m_colorR; g = m_colorG; b = m_colorB; }

  void setColor(float r, float g, float b)
  { m_colorR = r; m_colorG = g; m_colorB = b; }

  float inqRelaxationFactor() const { return m_relaxationFactor; }
  void setRelaxationFactor(float f) { m_relaxationFactor = f; }

  int inqIterations() const { return m_iterations; }
  void setIterations(int i) { m_iterations = i; }

  float inqAmbient() const { return m_ambient; }
  void setAmbient(float f) { m_ambient = f; }

  float inqDiffuse() const { return m_diffuse; }
  void setDiffuse(float f) { m_diffuse = f; }

  float inqOpacity() const { return m_opacity; }
  void setOpacity(float f) { m_opacity = f; }

  float inqSpecular() const { return m_specular; }
  void setSpecular(float f) { m_specular = f; }

  float inqSpecularPower() const { return m_specularPower; }
  void setSpecularPower(float f) { m_specularPower = f; }

  float inqFeatureAngle() const { return m_featureAngle; }
  void setFeatureAngle(float fa) { m_featureAngle = fa; }

private:
  int m_upperThreshold, m_lowerThreshold;
  float m_mcThreshold;
  int m_interpMode;
  int m_kernelX;
  int m_kernelY;
  int m_kernelZ;
  int m_iterations;
  float m_relaxationFactor;
  float m_ambient, m_diffuse, m_opacity, m_specular, m_specularPower;
  float m_featureAngle;
  float m_stdDev[3];
  float m_radii[3];
  float m_colorR, m_colorG, m_colorB;
};


class VTKWidget : public ImageWidget  
{
  Q_OBJECT
public:
  VTKWidget(QWidget *parent, 
	    ImageGroup::Handle i, 
	    OverlayList::Handle ol,
	    Cursor::Handle c);
  virtual ~VTKWidget();
  virtual void update(const Cursor::Handle);
  virtual void update(const OverlayList*, OverlayListMsg);

private slots:
  void print();
  void options();

private:
  struct Implementation;  
  const std::auto_ptr<Implementation> m_impl;

  ImageGroup::Handle m_image;
  QVTKWidget        *m_vtkwidget;
};

#endif
