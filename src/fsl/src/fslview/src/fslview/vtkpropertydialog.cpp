/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "vtkpropertydialog.h"
#include "vtkwidget.h"

#include <qspinbox.h>
#include <qlineedit.h>
#include <qvalidator.h>
#include <qcolordialog.h>
#include <qpushbutton.h>

VTKPropertyDialog::VTKPropertyDialog(QWidget* parent, VTKProperties& props): 
  VTKPropertyDialogBase(parent), m_props(props)
{
  m_lowerThreshold->setValue(m_props.inqLowerThreshold());
  m_upperThreshold->setValue(m_props.inqUpperThreshold());
  m_mcThreshold->setValidator(new QDoubleValidator(this));
  m_mcThreshold->setText(tr("%1").arg(m_props.inqMcThreshold()));
  int x, y, z;
  m_props.inqKernel(x, y, z);
  m_kernelX->setValue(x); m_kernelY->setValue(y); m_kernelZ->setValue(z);

  m_sdX->setValidator(new QDoubleValidator(this));
  m_sdY->setValidator(new QDoubleValidator(this));
  m_sdZ->setValidator(new QDoubleValidator(this));
  float sx, sy, sz;
  m_props.inqStdDev(sx, sy, sz);
  m_sdX->setText(tr("%1").arg(sx)); m_sdY->setText(tr("%1").arg(sy)); m_sdZ->setText(tr("%1").arg(sz));
  m_radiiX->setValidator(new QDoubleValidator(this));
  m_radiiY->setValidator(new QDoubleValidator(this));
  m_radiiZ->setValidator(new QDoubleValidator(this));
  float rx, ry, rz;
  m_props.inqRadii(rx, ry, rz);
  m_radiiX->setText(tr("%1").arg(rx)); m_radiiY->setText(tr("%1").arg(ry)); m_radiiZ->setText(tr("%1").arg(rz));

  m_iterations->setValue(m_props.inqIterations());
  m_relaxFactor->setValidator(new QDoubleValidator(this));
  m_relaxFactor->setText(tr("%1").arg(m_props.inqRelaxationFactor()));

  m_interpMode->setCurrentItem(m_props.inqInterpMode());
  m_ambient->setValidator(new QDoubleValidator(this));
  m_diffuse->setValidator(new QDoubleValidator(this));
  m_opacity->setValidator(new QDoubleValidator(this));
  m_specular->setValidator(new QDoubleValidator(this));
  m_specularPower->setValidator(new QDoubleValidator(this));
  m_ambient->setText(tr("%1").arg(m_props.inqAmbient()));
  m_diffuse->setText(tr("%1").arg(m_props.inqDiffuse()));
  m_opacity->setText(tr("%1").arg(m_props.inqOpacity()));
  m_specular->setText(tr("%1").arg(m_props.inqSpecular()));
  m_specularPower->setText(tr("%1").arg(m_props.inqSpecularPower()));

  m_featureAngle->setValidator(new QDoubleValidator(this));
  m_featureAngle->setText(tr("%1").arg(m_props.inqFeatureAngle()));

  float cr, cg, cb;
  m_props.inqColor(cr, cg, cb);
  m_colorR->setValidator(new QDoubleValidator(this));
  m_colorR->setText(tr("%1").arg(cr));
  m_colorG->setValidator(new QDoubleValidator(this));
  m_colorG->setText(tr("%1").arg(cg));
  m_colorB->setValidator(new QDoubleValidator(this));
  m_colorB->setText(tr("%1").arg(cb));

  //  connect(m_selectButton, SIGNAL(clicked()), this, SLOT(selectColor()));
}

VTKProperties& VTKPropertyDialog::getProperties()
{
  m_props.setLowerThreshold(m_lowerThreshold->value());
  m_props.setUpperThreshold(m_upperThreshold->value());
  m_props.setMcThreshold(m_mcThreshold->text().toDouble());
  m_props.setInterpMode(m_interpMode->currentItem());
  m_props.setKernel(m_kernelX->value(), m_kernelY->value(), m_kernelZ->value());
  m_props.setRelaxationFactor(m_relaxFactor->text().toDouble());
  m_props.setIterations(m_iterations->value());
  m_props.setAmbient(m_ambient->text().toDouble());
  m_props.setDiffuse(m_diffuse->text().toDouble());
  m_props.setOpacity(m_opacity->text().toDouble());
  m_props.setSpecular(m_specular->text().toDouble());
  m_props.setSpecularPower(m_specularPower->text().toDouble());
  m_props.setFeatureAngle(m_featureAngle->text().toDouble());
  m_props.setStdDev(m_sdX->text().toDouble(), m_sdY->text().toDouble(), m_sdZ->text().toDouble());
  m_props.setRadii(m_radiiX->text().toDouble(), m_radiiY->text().toDouble(), m_radiiZ->text().toDouble());
  m_props.setColor(m_colorR->text().toDouble(), m_colorG->text().toDouble(), m_colorB->text().toDouble());
  return m_props;
}

void VTKPropertyDialog::selectColor()
{
  float cr, cg, cb;
  m_props.inqColor(cr, cg, cb);
  QColor orig = QColor(int(cr * 255), int(cg * 255), int(cb * 255));
  QColor color = QColorDialog::getColor(orig, this);
  if (color.isValid()) {
//     colorLabel->setText(color.name());
//     colorLabel->setPalette(QPalette(color));
    m_colorR->setText(tr("%1").arg(color.red()/255.0));
    m_colorG->setText(tr("%1").arg(color.green()/255.0));
    m_colorB->setText(tr("%1").arg(color.blue()/255.0));
  }
}
