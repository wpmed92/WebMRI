#include "propertiesdialogimpl.h"
#include "properties.h"

PropertiesDialogImpl::PropertiesDialogImpl(QWidget *parent)
{
}

PropertiesDialogImpl::~PropertiesDialogImpl()
{
}

Properties::Handle PropertiesDialogImpl::getProperties(QWidget *parent)
{
  Properties::Handle p = Properties::create();
  PropertiesDialogImpl propertiesDialog(parent);

  if(propertiesDialog.exec() == QDialog::Accepted)
    {
    }

  return p;
}
