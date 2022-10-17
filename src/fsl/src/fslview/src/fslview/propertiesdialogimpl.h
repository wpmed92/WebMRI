#if !defined(PROPERTIESDIALOGIMPL_H)
#define PROPERTIESDIALOGIMPL_H

/* CCOPYRIGHT */

#include "propertiesdialog.h"
#include "properties.h"

class PropertiesDialogImpl : public PropertiesDialog
{
 public:
  static Properties::Handle getProperties(QWidget *);

 private:
  PropertiesDialogImpl(QWidget *);
  ~PropertiesDialogImpl();
};
  
#endif
