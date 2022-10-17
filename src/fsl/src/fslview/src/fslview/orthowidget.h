/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(ORTHOWIDGET_H)
#define ORTHOWIDGET_H

#include "imagewidget.h"
#include "slicewidget.h"
#include "sliceview.h"
#include "overlaylist.h"

class OverlayWidget;
class QGridLayout;
class QTimer;
class QWidget;

class OrthoWidget : public ImageWidget
{
  Q_OBJECT
public:
  OrthoWidget(QWidget *parent, ImageGroup::Handle& i, 
              OverlayList::Handle ol, Cursor::Handle& c );
  virtual ~OrthoWidget();
 
  void update(const OverlayList*, OverlayListMsg);

signals:
  void  volChanged(int n);

private:
  void setLabels(const OverlayList*);

  QTimer            *m_timer;
  //SliceListHandle    m_slices;
  SliceView         *m_coronal;
  SliceView         *m_sagittal;
  SliceView         *m_axial;
  QGridLayout       *m_grid; 
  ImageGroup::Handle m_image;
  QWidget           *m_centralWidget;
};

#endif
