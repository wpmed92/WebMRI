/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "imagewidget.h"
#include "slicewidget.h"
#include <qwidget.h>
#include <qcheckbox.h>
#include <qtoolbar.h>
#include <qdockwindow.h>
#include <qspinbox.h>
#include <qslider.h>
#include <qtooltip.h>
#include <qlayout.h>
#include <qpixmap.h>
#include <qtoolbutton.h>
#include <qbuttongroup.h>
#include <qmessagebox.h>
#include <qstatusbar.h>

#include "maintoolbar.h"
#include "modetoolbar.h"
#include "overlaylist.h"
#include "overlaywidget.h"
#include "overlayinfodialog.h"
#include "briconwidget.h"
#include "drawwidget.h"
#include "cursorwidget.h"

#include "tracker.h"

ImageWidget::ImageWidget(QWidget* parent,        ImageGroup::Handle i,
                         OverlayList::Handle ol, Cursor::Handle c): 
  ViewWidget(parent), m_imageGroup(i), m_cursor(c), m_overlayList(ol)
{
  setFocusPolicy(ClickFocus);

  font().setPointSize(8);

  if(!m_overlayList)
      m_overlayList = OverlayList::create(m_imageGroup);
  m_drawSettings = DrawSettings::create();
  
  m_overlayDialog = NULL;
 
  constructToolBar();

  m_overlayList->attach(this);
  m_cursor->attach(this);

  connect(parent, SIGNAL(windowActivated(QWidget*)), this, SLOT(windowActivated(QWidget*)));

  m_mainToolbarWidget->setCursorMode();
}

ImageWidget::~ImageWidget()
{

  TRACKER("ImageWidget::~ImageWidget()");
  m_overlayWidget->close();
  m_overlayList->detach(this); 
  m_cursor->detach(this);
  m_movieTimer->stop();
}

void ImageWidget::windowActivated(QWidget *w) 
{
  TRACKER("ImageWidget::windowActivated(QWidget *w)");
  if(w == this) {
    setFocusPolicy(QWidget::StrongFocus);
    MESSAGE("setFocusPolicy(QWidget::StrongFocus)");
  } else {
    setFocusPolicy(QWidget::NoFocus);
    MESSAGE("setFocusPolicy(QWidget::NoFocus)");
  }
}

QSize ImageWidget::sizeHint() const
{
  return QSize(600, 600);
}

void ImageWidget::constructToolBar()
{
  m_toolbar = new QToolBar(this, "Display mode tools");
  QToolBar *modebar = new QToolBar(this, "Advanced mode tools");
  QDockWindow *cursordock = new QDockWindow(QDockWindow::InDock, this);
  cursordock->setCloseMode(QDockWindow::Always);
  cursordock->setResizeEnabled(true);
  cursordock->setCaption( tr("Cursor tools") );
  setAppropriate(cursordock, true);

  m_briconToolbar = new QToolBar(this, "Brightness/contrast controls");

  m_overlayDock = new QDockWindow(QDockWindow::InDock, this);
  m_overlayDock->setCloseMode(QDockWindow::Always);
  m_overlayDock->setResizeEnabled(true);
  m_overlayDock->setCaption( tr("Overlay settings") );
  setAppropriate(m_overlayDock, true);

  m_drawToolbar = new QToolBar(this, "Draw toolbar");

  ImageInfo::Handle info(m_imageGroup->getMainImage()->getInfo());

  m_modeWidget     = new ModeToolBarWidget(modebar);
  m_briconWidget   = new BriConWidget(m_briconToolbar,m_overlayList);
  m_overlayWidget  = new OverlayWidget(m_overlayDock, m_overlayList); 
  m_drawWidget     = new DrawWidget(m_drawToolbar, m_overlayList, m_drawSettings); m_drawWidget->setEnabled(false);
  m_cursorWidget   = new CursorWidget(cursordock, m_cursor, m_overlayList);

  m_overlayDock->setWidget(m_overlayWidget);
  cursordock->setWidget(m_cursorWidget);

  connect(m_overlayWidget,SIGNAL(infoButtonAction()),
          this,           SLOT(openOverlayDialog()));

  m_mainToolbarWidget = new MainToolBarWidget(m_toolbar, 0, info->inqNumVolumes() - 1);
  m_mainToolbarWidget->enableMaskMode(false);

  connect(m_mainToolbarWidget,  SIGNAL(zoomValueChanged(int)),      SIGNAL(zoomValueChanged(int)));
  connect(m_mainToolbarWidget,  SIGNAL(crossHairStateChanged(int)), SLOT(crossHairModeChanged(int)));
  connect(m_mainToolbarWidget,  SIGNAL(modeChanged(SliceWidget::Mode)), 
	  SLOT(changeMode(SliceWidget::Mode)));
  connect(m_mainToolbarWidget, SIGNAL(resetZoomClicked()), SIGNAL(resetZoom()));
  connect(m_modeWidget, SIGNAL(movieStateChanged(int)), SLOT(toggleMovie(int)));
  connect(m_drawWidget, SIGNAL(undoButtonClicked()), SLOT(undoGraphics()));

  connect(m_cursorWidget,  SIGNAL(volumeValueChanged(int)),    SLOT(setVolumeValue(int)));

  addDockWindow(m_toolbar,       tr("Main mode tools"), Top, FALSE);
  addDockWindow(modebar,         tr("Animation controls"), Top, FALSE);
  addDockWindow(m_briconToolbar, tr("Brightness/contrast tools"), Top, FALSE);  
  addDockWindow(m_drawToolbar,   tr("Pen/drawing palette"), Top, FALSE);
  addDockWindow(cursordock,      tr("Cursor tools"), Bottom, FALSE);
  addDockWindow(m_overlayDock,   tr("Overlay settings"), Bottom, FALSE);

  m_drawToolbar->hide();

  m_movieTimer = new QTimer(this);
  connect( m_movieTimer, SIGNAL(timeout()), SLOT(nextFrame()) );
   
  m_modeWidget->enableMovieMode(info->inqNumVolumes() > 1);
}

void ImageWidget::changeMode(SliceWidget::Mode m)
{
  switch(m) {
  case SliceWidget::Masking:
    m_drawToolbar->show();
    m_drawWidget->setEnabled(true);
    m_drawWidget->updateControls();
    break;
  default:
    m_drawToolbar->hide();
    m_drawWidget->setEnabled(false);
    break;
  }

  emit modeChanged(m);
}

void ImageWidget::crossHairModeChanged(int state) 
{ 
  emit crossHairModeChanged(state == QButton::On); 
}


void ImageWidget::setVolumeValue(int n)
{
  TRACKER("ImageWidget::setVolume");
  m_cursor->setCursor(m_cursor->inqX(),m_cursor->inqY(),m_cursor->inqZ(),n);
}

void ImageWidget::setZoomValue(int f)
{
   TRACKER("ImageWidget::setZoomValue");
   m_mainToolbarWidget->setZoomValue(f);
}

void ImageWidget::openOverlayDialog()
{
  TRACKER("ImageWidget::openOverlayDialog");
  if(m_overlayDialog)delete m_overlayDialog;
  m_overlayDialog = new OverlayInfoDialog(this,m_overlayList,m_imageGroup);
  m_overlayDialog->show();
  connect(m_overlayDialog, SIGNAL(message(const QString&, int)), SIGNAL(message(const QString&, int)));
}

void ImageWidget::update(const Cursor::Handle& c)
{
  TRACKER("ImageWidget::update(const Cursor::Handle& c)");

  m_cursorWidget->setVolumeValue(c->inqV());
}
  
void ImageWidget::update(const OverlayList* i, OverlayListMsg msg)
{
  TRACKER("ImageWidget::update(const OverlayList* i, OverlayListMsg msg)");

  if(OverlayListMsg(Select) == msg){clearUndoList();}
  if(OverlayListMsg(DtiMode) == msg)
    {
      MetaImage::Handle m = i->getActiveMetaImage();
      if(m.get())
      {
        dtiDisplayMode(m->getDs()->inqDtiDisplay());
      }
    }
  if(OverlayListMsg(Select) == msg || OverlayListMsg(DtiMode) == msg || 
     OverlayListMsg(Visibility) == msg || OverlayListMsg(Security) == msg)
    {
   
    bool state(i->getActiveMetaImage()->inqVisibility() && !i->getActiveMetaImage()->inqReadOnly());
    m_mainToolbarWidget->enableMaskMode(state);
    if(!state)
      {
	m_drawWidget->setEnabled(state);
      }
    else if(m_mainToolbarWidget->inqMaskMode())
      {
        m_drawWidget->setEnabled(true);
	m_drawSettings->setMode(DrawSettings::FreeHand);
      }
    else 
      {
	m_drawWidget->setEnabled(false);
	m_mainToolbarWidget->setCursorMode();
      }
    }
  MESSAGE("finish");
}  


void ImageWidget::clearUndoList()
{
  m_undoList.clear();
}


void ImageWidget::undoGraphics()
{
  if(!m_undoList.empty())
    {
      Shape::Handle undoBuffer = m_undoList.back();
      m_undoList.pop_back();
      undoBuffer->commit();
    } 
  m_cursor->repaint();
}

OverlayList::Handle ImageWidget::getOverlayList()
{
 return  m_overlayList;
}

void ImageWidget::nextFrame()
{
  ImageInfo::Handle info(m_imageGroup->getMainImage()->getInfo());
  int n = m_cursor->inqV();
  if(++n >= info->inqNumVolumes())n = 0;
 
  m_cursor->setCursor(m_cursor->inqX(), m_cursor->inqY(), m_cursor->inqZ(), n);

  setVolumeValue(n);  
}

void ImageWidget::toggleMovie(int state)
{
  if(!m_movieTimer->isActive())
    m_movieTimer->start(50, false);
  else
    m_movieTimer->stop();
}  

void ImageWidget::dtiDisplayMode(int dtiMode)
{
  if(DtiDisplay(None) == dtiMode)
    {  
      ImageInfo::Handle info(m_imageGroup->getMainImage()->getInfo());
      m_modeWidget->enableMovieMode(info->inqNumVolumes() > 1);
      m_cursorWidget->enableVolumeSpinBox(true); 
    }
  else
    { 
      m_cursorWidget->setVolumeValue(0);
      m_cursorWidget->enableVolumeSpinBox(false); 
      m_modeWidget->enableMovieMode(false); 
      m_movieTimer->stop();
    }
}
