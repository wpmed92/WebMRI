/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "timeserieswidget.h"
#include "cursor.h"
#include "tracker.h"
#include "modelfit.h"
#include "storage/timeseries.h"
#include <qlayout.h>
#include <qpushbutton.h>
#include <qspinbox.h>
#include <qbuttongroup.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qtooltip.h>
#include <list>
#include <qlabel.h>
#include <qcheckbox.h>

#include "singleserieswidget.h"

#include "timeseriestoolbar.h"

#include <qwt_plot.h>

class Singleserieswidget;


TimeSeriesWidget::TimeSeriesWidget(QWidget *parent, 
                                   Image::Handle& image, 
                                   Cursor::Handle& cursor):
  ViewWidget(parent),m_image(image),
  m_cursor(cursor),m_viewNumber(0)
{  
  m_options = PlotOptions::create();

  constructor();
}

TimeSeriesWidget::TimeSeriesWidget(QWidget *parent, 
                                   Image::Handle& image, 
                                   Cursor::Handle& cursor,
                                   ModelFit::Handle& modelFit):
  ViewWidget(parent),m_image(image),
  m_cursor(cursor),m_viewNumber(0)
{  
  m_options = PlotOptions::create();
  m_options->setModelFit(modelFit); 
  m_options->setFeedBack(true);
  constructor();
}

void TimeSeriesWidget::constructor()
{  
  TRACKER("TimeSeriesWidget::constructor");  

  QToolBar *tb = new QToolBar(this);
  m_toolbar = new TimeSeriesToolbar(tb);
  m_contrListIndex=0;
  
  m_displayWidget = TimeSeriesDisplay::Handle(new SingleSeriesWidget(this,m_image,m_cursor,m_options));

  connect(m_toolbar->m_addButton,SIGNAL(pressed()),
          m_displayWidget.get(),SLOT(addTimeSeries()));
  connect(m_toolbar->m_removeButton,SIGNAL(pressed()),
          m_displayWidget.get(),SLOT(remTimeSeries()));
  connect(m_toolbar->m_demeanButton,SIGNAL(toggled(bool)),
          m_displayWidget.get(),SLOT(demeanButtonToggle(bool)));
  connect(m_toolbar->m_percentButton,SIGNAL(toggled(bool)),
	  m_displayWidget.get(),SLOT(percentButtonToggle(bool)));
  connect(this,SIGNAL(setEnabled(bool)),
          m_displayWidget.get(),SLOT(setEnabled(bool)));
  connect(m_toolbar->m_axisButton,SIGNAL(pressed()),
          m_displayWidget.get(),SLOT(axisDisplay()));
  connect(m_toolbar->m_contrastComboBox, SIGNAL(activated(int)), SLOT(setCurContrItem(int)));
  if(m_options->inqFeatMode())
    {
      m_toolbar->m_addButton->setEnabled(false);
      m_toolbar->m_removeButton->setEnabled(false);
      m_toolbar->m_contrastComboBox->setEnabled(true);
      m_toolbar->m_contrastComboBox->clear();
      m_toolbar->m_contrastComboBox->insertItem("No model");
      m_toolbar->m_contrastComboBox->insertItem("Full model only");
      for(unsigned int i=0; i<m_options->getModelFit()->numFits(); i++)
	m_toolbar->m_contrastComboBox->insertItem(m_options->getModelFit()->getConName(i));
    }

  addToolBar(tb,Top,FALSE);
  setCentralWidget(m_displayWidget.get());  
  m_displayWidget->show();
  m_cursor->setCursor(m_cursor->inqX(),m_cursor->inqY(),m_cursor->inqZ());
}

TimeSeriesWidget::~TimeSeriesWidget()
{
}

void TimeSeriesWidget::singleView()
{
    m_viewNumber = TimeSeriesWidget::Single;
    newSlice();
}

void TimeSeriesWidget::newSlice()
{
  TRACKER("TimeSeriesWidget::newSlice");

  setCentralWidget(NULL);

  if(m_displayWidget)
    {
    m_displayWidget->hide();  
    m_displayWidget->setEnabled(false);
     
    disconnect(m_toolbar->m_demeanButton,SIGNAL(toggled(bool)),
               m_displayWidget.get(),SLOT(demeanButtonToggle(bool)));
    disconnect(m_toolbar->m_addButton,SIGNAL(pressed()),
               m_displayWidget.get(),SLOT(addTimeSeries()));
    disconnect(m_toolbar->m_removeButton,SIGNAL(pressed()),
               m_displayWidget.get(),SLOT(remTimeSeries()));
    disconnect(this,SIGNAL(setEnabled(bool)),
               m_displayWidget.get(),SLOT(setEnabled(bool)));
    disconnect(m_toolbar->m_axisButton,SIGNAL(pressed()),
               m_displayWidget.get(),SLOT(axisDisplay()));
    }

  switch(m_viewNumber)
    {
    case TimeSeriesWidget::Single:

      if(!m_singleWidget)
        {
          m_singleWidget =  TimeSeriesDisplay::Handle
	    (new SingleSeriesWidget(this,m_image,m_cursor,m_options));
        }
      m_displayWidget = m_singleWidget;
      m_toolbar->m_addButton->setEnabled(m_options->inqAddRemEnabled());
      m_toolbar->m_removeButton->setEnabled(m_options->inqAddRemEnabled());      
      m_toolbar->m_demeanButton->setEnabled(true);
      break;
    
//     case TimeSeriesWidget::Grid:
//       if(!m_gridWidget)
//         {       
//           m_gridWidget =  TimeSeriesDisplay::Handle
// 	    (new GridSeriesWidget(this,m_image,m_cursor,m_options,0)); 
//         } 
//       m_displayWidget = m_gridWidget;
//       m_toolbar->m_addButton->setEnabled(false);      
//       m_toolbar->m_removeButton->setEnabled(false);
//       m_toolbar->m_demeanButton->setEnabled(false);
//       break;
  
//     case TimeSeriesWidget::Cube:
//       if(!m_cubeWidget)
//         {      
//           m_cubeWidget = TimeSeriesDisplay::Handle
// 	    (new CubeSeriesWidget(this,m_image,m_cursor,m_options));
//         }
//       m_displayWidget = m_cubeWidget;
//       m_toolbar->m_addButton->setEnabled(false);      
//       m_toolbar->m_removeButton->setEnabled(false);
//       m_toolbar->m_demeanButton->setEnabled(false);      
//       break;
     }
  
  connect(m_toolbar->m_demeanButton,SIGNAL(toggled(bool)),
          m_displayWidget.get(),SLOT(demeanButtonToggle(bool)));
  connect(m_toolbar->m_addButton,SIGNAL(pressed()),
          m_displayWidget.get(),SLOT(addTimeSeries()));
  connect(m_toolbar->m_removeButton,SIGNAL(pressed()),
          m_displayWidget.get(),SLOT(remTimeSeries()));
  connect(this,SIGNAL(setEnabled(bool)),
          m_displayWidget.get(),SLOT(setEnabled(bool)));
  connect(m_toolbar->m_axisButton,SIGNAL(pressed()),
          m_displayWidget.get(),SLOT(axisDisplay()));

  m_displayWidget->setEnabled(true);

  setCentralWidget(m_displayWidget.get());  
  m_displayWidget->show();
  m_cursor->setCursor(m_cursor->inqX(),m_cursor->inqY(),m_cursor->inqZ());
}

void TimeSeriesWidget::setCurContrItem(int curItem)
{
  TRACKER("TimeSeriesWidget::setCurContrItem(int curItem)");
  
  MESSAGE(QString("curItem = %1").arg(curItem));

  ModelFit::Handle m(m_options->getModelFit());
  m_contrListIndex=curItem-2;
  switch(curItem) {
  case 0:
    m->showFull(false);
    m->showPartial(false);
    break;
  case 1:
    m->showFull(true);
    m->showPartial(false);
    break;
  default:
    m->showFull(true);
    m->showPartial(true);
    m->curFit(m_contrListIndex);
    break;
  }
}

void TimeSeriesWidget::setCurText(const QString &str)
{
//m_contrList[m_contrListIndex]=str;
}

