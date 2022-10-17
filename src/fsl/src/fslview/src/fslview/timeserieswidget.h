/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(TIMESERIESWIDGET_H)
#define TIMESERIESWIDGET_H

#include <qtoolbar.h>
#include <qstringlist.h>
#include <qcheckbox.h>
#include "viewwidget.h"
#include "singleserieswidget.h"
#include "cursor.h"
#include "qcombobox.h"
#include "qstring.h"

class QSpinBox;
class QToolButton;
class TimeSeriesPlot;
class TimeSeriesToolbar;

class TimeSeriesWidget : public ViewWidget
{
  Q_OBJECT
public:
  TimeSeriesWidget(QWidget *parent,
                   Image::Handle& image,
                   Cursor::Handle& cursor);  

  TimeSeriesWidget(QWidget *parent,
                   Image::Handle& image,
                   Cursor::Handle& cursor,
                   ModelFit::Handle& modelFit);

  virtual ~TimeSeriesWidget();
  
  void addFeatComboBox(QToolBar *);
  
private:

  void constructor();

  Image::Handle    m_image;
  Cursor::Handle   m_cursor;
  PlotOptions::Handle m_options;
  int              m_viewNumber;
  int              m_contrListIndex;
  

  TimeSeriesDisplay::Handle  m_displayWidget;
  TimeSeriesDisplay::Handle  m_singleWidget;
  //  TimeSeriesDisplay::Handle  m_gridWidget;
  //  TimeSeriesDisplay::Handle  m_cubeWidget;
  
  typedef enum {Single, Grid, Cube} GraphView;

  TimeSeriesToolbar*  m_toolbar;
  void                newSlice();

public  slots:
  void singleView(); 
  //  void gridView();
  //  void cubeView();
  void setCurContrItem(int);
  void setCurText(const QString&);
  signals:
  void setEnabled(bool);
  
};

#endif
