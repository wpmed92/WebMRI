/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney, Stephen Smith and V Rama Aravind, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

/****************************************************************************
** $Id: application.h,v 1.51 2006/03/17 19:05:22 flitney Exp $
**
** Copyright (C) 2002 University of Oxford.  All rights reserved.
**
*****************************************************************************/

#if !defined (APPLICATION_H)
#define APPLICATION_H

#if defined(WIN32)
#pragma warning(disable:4786)
#endif

#include <qmainwindow.h>
#include "cursor.h"
#include <qlist.h>
#include "imagegroup.h" 
#include "overlaylist.h"
#include "properties.h"
#include "options.h"
#include "command.h"
#include "modelfit.h"

//#include "boost/smart_ptr.hpp"
//#include "featmodel.h"

class QLabel;
class QString;
class QStatusBar;
class QScrollView;
class QMultiLineEdit;
class QToolBar;
class QPopupMenu;
class QWorkspace;
class FileError;
class ViewWidget;
class FileOpen;

#include "applicationwindowbase.h"

class ApplicationWindow: public ApplicationWindowBase, 
			 public CursorObserver
{
  Q_OBJECT
public:
  ApplicationWindow(const OverlayOptionList&);
  virtual ~ApplicationWindow();
  virtual void update(const Cursor::Handle&);

  bool loadFile(const QString & absFilePath);
  void setCurrentDir(QString);
  void loadFeat(QString &);

  //public slots:
  //  void saveOverlay();
  //  void addOverlay();  
  //  void remOverlay();
  //  void fileCreateMask();
  //  void fileCloseWindow();
  
private slots:
  void fileOpen();
  void fileOpen152();
  void fileSaveAs();
  void fileAdd();
  void fileRemove();
  void fileCreateMask();
  void fileCloseWindow();

//   void windowCascade();
//   void windowTile();

  void viewOrthographic();
  void viewLightbox();
  void viewSingle();
  void view3d();
  void viewTimeseries();  
  void viewImageHistogram();
  void viewClusterBrowser();

  void applicationProperties();
  void helpAbout();
  void helpAboutQt();

  void fileMenuAboutToShow();  
  void viewMenuAboutToShow();
  void windowMenuAboutToShow();
  void windowMenuActivated( int id );
  void displayMessage(const QString &, int);
  void setMenuItems_NoImages(void);

public slots:
  void addLookUpTable();
  void childWindowClose(QCloseEvent*);
  
signals:
  void message(const QString&, int );
  // emits a signal to inform that the work space is empty; used to set/reset menu options
  void workSpaceEmpty(void);

private:
  QWorkspace* m_ws;
//   QMenuBar* m_mb;
//   QPopupMenu* m_fileMenu;
//   QPopupMenu* m_windowsMenu;
//   QPopupMenu* m_toolsMenu;
  
//   int m_fileOpenID,m_fileCloseID,m_fileAddOverlayID,m_fileRemOverlayID,
//     m_fileCreateMaskID, m_fileSaveOverlayID,m_fileOpenFeatID;
//   int m_imageAddLutID;
//   int m_viewOrthoID,m_viewLightID,m_viewSingleID,m_viewVtkID,
//       m_viewHistID,m_viewTimeID,m_viewTimeGridID;
  bool loadOverlay(const QString & absFilePath);
  void buildMenus();
  bool checkAbsFilePath(const QString & absFilePath, QString ext);
  void setupStatusBar();
  void constructToolBar();
  void connectControls();
  bool windowListEmpty();
  int  windowListCount();
  void catchFileError(FileError);
  void viewShow(ViewWidget* v);
  bool tarnishCheck(Image::Handle &);
  bool tarnishCheck();
  void removeExtensions(QString & fileName);
  bool checkFilesExist(const QString & baseName, bool justImg); 
  bool checkSpecificFilesExist(const QString & baseName);
  bool checkForDuplicates(const QString & fn);
  OverlayList::Handle activeOverlayList();
  OverlayList::Handle copyActiveOverlayList();

  //helper functions for initializing PopUpMenus items; Rama - 3/11/04
  void initFileMenuItems(bool);
  void initViewMenuItems(bool);
  void setFileMenuItemsState(void);
  void setViewMenuItemsState(void);

  QString            m_filename;
  QStatusBar        *m_statusBar;

  QPopupMenu        *m_viewMenu;
  QScrollView       *m_sv;

  ImageGroup::Handle m_imageGroup;
  Cursor::Handle     m_cursor;

  Properties::Handle m_properties;
  
  ModelFit::Handle m_modelFit;
  
  FileOpen* m_fileOpen;

protected:
  virtual void closeEvent(QCloseEvent*);
};


#endif
