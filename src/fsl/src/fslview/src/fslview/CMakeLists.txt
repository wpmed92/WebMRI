SET(srcs
  main.cpp
  logger.cpp
  preferences.cpp
  application.cpp
  splashscreen.cpp 
  version.cpp
  tracker.cpp
  cursor.cpp
  bricon.cpp
  imagegroup.cpp 
  overlaylist.cpp 
  metaimage.cpp 
  imagedata.cpp 
  drawsettings.cpp 
  imagedisplaysetting.cpp 
  lookuptable.cpp 
  curvedatalist.cpp 
  imagedatastore.cpp 
  imagebuffer.cpp 
  filemanager.cpp 
  graphmanager.cpp 
  briconwidget.cpp 
  cursorwidget.cpp 
  drawwidget.cpp 
  imagewidget.cpp 
  viewwidget.cpp 
  orthowidget.cpp 
  singlewidget.cpp 
  lightboxwidget.cpp 
  slicewidget.cpp 
  sliceview.cpp 
  timeserieswidget.cpp 
  singleserieswidget.cpp 
  overlaywidget.cpp 
  histogramwidget.cpp 
  clusterbrowser.cpp
  histogramtoolbar.cpp 
  overlayinfodialog.cpp
  histogramoptionsdialogimpl.cpp 
  properties.cpp 
  propertiesdialogimpl.cpp 
  vtkpropertydialog.cpp 
  vtktoolbar.cpp 
  timeseriestoolbar.cpp
  shape.cpp 
  rect.cpp 
  command.cpp 
  modelfit.cpp
  clusterdata.cpp
  tsplotcode.cpp 
  maintoolbar.cpp 
  modetoolbar.cpp 
  vtkwidget.cpp)

QT_WRAP_CPP( ui moc_srcs
  application.h
  overlayinfodialog.h
  imagewidget.h
  viewwidget.h
  sliceview.h
  orthowidget.h
  singlewidget.h
  lightboxwidget.h
  briconwidget.h
  cursorwidget.h
  drawwidget.h
  slicewidget.h
  splashscreen.h
  timeserieswidget.h
  clusterbrowser.h
  singleserieswidget.h
  maintoolbar.h
  modetoolbar.h
  histogramwidget.h
  vtkwidget.h
  overlaywidget.h)

QT_WRAP_UI(ui ui_hdrs ui_srcs
  sliceviewbase.ui
  applicationwindowbase.ui
  clusterbrowserbase.ui
  histogramoptionsdialog.ui 
  histogramtoolbarbase.ui 
  propertiesdialog.ui 
  overlaywidgetbase.ui 
  createmaskdialog.ui 
  overlayinfodialogbase.ui 
  timeseriestoolbarbase.ui 
  maintoolbarbase.ui 
  modetoolbarbase.ui 
  vtktoolbarbase.ui 
  vtkpropertydialogbase.ui 
  drawtoolbarbase.ui 
  cursorwidgetbase.ui
  briconwidgetbase.ui)

FIND_PACKAGE( Qt REQUIRED )
FIND_PACKAGE( VTK REQUIRED )
FIND_PATH( BOOST_DIR boost/shared_ptr.hpp
	   /Developer/boost
	   /usr/local/boost
	   /usr/local/boost_1_33_0 )
FIND_PATH( QWT_DIR include/qwt_wheel.h 
	   /Developer/qwt 
	   /usr/local/qwt
	   /usr/local/qwt-4.2.0)
IF(NOT WINDOWS)
  FIND_PATH( FSL_DIR etc/fslconf/fsl.sh
	   /usr/local/fsl )
ENDIF(NOT WINDOWS)

ADD_DEFINITIONS( ${QT_DEFINITIONS} )

IF(VTK_DIR)
  INCLUDE( ${VTK_USE_FILE} )
  ADD_DEFINITIONS( -DHAVE_VTK )
ENDIF(VTK_DIR)

SET(HAS_SYSLOG ON CACHE BOOL "Use syslog found on many Unix systems")

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h )

INCLUDE_DIRECTORIES( ${QT_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${QWT_DIR}/include )
INCLUDE_DIRECTORIES( ${BOOST_DIR} )
INCLUDE_DIRECTORIES( ${FSL_DIR}/include )
INCLUDE_DIRECTORIES( ${FSL_DIR}/extras/include/newmat )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

LINK_DIRECTORIES( ${VTK_LIBRARY_DIRS} )
LINK_DIRECTORIES( ${QWT_DIR}/lib )
LINK_DIRECTORIES( ${FSL_DIR}/lib ${FSL_DIR}/extras/lib )
LINK_DIRECTORIES( ${FSLVIEW_BINARY_DIR}/lib )

ADD_LIBRARY( ui ${ui_srcs} ${moc_srcs} )
IF(APPLE)
  SET( MACOSX_BUNDLE_ICON_FILE application.icns )
  SET( MACOSX_BUNDLE_COPYRIGHT "Copyright (c) 2005, University of Oxford" )
  ADD_EXECUTABLE( fslview MACOSX_BUNDLE ${srcs} )
  INSTALL_TARGETS( /bin fslview.app )
ELSE (APPLE)
  ADD_EXECUTABLE( fslview ${srcs} )
  INSTALL_TARGETS( /bin fslview )
ENDIF(APPLE)

TARGET_LINK_LIBRARIES( fslview ui storage)

IF(UNIX)
  TARGET_LINK_LIBRARIES( fslview ${FSLVIEW_BINARY_DIR}/lib/libfslio.a )
ELSE(UNIX)
  IF(WINDOWS)
    LINK_DIRECTORIES( ${FSLVIEW_BINARY_DIR}/lib/Debug )
  ENDIF(WINDOWS)
  TARGET_LINK_LIBRARIES( fslview fslio )
ENDIF(UNIX)

TARGET_LINK_LIBRARIES( fslview niftiio znz z )
TARGET_LINK_LIBRARIES( fslview miscmaths newmat utils )
TARGET_LINK_LIBRARIES( fslview qwt )
TARGET_LINK_LIBRARIES( fslview QVTK vtkHybrid vtkWidgets)
