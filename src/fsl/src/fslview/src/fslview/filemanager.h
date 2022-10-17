/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(FILEMANAGER_H)
#define FILEMANAGER_H

#include <string>

#include "modelfit.h"
#include "clusterdata.h"

class LookUpTable;

//! @brief Collection of file utilities
class FileManager
{
public:

   FileManager();
   virtual ~FileManager();
   static void readLutFile(const std::string&, LookUpTable*);
   static void readRgbFile(const std::string&, LookUpTable*);
   static void readClusters(const std::string&, ClusterList&);
   static void readTalairachClusters(const std::string&, ClusterList&);

   static bool checkFileExists(const std::string& path);
};
#endif
