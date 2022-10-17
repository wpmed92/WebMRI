/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "tracker.h"

#include "filemanager.h"
#include "lookuptable.h"

#include <cstdio>
#include <ctype.h>

#include <iostream>
#include <fstream>
#include <string>
#include <limits>

FileManager::FileManager()
{ 
}

FileManager::~FileManager()
{
}

struct Output
{
  Output(ostream& os) : m_os(os) {}
  void operator()(const BaseCluster::Handle& h)
  {
    m_os << h << std::endl;
  }
  ostream& m_os;
};

//! @brief Read clusters from cluster_*.txt files in a FEAT directory
//! @param filename The filename in which the clusters should be found
//! @param clusters A handle for the @ref ClusterDataList which will hold all the data
void FileManager::readClusters(const std::string& filename, ClusterList& clusters)
{
  STATIC_TRACKER("FileManager::readClusters(const std::string&, ClusterList&)");
  clusters.clear();

  MESSAGE("Reading " + filename);
  std::ifstream cf(filename.c_str());

  if(!cf)
    throw std::ios::failure("FileManager::readClusters failed to open file: " + filename);

  std::list<string> headers;
				// Accumulate a list of header names
  string header;
  getline(cf, header);
  istringstream ss(header);
  string name;
  while(getline(ss, name, '\t'))
    {
      //      cout << name << endl;
      headers.push_back(name);
    }

  				// ...read the remaining cluster data.
  BaseCluster::Handle c = Cluster::create(headers);
  cf >> c;
  while(!cf.eof())
    {
      //     std::cout << c << std::endl;
      clusters.push_back(c);

      c = Cluster::create(headers);
      cf >> c;
   }
}

//! @brief Read clusters from cluster_*.txt files in a FEAT directory
//! @param filename The filename in which the clusters should be found
//! @param clusters A handle for the @ref ClusterDataList which will hold all the data
void FileManager::readTalairachClusters(const std::string& filename, ClusterList& clusters)
{
  STATIC_TRACKER("FileManager::readTalairachClusters(const std::string&, ClusterList&)");
  clusters.clear();

  MESSAGE("Reading " + filename);
  std::ifstream cf(filename.c_str());

  if(!cf)
    throw std::ios::failure("FileManager::readClusters failed to open file: " + filename);

  std::list<string> headers;
				// Accumulate a list of header names
  string header;
  getline(cf, header);
  istringstream ss(header);
  string name;
  while(getline(ss, name, '\t'))
    {
      //      cout << name << endl;
      headers.push_back(name);
    }

				// ...read the remaining cluster data.
  BaseCluster::Handle c = TalairachCluster::create(headers);
  cf >> c;
  while(!cf.eof())
    {
      //      std::cout << c << std::endl;
      clusters.push_back(c);

      c = TalairachCluster::create(headers);
      cf >> c;
    }
}

//! @brief Read a LUT style LookUpTable from file
//!
//! @param filename Name of the file in which the lookuptable can
//!        be found.
//! @param lut Pointer to a LookUpTable to be populated from the file
void  FileManager::readLutFile(const std::string& filename, LookUpTable* lut)
{
  float red, blue, green;
  char chLine[64] = "";
  int linenumber = 0,count = 0;
  
  FILE* fp = fopen(filename.c_str(),"r");
  // FILE* out = fopen("render3t.ldt","w");
  if(!fp)
    throw std::ios::failure("Couldn't open lut file!");

  fseek(fp,0,SEEK_SET);

  // Should start %!VEST-LUT

  fgets(chLine,64,fp);
  if(!strcmp(chLine, "%!VEST-LUT"))
    {
      fclose(fp);
      throw std::ios::failure("File is not a valid lut file! Should start: %!VEST-LUT");
    }

  while(fgets(chLine,64,fp)){if(strncmp(chLine,"<-color{",8) == 0) ++count;} 

  lut->allocateMemory(count);

  fseek(fp,0,SEEK_SET);

  while(fgets(chLine,64,fp))
  {
    if(strncmp(chLine,"<-color{",8) == 0)
      {
      sscanf(chLine + 8,"%f , %f , %f",&red,&green,&blue);
      lut->pushValue((unsigned char)(red*255.0),(unsigned char)(green*255.0), (unsigned char)(blue*255.0),
                     linenumber);
      //    fprintf(out,"\"%f,%f,%f\",\n",red,green,blue);
      ++linenumber;
      }    
  }
  fclose(fp);
  // fclose(out);
}  

//! @brief Read an RGB style LookUpTable from file
//!
//! @param filename Name of the file in which the lookuptable can
//!        be found.
//! @param lut Pointer to a LookUpTable to be populated from the file
void  FileManager::readRgbFile(const std::string& filename, LookUpTable* lut)
{
  int  index, red, blue, green,
       idxVal(0),idxValMax(0);
  char name[64]   = "";
  char chLine[128] = "";
  char idxStr[10] = "";
  char firstChar;

  FILE* fp = fopen(filename.c_str(),"r");
  if(!fp)
    throw std::ios::failure("Couldn't open rgb file!");

  fseek(fp,0,SEEK_SET);

  fgets(chLine,64,fp);
  if(!strcmp(chLine, "%!VEST-LUT"))
    {
      fclose(fp);
      throw std::ios::failure("File is not a valid rgb file!");
    }

  while(fgets(chLine,128,fp))
    {
      sscanf(chLine," %c",&firstChar);

      if(isdigit(firstChar))
        {      
          sscanf(chLine," %s",idxStr);
          idxVal = atoi(idxStr);
          if (idxVal > idxValMax)idxValMax = idxVal;
        }

    }

  lut->allocateMemory(idxValMax + 1);

  fseek(fp,0,SEEK_SET);

  while(fgets(chLine,128,fp))
  {
    sscanf(chLine," %c",&firstChar);
    if(isdigit(firstChar))
      {
        sscanf(chLine,"%i %s %i %i %i",&index,name,&red,&green,&blue);
        lut->pushValue((unsigned char)red,(unsigned char)green,(unsigned char)blue,index);
      }
  }

  fclose(fp);
}  

bool FileManager::checkFileExists(const std::string& path)
{
  bool exists(false);
  std::ifstream is(path.c_str());
  if(is)
    exists = true;
  return exists;
}
