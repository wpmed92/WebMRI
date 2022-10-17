/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <algorithm>
#include "image.h"
#include <stdio.h>
#include "error.h"
#include <assert.h>

Image::Image(const std::string& filename):
  m_avw(NULL)
{
  if(!(m_avw = FslOpen(filename.c_str(), "r")))
    throw Exception("Failed to open file " + filename);
  m_imageInfo = ImageInfo::init(m_avw,filename); 
  m_cachedVolumes.clear();
  m_cachedVolumes.resize(m_imageInfo->inqNumVolumes());
}

Image::Image(ImageInfo::Handle info):
  m_avw(NULL), m_imageInfo(info)
{
  m_cachedVolumes.resize(info->inqNumVolumes());
 
  int x = m_imageInfo->inqX();
  int y = m_imageInfo->inqY();
  int z = m_imageInfo->inqZ();
  int v = m_imageInfo->inqNumVolumes();

  m_imageInfo->setMin(0);
  m_imageInfo->setMax(4);

  for(int n = 0; n < v; ++n)
    {
      switch(m_imageInfo->inqDt()) 
      {
       case DT_UNSIGNED_CHAR: m_cachedVolumes[n] = VolumeB::createBlank(x,y,z); break;
       case DT_SIGNED_SHORT:  m_cachedVolumes[n] = VolumeS::createBlank(x,y,z); break;
       case DT_SIGNED_INT:    m_cachedVolumes[n] = VolumeI::createBlank(x,y,z); break;
       case DT_FLOAT:         m_cachedVolumes[n] = VolumeF::createBlank(x,y,z); break;
       case DT_DOUBLE:        m_cachedVolumes[n] = VolumeD::createBlank(x,y,z); break;
      }
    }
}

Image::~Image()
{  
  if(m_avw != NULL){FslClose(m_avw);}
}

Image::Handle Image::load(const std::string& filename)
{
  Image::Handle im(new Image(filename));
//   if(im->getInfo()->inqLutName() == "")
//     im->getVolume(0)->calculateMinMax();
//   else if(im->getInfo()->inqMax() <= im->getInfo()->inqMin())
//     im->getVolume(0)->calculateMinMax();

  float min(im->getInfo()->inqMin());
  float max(im->getInfo()->inqMax());
  
  if(min == max) {
    switch(im->getInfo()->inqDt()) 
      {
      case DT_UNSIGNED_CHAR:
      case DT_SIGNED_SHORT:
      case DT_SIGNED_INT:
	if((im->getVolume(0)->minValue() >= -100) &&
	   (im->getVolume(0)->maxValue() <=  100))
	  {
	    im->getVolume(0)->calculateMinMax();
	    break;
	  }
      case DT_FLOAT:
      case DT_DOUBLE:
	if((min == 0) || (max == 0))
	  im->getVolume(0)->calculateRobustMinMax();
	else
	  im->getVolume(0)->calculateRobustMinMax(true);
	break;
      }
  }

  return im;
}

bool Image::save(const std::string& filename)
{
  FSLIO *avw = FslOpen(filename.c_str(), "w");

  if( !avw ) return false;

  m_imageInfo->saveAvwHeader(avw);

  int v = m_imageInfo->inqNumVolumes();
  for(int n = 0; n < v; ++n)
    {
      Volume::Handle v = getVolume(n);

      v->saveVolume(avw);
    }

  FslClose(avw);

  return true;
}

Image::Handle Image::cloneStructure()
{
  Image::Handle h = Image::Handle(new Image(m_imageInfo->clone()));
  h->getInfo()->setLutName("Red-Yellow");

  return h;
}

Image::Handle Image::clone3dStructure()
{
  ImageInfo::Handle i(m_imageInfo->clone());

  i->setNumVolumes(1);

  Image::Handle h = Image::Handle(new Image(i));
  h->getInfo()->setLutName("Red-Yellow");

  return h;
}

Volume::Handle Image::getVolume(short n) const
{
  n = std::max(0, std::min(int(n), m_imageInfo->inqNumVolumes() - 1));

  float min = m_imageInfo->inqMin();
  float max = m_imageInfo->inqMax();

  if(!m_cachedVolumes[n].get())
    {
      std::string name;
      
      switch(m_imageInfo->inqDt()) 
      {
      case DT_UNSIGNED_CHAR: 
	name="unsigned char"; 
	m_cachedVolumes[n] = VolumeB::getVolume(m_avw, min, max, n);
	break;
      case DT_SIGNED_SHORT:
	name="signed short"; 
	m_cachedVolumes[n] = VolumeS::getVolume(m_avw, min, max, n);
	break;
      case DT_SIGNED_INT: 
	name="signed int"; 
	m_cachedVolumes[n] = VolumeI::getVolume(m_avw, min, max, n);
	break;
      case DT_FLOAT:         
	name="float";
	m_cachedVolumes[n] = VolumeF::getVolume(m_avw, min, max, n);
	break;
      case DT_DOUBLE:
        name="double";
	m_cachedVolumes[n] = VolumeD::getVolume(m_avw, min, max, n);
	break;
      default:
	throw FileError(m_imageInfo->inqFileName(),"Data type not supported.");
      }
    }

  return m_cachedVolumes[n];
}

TimeSeries::Handle Image::getTimeSeries(short x, short y, short z) const
{
  short xDim,yDim,zDim;
  long  voxel;
  TimeSeriesMap::iterator i;
  
  xDim = m_imageInfo->inqX();
  yDim = m_imageInfo->inqY();
  zDim = m_imageInfo->inqZ();

  voxel = xDim*(y + z * yDim) + x;
  i =  m_cachedTimeSeries.find(voxel);

  if(m_cachedTimeSeries.end() == i)
  {     
      switch(m_imageInfo->inqDt()) 
      {
       case DT_UNSIGNED_CHAR: m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesB::getTimeSeries(m_avw,x,y,z))); break;
       case DT_SIGNED_SHORT:  m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesS::getTimeSeries(m_avw,x,y,z))); break;
       case DT_SIGNED_INT:    m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesI::getTimeSeries(m_avw,x,y,z))); break;
       case DT_FLOAT:         m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesF::getTimeSeries(m_avw,x,y,z))); break;
       case DT_DOUBLE:        m_cachedTimeSeries.insert(TimeSeriesMap::value_type(voxel, TimeSeriesD::getTimeSeries(m_avw,x,y,z))); break;
      }
  }
 return m_cachedTimeSeries[voxel];
}

Volume::Handle Image::blankDraw()
{
  VolumeS::Handle v = VolumeS::create(m_imageInfo->inqX(),
                                      m_imageInfo->inqY(),
                                      m_imageInfo->inqZ(),
                                      new short[m_imageInfo->inqX()*
                                                m_imageInfo->inqY()*
                                                m_imageInfo->inqZ()]);
  v->setMin(0);
  v->setMax(255);
  

  unsigned int width  = v->inqX();
  unsigned int height = v->inqY();
  unsigned int depth  = v->inqZ();

  for(unsigned int z = 0; z< depth; ++z)    
  {
    for(unsigned int y = 0; y < height; ++y)
    {      
      for(unsigned int x = 0; x < width; ++x)
        {
          (*v)(x,y,z) = 0;
        }
    }
  }
  
  return v;
}
