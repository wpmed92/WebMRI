/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "shape.h"
#include "slicewidget.h"
#include <algorithm>
#include <math.h>
#include <qobject.h>
#include <qpainter.h>
class DrawVertex {
public:
  DrawVertex(QPainter* p):m_paint(p)  {}
  void operator()(Pixel::Handle vo) 
  { 
        m_paint->setPen(QColor(255,255,255) );
        m_paint->drawPoint(vo->inqX(),vo->inqY());
        m_paint->drawRect(vo->inqX(),vo->inqY(),2,2);
  }
private:
  QPainter* m_paint;
};


class CommitVertex {
public:
  CommitVertex(Volume::Handle vol,int orient,int s) : m_volume(vol),
                                                      m_orient(orient),
                                                      m_slice(s)  {}
  void operator()(Pixel::Handle vo) 
  {  
    if(m_orient == SliceWidget::Axial)   m_volume->setValue(vo->inqX(),vo->inqY(),m_slice,vo->inqVal());
    if(m_orient == SliceWidget::Sagittal)m_volume->setValue(m_slice,vo->inqX(),vo->inqY(),vo->inqVal());
    if(m_orient == SliceWidget::Coronal) m_volume->setValue(vo->inqX(),m_slice,vo->inqY(),vo->inqVal());
  }
private:
  Volume::Handle m_volume;
  int m_orient;
  int m_slice;
};

class PixelRead {
public:
  PixelRead(Volume::Handle vol,int orient,int s,
            std::vector<Pixel::Handle>& v) :  m_volume(vol),
                                            m_orient(orient),
                                            m_slice(s),
                                            m_volPixels(v)
                                            {}
  void operator()(Pixel::Handle pix) 
  {  
      Pixel::Handle newPixel;

      if(m_orient == SliceWidget::Axial)
        {   
  
            if(pix->inqX() >= 0 && pix->inqX() < m_volume->inqX() &&
               pix->inqY() >= 0 && pix->inqY() < m_volume->inqY() &&
               m_slice     >= 0 && m_slice     < m_volume->inqZ())
              {
               newPixel = Pixel::create(pix->inqX(),pix->inqY(),m_volume->value(pix->inqX(),pix->inqY(), m_slice));  
               m_volPixels.push_back(newPixel);
              }
        }
      if(m_orient == SliceWidget::Sagittal)
        {
         
            if(m_slice     >= 0 && m_slice     < m_volume->inqX() &&
               pix->inqX() >= 0 && pix->inqX() < m_volume->inqY() &&
               pix->inqY() >= 0 && pix->inqY() < m_volume->inqZ()) 
              {
               newPixel = Pixel::create(pix->inqX(),pix->inqY(),m_volume->value(m_slice,    pix->inqX(), pix->inqY()));
               m_volPixels.push_back(newPixel);
              }
        }
      
      if(m_orient == SliceWidget::Coronal) 
        {            
            if(pix->inqX() >= 0 && pix->inqX() < m_volume->inqX() &&
               m_slice     >= 0 && m_slice     < m_volume->inqY() &&
               pix->inqY() >= 0 && pix->inqY() < m_volume->inqZ()) 
              { 
                newPixel = Pixel::create(pix->inqX(),pix->inqY(),m_volume->value(pix->inqX(),m_slice,     pix->inqY()));
                m_volPixels.push_back(newPixel);
              }
        }

    
  }
private:
  Volume::Handle m_volume;
  int m_orient;
  int m_slice;
  std::vector<Pixel::Handle>& m_volPixels;
};


class PixelSearch 
{
public:
  PixelSearch(Pixel::Handle & pix) : 
    m_found(false),m_newPix(pix){}
  void operator()(Pixel::Handle pix) 
  {
    if ((pix->inqX() == m_newPix->inqX()) &&
        (pix->inqY() == m_newPix->inqY()))
      {
        m_found = true; m_foundPix = pix;
      }
  }
  bool m_found;
  Pixel::Handle m_foundPix;
private:
  Pixel::Handle m_newPix;
};


Shape::Shape(QPainter* p,Volume::Handle vol, int orient, int slice):
             m_paint(p),       m_volume(vol),
             m_orient(orient), m_slice(slice){}
Shape::~Shape(){}

Shape::Handle Shape::create(QPainter* painter, Volume::Handle vol, int orient, int slice)
{
  Handle dst(new Shape(painter, vol,orient,slice));
  return dst;
}

void Shape::draw()
{
	std::copy(m_pixels.begin(),m_pixels.end(),std::back_inserter(m_commitPixels));
  std::for_each(m_pixels.begin(),m_pixels.end(),DrawVertex(m_paint));
  m_pixels.clear();
}

void Shape::commit()
{
  std::for_each(m_commitPixels.begin(),
                m_commitPixels.end(),
                CommitVertex(m_volume,m_orient,m_slice));

  m_commitPixels.clear();
}

void Shape::floodFill(int x, int y,float newVal)
{
  /*
    This flood algorithm has been taken from "CVu The Journal of the ACCU" 
    (August 2003 Volume 15 No 4) See www.accu.org for more details.
    The article was by James Holland.  The algorithm is very similar to the
    spans fill algorithm mentioned in "Foley and van Dam"(Computer Graphics
    Section 19.5.2 The Basic Filling Algorithms.
  */
  
  if(inRange(x,y))
  {
    float oldVal = readPixel(x,y);
  if(oldVal != newVal)
  {
    Location seed_location = {x,y};
    m_seedStack.push(seed_location);
    
    while(!m_seedStack.empty())
      {
        Location location = m_seedStack.top();
        Location locationOrig = m_seedStack.top();
        //Push pixel so that fill can be undone
        pushFloodUndoPixel(location.column, location.row,oldVal);
        writePixel(location.column, location.row,newVal);
        --location.column;
        while(location.column >= 0 && 
              readPixel(location.column,location.row) == oldVal)
          {            
            pushFloodUndoPixel(location.column, location.row,oldVal);
            writePixel(location.column, location.row,newVal);
            --location.column;
          }
        int extreme_left = location.column + 1;

        location.column = m_seedStack.top().column + 1;
        m_seedStack.pop();

        while(inRange(location.column,location.row) &&
              readPixel(location.column, location.row) == oldVal)
          {
            pushFloodUndoPixel(location.column, location.row,oldVal);
            writePixel(location.column, location.row,newVal);
            ++location.column;
          }
        int extreme_right = location.column - 1;
        
        //Scan above the seed row

        if(inRange(locationOrig.column,locationOrig.row + 1))
          {        
            location.row = locationOrig.row + 1;

            bool previous_pixel_is_border = true;
            for (location.column = extreme_right;
                 location.column>= extreme_left;
                 --location.column)
              {
                if(previous_pixel_is_border &&
                   readPixel(location.column, location.row) == oldVal)
                  {
                    m_seedStack.push(location);
                    previous_pixel_is_border = false;
                  }
                else
                  if(readPixel(location.column, location.row) != oldVal)
                    {
                      previous_pixel_is_border = true;
                    }
                  }
              }
        //Scan below the seed row
        location.row = locationOrig.row - 1;
        
        if(location.row >= 0)
          {
            bool previous_pixel_is_border = true;
            for (location.column = extreme_right;
                 location.column>= extreme_left;
                 --location.column)
              {
                if(previous_pixel_is_border &&
                   readPixel(location.column, location.row) == oldVal)
                  {
                    m_seedStack.push(location);
                    previous_pixel_is_border = false;
                  }
                else
                  if(readPixel(location.column, location.row) != oldVal)
                    {
                      previous_pixel_is_border = true;
                    }
              }
          }
      }
  }
  }
  
}




void Shape::addVertex(int x, int y, int size, float val)
{  
  /*
    When a vertex is added a line of pixels is automatically drawn
    between the last vertex and the current vertex.  The basic concept
    of how this is acheived is roughly described in Foley and van Dam 
    section 3.2.1 "The Basic Incremental Algoritm"

    Different processes occur depending on wether the gradient of the line
    is less than or bigger than 1. This test is used to decide wether to
    increment x and calculate y or increment y and calculate x.
  */

  Pixel::Handle cur = Pixel::create(x,y,val);
 
  if(m_commitPixels.empty())
   {
      addSurroundingPixels(cur,size,val);
   }
  else
  {
    Pixel::Handle prev = m_commitPixels.back();
    
    int diffX = cur->inqX() - prev->inqX();
    int diffY = cur->inqY() - prev->inqY();

    if((diffX != 0)||(diffY != 0))
      {
        float grad,gradrecip;
        if(diffY == 0)
          {
            grad = 0;gradrecip = 999;
          }
        else if(diffX == 0)
          {
            grad = 999;gradrecip = 0;
          }
        else
          {
            grad = float(diffY)/float(diffX);
            gradrecip = 1.0/grad;
          }
        //increment x if grad is less than 1.0
        if(fabs(grad) <= 1.0)
          {
            float yValue;
            int xStep(0),yValueInt,yValueSign;
            while(abs(xStep)<= abs(diffX))
              {
                yValue     = grad * xStep + prev->inqY();

                if(yValue < 0){yValueSign = -1;}else{yValueSign = 1;}

                yValueInt  = int(floor(fabs(yValue + 0.5))) * yValueSign;

                Pixel::Handle mid = Pixel::create(xStep + prev->inqX(),
                                                  yValueInt,
                                                  val);     
                addSurroundingPixels(mid,size,val);

                if(diffX < 0){--xStep;}
                else         {++xStep;}
              }

          }
        else
          {
            //increment y if grad is bigger than 1.0
            float xValue;
            int yStep(0),xValueInt,xValueSign;
            while(abs(yStep) <= abs(diffY))
              {                
                xValue     = gradrecip * yStep + prev->inqX();  
                if(xValue < 0){xValueSign = -1;}else{xValueSign = 1;}
                xValueInt  = int(floor(fabs(xValue + 0.5))) * xValueSign;                        
                Pixel::Handle mid = Pixel::create(xValueInt,
                                                  yStep + prev->inqY(),
                                                  val); 
                addSurroundingPixels(mid,size,val);

                if(diffY < 0){--yStep;}
                else         {++yStep;}
              }

          }
      }
  }
}


Shape::Handle Shape::getBuffer()
{
  Shape::Handle hnd; 
  hnd = Handle(new Shape(m_paint,m_volume,m_orient,m_slice));
  std::for_each(m_commitPixels.begin(),
                m_commitPixels.end(),
                PixelRead(m_volume,m_orient,m_slice,hnd->m_commitPixels));
  
  return hnd; 
}

Shape::Handle Shape::getFloodBuffer()
{
  Shape::Handle hnd; 
  hnd = Handle(new Shape(m_paint,m_volume,m_orient,m_slice));
  hnd->m_commitPixels = m_floodUndoPixels;

  return hnd; 
}






Pixel::Handle Pixel::create(int tx, int ty, float tval)
{
  Handle dst(new Pixel(tx,ty,tval));
  return dst;
}

Pixel::Pixel(int tx,int ty, float tval):m_x(tx),m_y(ty),
                                        m_val(tval),m_drawn(false)
{

}

float Shape::readPixel(int x, int y)
{
  float value(0.0);
  
  if(m_orient == SliceWidget::Axial)
             value = m_volume->value(x,y, m_slice); 
  if(m_orient == SliceWidget::Sagittal)
             value =  m_volume->value(m_slice, x,y);
  if(m_orient == SliceWidget::Coronal)            
             value = m_volume->value(x,m_slice,y);

  return value;
}

void Shape::writePixel(int x, int y,float newVal)
{ 
  if(m_orient == SliceWidget::Axial) 
    m_volume->setValue(x,y,m_slice,newVal);

  if(m_orient == SliceWidget::Sagittal)
    m_volume->setValue(m_slice,x,y,newVal);

  if(m_orient == SliceWidget::Coronal) 
    m_volume->setValue(x,m_slice,y,newVal);
}

void Shape::pushFloodUndoPixel(int x, int y,float oldVal)
{
  Pixel::Handle cur = Pixel::create(x,y,oldVal);
  m_floodUndoPixels.push_back(cur);
}

bool Shape::inRange(int x, int y)
{
 bool  result(false);
 if(m_orient == SliceWidget::Axial)
       {   
  
         if(x >= 0 && x < m_volume->inqX() &&
            y >= 0 && y < m_volume->inqY() &&
            m_slice     >= 0 && m_slice     < m_volume->inqZ())
           {
             result = true;
           }
        }
  if(m_orient == SliceWidget::Sagittal)
        {
         
          if(m_slice     >= 0 && m_slice     < m_volume->inqX() &&
             x >= 0 && x < m_volume->inqY() &&
             y >= 0 && y < m_volume->inqZ()) 
            {
              result = true;
           }
        }
      
  if(m_orient == SliceWidget::Coronal) 
    {            
      if(x >= 0 && x < m_volume->inqX() &&
         m_slice     >= 0 && m_slice     < m_volume->inqY() &&
         y >= 0 && y < m_volume->inqZ()) 
        { 
          result = true;
        }
    }

  return result;
}
  

void Shape::addSurroundingPixels(Pixel::Handle & pix, int size, float val)
{
  int x = pix->inqX();
  int y = pix->inqY();

  for(int nx = x; nx < x + size;nx++)
    {
    for(int ny = y; ny < y + size;ny++)
      {
        Pixel::Handle mid = Pixel::create(nx - size/2,ny - size/2,val);
        push_check(mid,size);
      }
    }

  Pixel::Handle mid = Pixel::create(pix->inqX(),pix->inqY(),val);
  m_pixels.push_back(mid); 
}

void Shape::push_check(Pixel::Handle& pix, int size)
{
  std::vector<Pixel::Handle>::iterator start;
  int listsize = m_pixels.size();
  int checkAmount = size * size;
  if(listsize < checkAmount)
    {
      start = m_pixels.begin();
    }
  else
    {
      start =(std::vector<Pixel::Handle>::iterator) &m_pixels[listsize - checkAmount];
    }

  PixelSearch search = std::for_each(start,
                                     m_pixels.end(),
                                     PixelSearch(pix));

  if(!search.m_found)m_pixels.push_back(pix);
}
