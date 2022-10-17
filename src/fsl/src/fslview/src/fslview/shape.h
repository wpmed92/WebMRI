/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(SHAPE_H)
#define SHAPE_H


#include <boost/shared_ptr.hpp>
#include <list>
#include <vector>

#include <stack>
#include "storage/volume.h"

class QPainter;
class PixelRead;


struct Location
{
  int column;
  int row;
};




class Pixel
{

public:   
  typedef boost::shared_ptr< Pixel > Handle;  

  static Handle create(int tx,int ty, float tval);   
  
  virtual ~Pixel(){};
  int inqX(){return m_x;}
  int inqY(){return m_y;}
  float inqVal(){return m_val;}
  void setDrawn(bool state){m_drawn = state;}
  bool inqDrawn(){return m_drawn;}

private: 
  Pixel(int tx, int ty, float tval);
  
  int m_x;
  int m_y;
  float m_val;
  bool m_drawn;
};

class Shape
{

public:

  typedef boost::shared_ptr< Shape > Handle;
  static Handle create(QPainter * p, Volume::Handle, int orient,int slice);
  virtual  ~Shape();
  void draw();   
  void commit();
  void list();
  bool empty(){return m_commitPixels.empty();}
  int  size() {return m_commitPixels.size();}
  Shape::Handle getBuffer();  
  Shape::Handle getFloodBuffer();
  virtual void addVertex(int x,int y, int size, float val);
  virtual void floodFill(int x, int y,float newVal);
  virtual float readPixel(int x,int y);
  virtual void writePixel(int x,int y,float newVal);
  virtual void pushFloodUndoPixel(int x, int y,float oldVal);
  virtual bool inRange(int x, int y);

private:
  Shape(QPainter* p, Volume::Handle vol, int orient, int slice);  
  void addSurroundingPixels(Pixel::Handle &,int size, float val);
  std::vector<Pixel::Handle> m_pixels;
  std::vector<Pixel::Handle> m_commitPixels;
  std::vector<Pixel::Handle> m_floodUndoPixels;
  QPainter* m_paint;
  void push_check(Pixel::Handle & pix, int size);
  Volume::Handle m_volume;
  int m_orient;
  int m_slice;
  int m_counter;
  std::stack<Location> m_seedStack;
};



#endif 
