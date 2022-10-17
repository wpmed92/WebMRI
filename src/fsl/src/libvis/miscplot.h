#if !defined (MISCPLOT_H)
#define MISCPLOT_H

#include <string>
#include "newmatap.h"
#include "miscmaths/histogram.h" 
#include "gd.h"

#ifndef FALSE
#define FALSE false
#endif
#ifndef TRUE
#define TRUE true
#endif


namespace MISCPLOT{
  
  class miscplot
    {
    public:
      
      //constructor
      miscplot(){req_xsize=0; req_ysize=0;explabel=string("");null_shift=0.0;};
      
      void timeseries(const NEWMAT::Matrix& mat, string filename, string title,
		      float tr = 0.0, int ysize = 150, int width = 4, int prec = 2, bool sci = false);

      void histogram(const NEWMAT::Matrix& mat, string filename, string title);

      void gmmfit(const NEWMAT::Matrix& mat, Matrix& mu, Matrix& sig, Matrix& pi, 
                  string filename, string title, bool mtype = false, 
                  float offset = 0.0, float detailfactor = 0.0); 
    
      inline void ggmfit(const NEWMAT::Matrix& mat, Matrix& mu, Matrix& sig, Matrix& pi, 
		       string filename, string title, float offset = 0.0, float detailfactor = 0.0)
	{this->gmmfit(mat, mu, sig, pi, filename, title, true, offset, detailfactor);}
      

      inline void add_label(string txt){ labels.push_back(txt);}
      inline void add_xlabel(string txt){ xlabels.push_back(txt);}
      inline void add_ylabel(string txt){ ylabels.push_back(txt);}
      inline void set_xysize(int xsize, int ysize){ req_xsize=xsize; req_ysize=ysize;}
      inline void set_nullshift(double val){null_shift=val;};

    private:
      
      vector<string> labels;
      vector<string> xlabels;
      vector<string> ylabels;
      string explabel;
      int req_xsize,req_ysize;
      double null_shift;
      
      gdImagePtr outim;

      void add_legend(void* ptr, unsigned long cmap[], bool inside=false);
    };
}
#endif
