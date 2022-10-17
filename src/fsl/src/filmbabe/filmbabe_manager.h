/*  filmbabe_manager.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  COPYRIGHT  */


#if !defined(filmbabe_manager_h)
#define filmbabe_manager_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "filmbabeoptions.h"
#include "newimage/newimageall.h"
#include "connected_offset.h"
#include "filmbabe_vb_flobs.h"

using namespace NEWIMAGE;

namespace Filmbabe {   
  
  // Give this class a file containing
  class Filmbabe_Manager
    {
    public:

      // constructor
      Filmbabe_Manager(FilmbabeOptions& popts) :
	opts(popts),
	connected_offsets()
	{ 
	  connected_offsets.push_back(Connected_Offset(-1,0,0,0,1));
	  connected_offsets.push_back(Connected_Offset(1,0,0,1,0));
	  connected_offsets.push_back(Connected_Offset(0,-1,0,2,3));
	  connected_offsets.push_back(Connected_Offset(0,1,0,3,2));
	  connected_offsets.push_back(Connected_Offset(0,0,-1,4,5));
	  connected_offsets.push_back(Connected_Offset(0,0,1,5,4));

	}

      // load data from file in from file and set up starting values
      void setup();

      // runs the chain
      void run();

      // saves results in logging directory
      void save();

      // getters
      const volume<int>& getMask() const {return mask;}

      // Destructor
      virtual ~Filmbabe_Manager();
 
    private:

      Filmbabe_Manager();
      const Filmbabe_Manager& operator=(Filmbabe_Manager& par);     
      Filmbabe_Manager(Filmbabe_Manager& des);

      volume<int> mask;
      volume4D<float> data;
      Matrix designmatrix;
      ColumnVector flobsregressors;

      Filmbabe_Vb_Flobs* filmbabe_vb_flobs;

      FilmbabeOptions& opts;

      volume4D<float> localweights;

      vector<Connected_Offset> connected_offsets;

      volumeinfo volinfo;

    };
}   
#endif

