/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if defined(WIN32) 
#pragma warning (disable:4786)
#endif
 
//#include "qwindowsstyle.h"

#include "application.h"
#include "splashscreen.h"

#include <qapplication.h>

namespace Utilities {

  bool string_to_T(std::pair<float,float> &p, const string& s) {
    string str(s), delin(",");
    std::vector<float> vf(0);
    if(str.find(":")!=string::npos)
      delin = ":";
    str=str+delin;
    vf.clear();
    while(str.size()) {
      float v = atof(str.substr(0,str.find(delin)).c_str());
      vf.push_back(v);
      str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
    }
    bool retval(false);
    if(vf.size() == 2) {
      p.first = vf[0];
      p.second = vf[1];
      retval = true;
    }
    return true;
  }
}

int main( int argc, char **argv )
{
  using namespace Utilities;

  try {
    QApplication::setColorSpec( QApplication::CustomColor );
    QApplication app(argc,argv);			

    Option<bool> verbose(string("-V,--verbose"), false, 
			 string("switch on diagnostic messages"), 
			 false, no_argument);
    Option<bool> help(string("-h,--help"), false,
		      string("display this message"),
		      false, no_argument);
    Option<string> mode(string("-m,--mode"), string("ortho"),
			string("Initial viewer mode. One of: 3d; ortho; lightbox"), false,
			requires_argument);
    string title("fslview (Version 2.4pre0)\n\nCopyright(c) 2005, University of Oxford\nDave Flitney");
    string usage("fslview [-m 3d|ortho|lightbox] <baseimage> [-l lutname] [-b low,hi] [ <overlay> [-l lutname] [-b low,hi] ] ...");
  
    OptionParser options(title, usage);

    OverlayOptionList overlays;

    options.add(verbose);
    options.add(help);
    options.add(mode);

				// The next two options are place
				// holders so that the usage and
				// option descriptions work. They
				// aren't actually used during command
				// line processing!
    Option<string> lutname(string("-l,--lut"), string("Unset"),
			   string("Lookup table name. As per GUI, one of: Greyscale; \"Red-Yellow\"; \"Blue-Lightblue\"; Red; Green; Blue; Yellow; Pink; Hot; Cool; Copper, etc."), 
			   false, requires_argument);
    Option< std::pair<float,float> > ibricon(string("-b,--bricon"), std::pair<float,float>(0.0,0.0), 
					     string("Initial bricon range, e.g., 2.3,6"), 
					     false, requires_argument);

    options.add(lutname);
    options.add(ibricon);

    for(unsigned int pos = options.parse_command_line(qApp->argc(), qApp->argv());
	int(pos) < qApp->argc(); ) {
      // Should be an image name followed by image sub options
      string filename(qApp->argv()[pos]);

      Option<string> lutname(string("-l,--lut"), string("Unset"),
			     string("Lookup table name. As per GUI, one of: Greyscale; \"Red-Yellow\"; \"Blue-Lightblue\"; Red; Green; Blue; Yellow; Pink; Hot; Cool; Copper, etc."), 
			     false, requires_argument);
      Option< std::pair<float,float> > ibricon(string("-b,--bricon"), std::pair<float,float>(0.0,0.0), 
					       string("Initial bricon range, e.g., 2.3,6"), 
					       false, requires_argument);

      OptionParser imageOptions("Per-image options", "image [-l GreyScale] [-b 2.3,6]");

      imageOptions.add(lutname);
      imageOptions.add(ibricon);

      pos += imageOptions.parse_command_line(qApp->argc() - pos, &(qApp->argv()[pos]));

      if(!imageOptions.check_compulsory_arguments())
	imageOptions.usage();

      overlays.push_back(OverlayOption(filename, lutname, ibricon));
    }

    if(help.value() || !options.check_compulsory_arguments())
      options.usage();

    else {
      app.connect( &app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()) );
      //  a.setStyle(new QWindowsStyle);
  
      SplashScreen *s = new SplashScreen(0, overlays);
      s->show();

      return app.exec();
    }

  } catch(X_OptionError& e) {
    //     options.usage();
    cerr << e.what() << endl;
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } catch (...) {
    cerr << "Unhandled exception!" << endl;
  }

  return -1;
}
