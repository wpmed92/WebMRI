/*  lslaplacemanager.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(lslaplacemanager_h)
#define lslaplacemanager_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "miscmaths/volume.h"
#include "miscmaths/volumeseries.h"
#include "miscmaths/minimize.h"
#include "model.h"
#include "bintoptions.h"

using namespace NEWMAT;
using namespace MISCMATHS;

namespace Bint {

  class SumSquaresEvalFunction : public EvalFunction
  {
  public:
    
    SumSquaresEvalFunction(ForwardModel& pmodel, const ColumnVector& pdata, int pdebuglevel, bool pupdateprec, float pprec, bool panalmargprec) : EvalFunction(), model(pmodel),data(pdata),ntpts(data.Nrows()),updateprec(pupdateprec), prec(pprec), debuglevel(pdebuglevel), analmargprec(panalmargprec){}
    
    float evaluate(const ColumnVector& x) const;

    virtual ~SumSquaresEvalFunction(){};

  protected:

    ForwardModel& model;
    const ColumnVector& data;
    mutable int ntpts;
    bool updateprec;
    float prec;
    int debuglevel;
    bool analmargprec;

  private:

    SumSquaresEvalFunction();
    const SumSquaresEvalFunction& operator=(SumSquaresEvalFunction& par);
    SumSquaresEvalFunction(const SumSquaresEvalFunction&);
  };

class SumSquaresgEvalFunction : public gEvalFunction
  {
  public:
    
    SumSquaresgEvalFunction(gForwardModel& pmodel, const ColumnVector& pdata, int pdebuglevel, bool pupdateprec, float pprec, bool panalmargprec) : gEvalFunction(),model(pmodel),data(pdata),ntpts(data.Nrows()),updateprec(pupdateprec), prec(pprec), debuglevel(pdebuglevel), analmargprec(panalmargprec){}
    
    float evaluate(const ColumnVector& x) const;
    ReturnMatrix g_evaluate(const ColumnVector& x) const;

    virtual ~SumSquaresgEvalFunction(){};

  protected:

    gForwardModel& model;
    const ColumnVector& data;
    mutable int ntpts;
    bool updateprec;
    float prec;
    int debuglevel;
    bool analmargprec;

  private:

    SumSquaresgEvalFunction();
    const SumSquaresgEvalFunction& operator=(SumSquaresgEvalFunction& par);
    SumSquaresgEvalFunction(const SumSquaresgEvalFunction&);
  };

  class LSLaplaceVoxelManager
  {
  public:

    LSLaplaceVoxelManager(ForwardModel& pmodel,int pdebuglevel, bool panalmargprec) : model(pmodel), nparams(0), debuglevel(pdebuglevel), analmargprec(panalmargprec) 
    {
      evalfunction = new SumSquaresEvalFunction(model,data,debuglevel,updateprec,prec,analmargprec);
    }

    LSLaplaceVoxelManager(gForwardModel& pmodel,int pdebuglevel, bool panalmargprec) : model(pmodel), nparams(0), debuglevel(pdebuglevel), analmargprec(panalmargprec) 
    {
      evalfunction = new SumSquaresgEvalFunction(pmodel,data,debuglevel,updateprec,prec,analmargprec);
    }

    virtual ~LSLaplaceVoxelManager(){delete evalfunction;}    

    void run();
    void setupparams(float prec);
    virtual void setdata(const ColumnVector& pdata);
    float geterrorprecisionmean() const { 
      if(updateprec)
	return parammeans(nparams); 
      else
	return prec;
    }

    const ColumnVector& getparammeans() const {return parammeans;}
    const SymmetricMatrix& getparaminvcovs() const {return paraminvcovs;}

    const string& getparamname(int p) const {return model.getparam(p).getname();}
    int getnparams() const {return nparams;}
    int getnvaryingparams() const {return nvaryingparams;}
    int getntpts() const {return ntpts;}
 
  protected:

    ForwardModel& model;    

    int ntpts;
    int nparams;
    int nvaryingparams;
    int debuglevel;
    bool analmargprec;

    ColumnVector parammeans;
    SymmetricMatrix paraminvcovs;

    ColumnVector data;

    bool updateprec;
    float prec;
    EvalFunction* evalfunction;
    
  private:

    LSLaplaceVoxelManager();
    const LSLaplaceVoxelManager& operator=(LSLaplaceVoxelManager& par);     
    LSLaplaceVoxelManager(LSLaplaceVoxelManager& des);
  };

  class LSLaplaceManager
  {
  public:

    LSLaplaceManager(int pdebuglevel, float pprecin, bool panalmargprec, ForwardModel& pmodel, VolumeSeries& pdata, Volume& pmask) : 
      data(pdata),
      mask(pmask),      
      debuglevel(pdebuglevel),
      precin(pprecin),
      analmargprec(panalmargprec)
    {
	voxelmanager = new LSLaplaceVoxelManager(pmodel,debuglevel,analmargprec);
    }

    LSLaplaceManager(BintOptions& opts, ForwardModel& pmodel, VolumeSeries& pdata, Volume& pmask) : 
      data(pdata),
      mask(pmask),      
      debuglevel(opts.debuglevel.value()),
      precin(opts.prec.value()),
      analmargprec(opts.analmargprec.value())
    {
	voxelmanager = new LSLaplaceVoxelManager(pmodel,debuglevel,analmargprec);
    }

    LSLaplaceManager(int pdebuglevel, float pprecin, bool panalmargprec, gForwardModel& pmodel, VolumeSeries& pdata, Volume& pmask) : 
      data(pdata),
      mask(pmask),      
      debuglevel(pdebuglevel),
      precin(pprecin),
      analmargprec(panalmargprec)
    {
	voxelmanager = new LSLaplaceVoxelManager(pmodel,debuglevel,analmargprec);
    }

    LSLaplaceManager(BintOptions& opts, gForwardModel& pmodel, VolumeSeries& pdata, Volume& pmask) : 
      data(pdata),
      mask(pmask),      
      debuglevel(opts.debuglevel.value()),
      precin(opts.prec.value()),
      analmargprec(opts.analmargprec.value())
    {

      voxelmanager = new LSLaplaceVoxelManager(pmodel,debuglevel,analmargprec);
    }

    // load data from file in from file and set up starting values
    void setup();

    // run
    void run();

    // saves results in logging directory
    void save();     

    // Destructor
    virtual ~LSLaplaceManager() {}
 
  protected:     

    VolumeSeries& data;
    Volume& mask;

    int ntpts;    
    int nvoxels;
    int nparams;

    VolumeSeries mns;
    VolumeSeries covs;

    Volume prec;

    int debuglevel;
    float precin;
    bool analmargprec;

    LSLaplaceVoxelManager* voxelmanager;
    
  private:

    LSLaplaceManager();
    const LSLaplaceManager& operator=(LSLaplaceManager& par);     
    LSLaplaceManager(LSLaplaceManager& des);

  };
}   
#endif







