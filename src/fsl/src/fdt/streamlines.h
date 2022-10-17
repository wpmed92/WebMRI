#include <fstream>
#include "newimage/newimageall.h"
#include "utils/log.h"
#include "meshclass/meshclass.h"
#include "probtrackxOptions.h"
#include "particle.h"
#include "tractvolsx.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace TRACTVOLSX;
using namespace mesh;
using namespace PARTICLE;

namespace TRACT{
  void read_masks(vector<string>& masks, const string& filename);


  class Streamliner{
    //Everything in DTI space is done INSIDE this class and lower level classes (particle and tractvolsx)
    //This class communicates with higher level classes in Seed voxels.
    //
    probtrackxOptions& opts;
    Log& logger;
    Particle m_part;
    vector<ColumnVector> m_path;
    volume<int> m_mask;
    volume<int> m_skipmask;
    volume<int> m_rubbish;
    volume4D<float> m_loopcheck;
    vector<volume<int>* > m_waymasks;
    vector<bool> m_passed_flags;
    vector<bool> m_own_waymasks;
    Matrix m_Seeds_to_DTI;
    Tractvolsx vols;
    float m_lcrat;
    float m_x_s_init;
    float m_y_s_init;
    float m_z_s_init;
  public:
    //Constructors
    Streamliner();
    ~Streamliner(){
      for(unsigned int i=0;i<m_waymasks.size();i++)
	if(m_own_waymasks[i]) delete m_waymasks[i];
    }
    void add_waymask(volume<int>& myway,const bool& ownership=false){
      //make sure that the waymask to add will not be deleted before
      //this streamliner goes out of scope!!
      m_waymasks.push_back(&myway);
      m_own_waymasks.push_back(ownership);
      m_passed_flags.push_back(false);
    }
    void pop_waymasks(){
      volume<int>* tmpptr=m_waymasks[m_waymasks.size()-1];
      m_waymasks.pop_back();
      m_passed_flags.pop_back();
      if(m_own_waymasks[m_own_waymasks.size()-1]){
	delete tmpptr;
      }
      m_own_waymasks.pop_back();
      
    }
    
    inline const float get_x_seed() const {return m_x_s_init;}
    inline const float get_y_seed() const {return m_y_s_init;}
    inline const float get_z_seed() const {return m_z_s_init;}
    inline const vector<ColumnVector>& get_path_ref() const{return m_path;}
    inline vector<ColumnVector> get_path() const{return m_path;}
    inline void reset(){
      m_part.reset();
    }
    inline void reverse(){
      m_part.restart_reverse();
    }
    bool streamline(const float& x_init,const float& y_init, const float& z_init,const ColumnVector& dim_seeds,const int& fibst);

  };


  class Counter{
    probtrackxOptions& opts;
    Log& logger;
    volume<int> m_prob;
    volume<int> m_beenhere;
    Matrix m_I;
    vector<ColumnVector> m_path;
    
    vector<volume<int> > m_seedcounts;
    vector<volume<int> > m_targetmasks;
    vector<string> m_targetmasknames;
    vector<int> m_targflags;
    //vector<vector<int> > m_particle_numbers;

    
    volume<int> m_ConMat;
    volume<int> m_CoordMat;
    int m_Conrow;

    volume<int> m_ConMat2;
    volume<int> m_CoordMat2;
    volume<int> m_CoordMat_tract2;
    volume<int> m_lrmask;
    volume4D<int> m_lookup2;
    volume<int> m_beenhere2;
    int m_Conrow2;
    ColumnVector m_lrdim;
    
    const volume<int>& m_seeds;
    ColumnVector m_seedsdim;
    const Streamliner& m_stline;
    Streamliner& m_nonconst_stline;
    
  public:
    Counter(const volume<int>& seeds,Streamliner& stline):opts(probtrackxOptions::getInstance()),
							  logger(LogSingleton::getInstance()),
							  m_seeds(seeds),m_stline(stline),
							  m_nonconst_stline(stline){
      //are they initialised to zero?
      m_beenhere.reinitialize(m_seeds.xsize(),m_seeds.ysize(),m_seeds.zsize());
      m_seedsdim.ReSize(3);
      m_seedsdim << m_seeds.xdim() <<m_seeds.ydim() <<m_seeds.zdim();
      m_I=Identity(4);
      
    }
    
    void initialise();
    
    void initialise_path_dist(){
      m_prob=m_seeds;
      m_prob=0;
    }
    void initialise_seedcounts();
    
    void initialise_matrix1(); //Need to make sure that initialise_path_dist is run first
    
    void initialise_matrix2();
    
    void initialise_maskmatrix(){} //not written yet
    
    inline void store_path(){ m_path=m_stline.get_path();}
    
    void count_streamline();
    void count_seed();
    void clear_streamline(const bool& forwardflag,const bool& backwardflag);
    
    
    void update_pathdist();
    void reset_beenhere(const bool& forwardflag,const bool& backwardflag);
    
    void reset_prob(){m_prob=0;}
    void update_seedcounts();
    void reset_targetflags(){
      for(unsigned int i=0;i<m_targflags.size();i++) m_targflags[i]=0;
    }
    
    
    void update_matrix1(); //update path_dist after each streamline, only run this after each voxel!!
    
    void update_matrix2_row(); //but run this one every streamline as with the others
    void next_matrix2_row(){m_Conrow2++;}//and then run this after each voxel..
    void reset_beenhere2(const bool& forwardflag,const bool& backwardflag);
  
    void update_maskmatrix(){} //not written yet
    
    void save();
    void save_pathdist();
    void save_pathdist(string add);
    void save_seedcounts();
    void save_matrix1();
    void save_matrix2();
    void save_maskmatrix(){}//not written yet
    

    inline const Streamliner& get_streamline() const {return m_stline;}
    inline Streamliner& get_nonconst_streamline() const {return m_nonconst_stline;}
    inline const volume<int>& get_seeds() const {return m_seeds;}

    
  };
  
  class Seedmanager{
    probtrackxOptions& opts;
    Log& logger;
    Counter& m_counter;    
    Streamliner& m_stline;
    const volume<int>& m_seeds;
    ColumnVector m_seeddims;
  public:
    Seedmanager(Counter& counter):opts(probtrackxOptions::getInstance()),
				  logger(LogSingleton::getInstance()),
				  m_counter(counter),
				  m_stline(m_counter.get_nonconst_streamline()),
				  m_seeds(m_counter.get_seeds()){
      m_seeddims.ReSize(3);
      m_seeddims<<m_seeds.xdim()<<m_seeds.ydim()<<m_seeds.zdim();
    }
    void run(const float& x,const float& y,const float& z,bool onewayonly=false, int fibst=-1);
  };

}
