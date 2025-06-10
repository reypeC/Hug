#ifndef _hugmh_
#define _hugmh_ 1

#include "data.hpp"
#include "sources.hpp"
#include "hugmodel.hpp"
#include "point.hpp"
#include "hull.hpp"

typedef std::vector<double> Coordinates;

class HugMH
{
private :
  
public :
  double timeMH;  // number of application of MH
  int t;       // Current time
  double timeUpdate; // nb of update in each iteration
  int numberiteration; // nb of update in each iteration

  double temperature;
  double temperaturemin;
  double coolingcoeff;
  int savingtime;

  char * sourcesfile;
  char *statisticsfile;
  char *sourcesfilefinal;
  char *statisticsfilefinal;
  
  double pbirth;
  double pdeath;
  double pchangein;
  double pchangeout;
  double radiuschange;
  
  // Model
  HugModel* theModel;
  
  Sources* theSources;
  vector<Coordinates> proposed_sources;
  
  // Data
  Data* theData;
  int number_elements;  // number of element in data
  int number_vectors;  // number of row of data

  //vector <Point> elt_usefulDatas;  // X_a Y_b useful current data
  vector<vector <Point>> used_sources;  // X_a Y_b used sources
  
  // Convexes hulls
  vector<vector <Point>> current_convex_hull; 
  vector<vector <Point>> next_convex_hull;
  vector<vector <Point>> auxiliary_convex_hull; 
  
  HugMH (Data* p_datas, HugModel* p_themodel, Sources* p_sources,int p_numberiteration,
                  double p_timeMH,
                  double p_pbirth, double p_pdeath,double p_pchange_in, double p_pchange_out, double p_radiuschange,
                  double p_temperature, double p_temperaturemin, double p_coolingcoeff, int p_savingtime,
                  char *p_sourcesfile, char *p_statisticsfile, char *p_sourcesfilefinal, char *p_statisticsfilefinal);
  void birth();
  void death();
  void change();
  void changeInSource();
  void changeOutSource();
  void sim();
  void update();
  
};

#endif

