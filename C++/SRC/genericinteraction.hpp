#ifndef _genericinteraction_
#define _genericinteraction_ 1

#include "point.hpp"
// #include "hugmh.hpp"
//class HugMH;
//#include "hugmodel.hpp"
class HugModel;
typedef std::vector<double> Coordinates;

class GenericInteraction
{
public:
  double parameter_model;
  double exp_parameter_model;
  int rank;

  GenericInteraction(double parameter_model_parameter, int rank_parameter);

  void setParameterInteraction(double parameter_interaction_parameter) { parameter_model = parameter_interaction_parameter; }
  double getParameterInteraction() { return parameter_model; }

  void setExpParameterInteraction(double exp_parameter_interaction_parameter) { exp_parameter_model = exp_parameter_interaction_parameter; }
  double getEXpParameterInteraction() { return exp_parameter_model; }

  virtual double cond_intens_birth(vector<Coordinates> proposed_sources,
                                   HugModel *the_hugmodel,
                                   vector<vector<Point>> *p_used_sources,
                                   vector<vector<Point>> *p_current_convex_hull,
                                   vector<vector<Point>> *p_next_convex_hull,
                                   vector<vector<Point>> *elt_usefulDatas) = 0;
  virtual double cond_intens_death(vector<Coordinates> proposed_sources,
                                   HugModel *the_hugmodel,
                                   vector<vector<Point>> *p_used_sources,
                                   vector<vector<Point>> *p_current_convex_hull,
                                   vector<vector<Point>> *p_next_convex_hull,
                                   vector<vector<Point>> *elt_usefulDatas) = 0;
  virtual double cond_intens_change(vector<Coordinates> proposed_sources,
                                    HugModel *the_hugmodel,
                                    vector<vector<Point>> *p_used_sources,
                                    vector<vector<Point>> *p_current_convex_hull,
                                    vector<vector<Point>> *p_next_convex_hull,
                                    vector<vector<Point>> *elt_usefulDatas) = 0;
  virtual void validationOfBirthEvolution() = 0;
  virtual void noValidationOfBirthEvolution() = 0;
  virtual void validationOfDeathEvolution() = 0;
  virtual void noValidationOfDeathEvolution() = 0;
  virtual void validationOfChangeEvolution() = 0;
  virtual void noValidationOfChangeEvolution() = 0;

  virtual void computeStatistic(HugModel *the_hugmodel,
                                vector<vector<Point>> *p_used_sources,
                                vector<vector<Point>> *p_current_convex_hull,
                                vector<vector<Point>> *elt_usefulDatas) = 0;

  virtual double getCurrentStatistic(int indice_plan) = 0;

  // virtual void decoreNewEvent( DecoratedEvent* a_decorated_event,int i ) { cout << "HUM ModelInteraction::decoreNewEvent !!!! " ; };

  // virtual double getCurrentStatistic() = 0;
  // virtual void computeStatistic( HugPattern& p ) = 0;
  // virtual void computeStatisticInCropWindow( HugPattern& p ) = 0;
  // virtual void printStatistic( HugPattern& p, char *name_file_parameter ) = 0;
  // virtual void printStatisticInCropWindow( HugPattern& p, char *name_file_parameter ) = 0;
};
#endif
