#ifndef _abcshadow_
#define _abcshadow_ 1

#include "mydefines.hpp"
#include "hugmodel.hpp"
// #include"standardpattern.hpp"

class ABCshadowHug
{
public:
  // log(theta) search window
  double *list_bound_min_log; // numberOfComponent
  double *list_bound_max_log; // numberOfComponent
  // step maximal variation of log(theta)
  double *list_delta_log; // numberOfComponent
  // Current theta ans potential next theta
  double logtheta[numberOfComponent];           // numberOfComponent
  double logtheta_candidate[numberOfComponent]; // numberOfComponent
  // ABCshadow parameters
  int number_of_loop; // number of time the loop is executed
  int number_of_step; // number of steps in each generation of control realization;
  // MH parameters
  int time_mh;
  double probability_birth;
  double probability_death;
  double probability_changein;
  double probability_changeout;
  double radiuschange;

  // Statisitics
  double observed_statistics[numberOfComponent];        // numberOfComponent
  double current_pattern_statistics[numberOfComponent]; // numberOfComponent

  HugModel *desired_model; // desired model for which the parameters are computed
  Data *desired_datas;
  Sources *current_sources;
  int savingtime;

  char *sourcesfile;
  char *statisticsfile;
  char *sourcesfilefinal;
  char *statisticsfilefinal;
  char *parametersfile;

  ABCshadowHug(double list_bound_min_log_parameter[],                                                                                                                                 // numberOfComponent
               double list_bound_max_log_parameter[],                                                                                                                                 // numberOfComponent
               double list_delta_log_parameter[],                                                                                                                                     // numberOfComponent
               double list_initial_log_parameter[],                                                                                                                                   // numberOfComponent
               double observed_statistics_parameter[],                                                                                                                                // numberOfComponent
               int number_of_loop_parameter,                                                                                                                                          // ABCshadow main loop
               int number_of_step_parameter,                                                                                                                                          // ABCshadow nb parameter updates by loop
               int time_mh_parameter,                                                                                                                                                 // MH nb iteration in MH to generate a realization
               double probability_birth_parameter, double probability_death_parameter, double probability_changein_parameter, double probability_changeout_parameter, double p_radiuschange, // MH birth, death
               Data *p_datas, HugModel *p_themodel, Sources *p_sources,
               int p_savingtime,
               char *p_sourcesfile, char *p_statisticsfile, char *p_sourcesfilefinal, char *p_statisticsfilefinal, char *p_parameters);

  void compute_parameters();
  void generation_realization();
  void generation_new_parameter_candidates();
  double acceptance_threshold_computation();
  void update_parameters();
  void save_current_theta();
  void save_in_file_last_realization();
};

#endif
