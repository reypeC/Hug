#include "sources.hpp"
#include "dvector.hpp"
#include "data.hpp"
#include "hugmodel.hpp"
#include "hugmh.hpp"
#include "hull.hpp"
#include "intruesourcesinteraction.hpp"
#include "straussinteraction.hpp"
#include "inconvexhullinteraction.hpp"
#include "abcshadow_hug.hpp"

#include "mydefines.hpp" // numberOfComponent

#include <cstdlib>
#include <cstdio>
#include <random>
#include <iostream>
#include <string>
#include <fstream>

int main(int argc, char *argv[])
{
  // Alea seeds
  int seed_a, seed_b, seed_c;
  // dimension window
  double K_x, K_y;
  // Observed statistics
  double nt_observed_statistics;
  double sr_observed_statistics;
  double ar_observed_statistics;
  double ne_observed_statistics;
  // Area
  double min_theta_e;
  double max_theta_e;
  double delta_theta_e;
  double theta_e_x_initial;
  // NotExplained
  double min_theta_a;
  double max_theta_a;
  double delta_theta_a;
  double theta_a_x_initial;
  double percent;
  // Poisson
  double min_theta_s;
  double max_theta_s;
  double delta_theta_s;
  double theta_s_x_initial;
  // Strauss
  double min_theta4;
  double max_theta4;
  double delta_theta4;
  double theta4_x_initial;
  double radius_fixed_Strauss;
  double radius_random_min_Strauss;
  double radius_random_max_Strauss;
  int alea_type_Strauss;
  double alea_parameter_Strauss;

  int dim;

  int number_of_loop_theta; // number of generations of realisation in ABC
  int nbiter_theta;         // number of theta updates between 2 realisations in ABC
  double probability_birth = 0.2;
  double probability_death = 0.2, probaChangeIn = 0.2, probaChangeOut = 0.05;
  double rayonChange = 0.3;
  int mh_time; // nb iterations in MH in the generation of realisation

  // sources, hull and other details;
  char selected_planes[1000] = "(1,2);(1,3);(2,3)",
       radius_planes[1000] = "0.1;0.1;0.1",
       data[200] = "../DATA/data.txt",
       truesources[200] = "../DATA/true_sources.txt",
       sources[200] = "../RESULTS/sources.txt",
       sources_final[200] = "../RESULTS/sources_final.txt",
       sources_detail[200] = "../RESULTS/sources_detail.txt",
       sources_detail_final[200] = "../RESULTS/sources_detail_final.txt",
       theta_detail[200] = "../RESULTS/theta.txt",
       parameters[200] = "../RESULTS/plan.txt";

  double observed_statistics[numberOfComponent];
  double list_bound_min_log[numberOfComponent];
  double list_bound_max_log[numberOfComponent];
  double list_delta_log[numberOfComponent];
  double list_initial_log[numberOfComponent];

  FILE *parameters_file;

  char inutil[200];
  int lval;

  clock_t timer1, timer2, timer3;

  timer1 = clock();

  // reading the parameters
  if (argc != 2)
  {
    cout << "ATTENTION : sim_abc.exe <parameters file> !!!\n";
    exit(1);
  }
  else
  {
    parameters_file = fopen(argv[1], "rt");
    if (parameters_file == NULL)
    {
      cout << "ERROR: Opening file " << argv[1] << " !\n";
      exit(1);
    }
    //===============================================
    // Seeds generator if 0 0 0 total random
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%d", &seed_a);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%d", &seed_b);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%d", &seed_c);

    //===============================================
    // Observed_statistics
    //   NotExplained ne_observed_statistics
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &ne_observed_statistics);
    //   NotInDisk nt_observed_statistics
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &nt_observed_statistics);
    //   Strauss sr_observed_statistics
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &sr_observed_statistics);
    //===============================================
    // ABC parameters
    //   time_mh_theta in ABCshadow main loop
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%d", &number_of_loop_theta);
    //   nbiter_theta in generation_new_parameter_candidates
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%d", &nbiter_theta);
    //===============================================
    // MH parameters
    //   probability_birth MH
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &probability_birth);
    //   probability_death MH
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &probability_death);

    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &probaChangeIn);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &probaChangeOut);

    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &rayonChange);
    //   mh_time_:_for_generating_an_x_sample
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%d", &mh_time);
    //===============================================
    // Research domain
    // Area  ===========
    //   NotExplained Component min_theta_e
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &min_theta_e);
    //   NotExplained Component max_theta_e
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &max_theta_e);
    //   NotExplained Component delta_theta_e
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &delta_theta_e);
    //   NotExplained Component theta_e_x_initial
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &theta_e_x_initial);
    // NotInDisk  ==========
    //   NotInDisk Component min_theta_a
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &min_theta_a);
    //   NotInDisk Component max_theta_a
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &max_theta_a);
    //   NotInDisk Component delta_theta_a
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &delta_theta_a);
    //   NotInDisk Component theta_a_x_initial
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &theta_a_x_initial);
    //   NotInDisk Component radius_planes
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", radius_planes);
    // STRAUSS  ===========
    //   Stauss Component min_theta_s
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &min_theta_s);
    //   Stauss Component max_theta_s
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &max_theta_s);
    //   Stauss Component delta_theta_s
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &delta_theta_s);
    //   Stauss Component theta_s_x_initial
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &theta_s_x_initial);
    //   = Fixed
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &radius_fixed_Strauss);

    //===============================================
    // RESULTS
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%d", &dim);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", selected_planes);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%lf", &percent);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", data);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", truesources);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", sources);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", sources_final);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", sources_detail);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", sources_detail_final);
    lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
    lval = fscanf(parameters_file, "%s", parameters);
    fclose(parameters_file);
  }

  // To facilitate the reproducibility of results seed_a seed_b seed_c can be fixed
  // else if (seed_a == 0) && (seed_b == 0) && (seed_c == 0) seed_$ are random
  if ((seed_a == 0) && (seed_b == 0) && (seed_c == 0))
  {
    default_random_engine defEngine(time(0));
    uniform_int_distribution<int> intDistro(1, 30000);
    seed_a = intDistro(defEngine);
    seed_b = intDistro(defEngine);
    seed_c = intDistro(defEngine);
  }
  whran.setSeed(seed_a, seed_b, seed_c);

  ////###############################################################################################################

  // ###################################### Hug complet

  //===============================================
  // Data conditioning
  //   Observed_statistics
  observed_statistics[0] = ne_observed_statistics;
  observed_statistics[1] = nt_observed_statistics;
  observed_statistics[2] = sr_observed_statistics;
  //   Min bounds of research domain
  list_bound_min_log[0] = -max_theta_e;
  list_bound_min_log[1] = -max_theta_a;
  list_bound_min_log[2] = -max_theta_s;
  //   Max bounds of research domain
  list_bound_max_log[0] = -min_theta_e;
  list_bound_max_log[1] = -min_theta_a;
  list_bound_max_log[2] = -min_theta_s;
  //   Delta max in log
  list_delta_log[0] = delta_theta_e;
  list_delta_log[1] = delta_theta_a;
  list_delta_log[2] = delta_theta_s;
  //   Initial positions of log in research domain
  list_initial_log[0] = -theta_e_x_initial;
  list_initial_log[1] = -theta_a_x_initial;
  list_initial_log[2] = -theta_s_x_initial;

  //===============================================
  // Component definitions
  GenericInteraction *list_interaction[numberOfComponent];

  // InConvexHull
  InConvexHullInteraction inconvexhullComponent(-theta_e_x_initial, percent, 0);
  list_interaction[0] = &inconvexhullComponent;

  // POISSON
  InTrueSourcesInteraction intruesourcesComponent(-theta_a_x_initial, 1);
  list_interaction[1] = &intruesourcesComponent;

  // STRAUSS
  StraussInteraction2 straussComponent(-theta_s_x_initial, radius_fixed_Strauss, 2);
  list_interaction[2] = &straussComponent;
  //===============================================

  ////###############################################################################################################

  // Model definition

  // Data
  Data thedata;

  thedata.loadDatasFromFile(data, dim, radius_planes);
  // thedata.display();

  timer2 = clock();

  // Sources (Pattern)
  Sources the_sources(&thedata);
  the_sources.saveHeaders(sources);
  the_sources.saveHeaders(sources_final);

  // the_sources.display();
  vector<vector<int>> planes = splitPlane(selected_planes, dim);
  thedata.loadSourcesFromFile(truesources, planes);

  // Model for fixed case
  //   Don't forget Event* HugModel :: newEvent() const redefinition
  HugModel the_model(list_interaction, &thedata, &the_sources, planes);
  the_model.saveStatisticsHeaders(sources_detail);
  the_model.saveStatisticsHeaders(sources_detail_final);
  // cout<<"model"<<endl;
  double temperature = 1.0, temperaturemin = 1.0, coolingcoeff = 1.0;
  int numberiteration = 1, savetime = 1000;

  HugMH algo_sampler(&thedata, &the_model, &the_sources, numberiteration,
                     mh_time,
                     probability_birth, probability_death, probaChangeIn, probaChangeOut,
                     rayonChange, temperature, temperaturemin, coolingcoeff,
                     savetime, sources, sources_detail, sources_final, sources_detail_final);

  //===============================================
  // ABCshadow definition
  ABCshadowHug shadowHugModel(list_bound_min_log,
                              list_bound_max_log,
                              list_delta_log,
                              list_initial_log,
                              observed_statistics,
                              number_of_loop_theta,
                              nbiter_theta,
                              mh_time,
                              probability_birth,
                              probability_death, probaChangeIn, probaChangeOut, rayonChange,
                              &thedata, &the_model, &the_sources, savetime,
                              sources, sources_detail, sources_final, sources_detail_final, parameters);

  //===============================================
  // Computation and results
  shadowHugModel.compute_parameters();
  // shadowStandardModel.report();
  //===============================================
}
