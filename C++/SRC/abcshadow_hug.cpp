#include "abcshadow_hug.hpp"
#include "sources.hpp"
#include "hugmh.hpp"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ABCshadowHug ::ABCshadowHug(double list_bound_min_log_parameter[],                                                                                                                                        // numberOfComponent
                            double list_bound_max_log_parameter[],                                                                                                                                        // numberOfComponent
                            double list_delta_log_parameter[],                                                                                                                                            // numberOfComponent
                            double list_initial_log_parameter[],                                                                                                                                          // numberOfComponent
                            double observed_statistics_parameter[],                                                                                                                                       // numberOfComponent
                            int number_of_loop_parameter,                                                                                                                                                 // ABCshadow main loop
                            int number_of_step_parameter,                                                                                                                                                 // ABCshadow nb parameter updates by loop
                            int time_mh_parameter,                                                                                                                                                        // MH nb iteration in MH to generate a realization
                            double probability_birth_parameter, double probability_death_parameter, double probability_changein_parameter, double probability_changeout_parameter, double p_radiuschange, // MH birth, death
                            Data *p_datas, HugModel *p_themodel, Sources *p_sources,
                            int p_savingtime,
                            char *p_sourcesfile, char *p_statisticsfile, char *p_sourcesfilefinal, char *p_statisticsfilefinal, char *p_parameters)
{
  list_bound_min_log = list_bound_min_log_parameter;
  list_bound_max_log = list_bound_max_log_parameter;
  list_delta_log = list_delta_log_parameter;
  number_of_loop = number_of_loop_parameter;
  number_of_step = number_of_step_parameter;
  time_mh = time_mh_parameter;
  probability_birth = probability_birth_parameter;
  probability_death = probability_death_parameter;
  probability_changein = probability_changein_parameter;
  probability_changeout_parameter = probability_changeout_parameter;
  radiuschange = p_radiuschange;
  desired_model = p_themodel;
  desired_datas = p_datas;
  current_sources = p_sources;
  savingtime = p_savingtime;
  sourcesfile = p_sourcesfile;
  statisticsfile = p_statisticsfile;
  sourcesfilefinal = p_sourcesfilefinal;
  statisticsfilefinal = p_statisticsfilefinal;
  parametersfile = p_parameters;

  // Initialization
  int i;
  for (i = 0; i < numberOfComponent; i++)
  {
    list_bound_min_log_parameter[i] = list_bound_min_log_parameter[i];
    list_bound_max_log_parameter[i] = list_bound_max_log_parameter[i];
    logtheta[i] = list_initial_log_parameter[i];
    logtheta_candidate[i] = list_initial_log_parameter[i];
    observed_statistics[i] = observed_statistics_parameter[i];
    current_pattern_statistics[i] = observed_statistics_parameter[i];
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ABCshadowHug ::compute_parameters()
{
  int step;
  double alea;
  double acceptance_threshold;

  int counterloop = number_of_loop;
  int counterstep = number_of_step;

  while (0 < counterloop)
  {
    cout << "counterloop = " << counterloop << "\n";
    generation_realization();
    int counterstep = number_of_step;
    while (0 < counterstep)
    {
      generation_new_parameter_candidates();
      alea = whran.uniform();
      acceptance_threshold = acceptance_threshold_computation();
      if (alea < acceptance_threshold)
      {
        update_parameters();
        counterstep = counterstep - 1;
      }
    }
    save_current_theta();
    counterloop = counterloop - 1;
  }
  save_in_file_last_realization();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ABCshadowHug ::generation_realization()
{
  int i;

  double temperature = 1.0, temperaturemin = 1.0, coolingcoeff = 1.0;
  int numberiteration = 1, savetime = 1000;

  HugMH algo_sampler(desired_datas, desired_model, current_sources, numberiteration, time_mh, probability_birth, probability_death,probability_changein,probability_changeout, radiuschange,
                     temperature, temperaturemin, coolingcoeff,
                     savetime, sourcesfile, statisticsfile, sourcesfilefinal, statisticsfilefinal);
  algo_sampler.sim();

  desired_model = algo_sampler.theModel;
  // cout<<"number sources = "<<desired_model->p_the_sources->getNumberSources()<<endl;
  current_sources = algo_sampler.theSources;

  for (i = 0; i < numberOfComponent; i++)
  {
    double statistic = 0.0;
    for (int indice_plan = 0; indice_plan < desired_model->number_planes; indice_plan++)
    {
      statistic += desired_model->list_interaction_model[i]->getCurrentStatistic(indice_plan);
    }

    current_pattern_statistics[i] = statistic / desired_model->number_planes;
  }
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ABCshadowHug ::generation_new_parameter_candidates()
{
  int pos; // random component of theta for evolution
  double logtheta_tmp;
  int i;

  for (i = 0; i < numberOfComponent; i++)
  {
    logtheta_candidate[i] = logtheta[i];
    // cout<<"candidate theta init = "<<logtheta_candidate[i]<<endl;
  }

  for (i = 0; i < number_of_step; i++)
  {
    pos = whran.discrete_uniform(0, numberOfComponent - 1);
    logtheta_tmp = logtheta_candidate[pos] + (list_delta_log[pos] * (whran.uniform() - 0.5));
    if ((list_bound_min_log[pos] <= logtheta_tmp) &&
        (logtheta_tmp <= list_bound_max_log[pos]))
    {
      logtheta_candidate[pos] = logtheta_tmp;
      // cout<<"candidate theta = "<<logtheta_candidate[i]<<endl;
    }
  }
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double ABCshadowHug ::acceptance_threshold_computation()
{
  double ratio = 1.0;
  int i;

  for (i = 0; i < numberOfComponent; i++)
  {

    ratio = ratio * exp((current_pattern_statistics[i] - observed_statistics[i]) *
                        (logtheta[i] - logtheta_candidate[i]));
  }

  return ratio;
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ABCshadowHug ::update_parameters()
{
  int i;

  for (i = 0; i < numberOfComponent; i++)
  {
    logtheta[i] = logtheta_candidate[i];
    desired_model->list_interaction_model[i]->setParameterInteraction(logtheta[i]);
    desired_model->list_interaction_model[i]->setExpParameterInteraction(exp(logtheta[i]));
    // cout<<"update theta = "<<logtheta[i]<<endl;
  }
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ABCshadowHug ::save_current_theta()
{
  int i;

  FILE *list_file;
  if ((list_file = fopen(parametersfile, "at")) == NULL)
  {
    cout << "ERROR : opening file in ABCshadowAreaIntFixed::save_current_theta()";
    cout << parametersfile << " \n";
  }
  else
  {
    for (i = 0; i < numberOfComponent; i++)
    {
      fprintf(list_file, "%f ", -logtheta[i]);
    }
    fprintf(list_file, "\n");
    fclose(list_file);
  }
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ABCshadowHug ::save_in_file_last_realization()
{
  current_sources->saveMatrixRaw(sourcesfile);
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
