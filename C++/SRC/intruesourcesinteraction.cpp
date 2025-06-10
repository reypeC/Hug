#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>

#include "intruesourcesinteraction.hpp"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InTrueSourcesInteraction ::InTrueSourcesInteraction(double parameter_model_parameter,
                                                    int rank_parameter) : GenericInteraction(parameter_model_parameter, rank_parameter)
{
  width = 1.0;
  height = 1.0;
  area = width * height;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InTrueSourcesInteraction ::InTrueSourcesInteraction(double parameter_model_parameter,
                                                    double width_parameter, double height_parameter,
                                                    int rank_parameter) : GenericInteraction(parameter_model_parameter, rank_parameter)
{
  width = width_parameter;
  height = height_parameter;
  area = width * height;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InTrueSourcesInteraction ::getCurrentStatistic(int indice_plan)
{
  return current_statistic[indice_plan];
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InTrueSourcesInteraction ::cond_intens_birth(vector<Coordinates> proposed_sources,
                                                    HugModel *the_hugmodel,
                                                    vector<vector<Point>> *p_used_sources,
                                                    vector<vector<Point>> *p_current_convex_hull,
                                                    vector<vector<Point>> *p_next_convex_hull,
                                                    vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet is the last of vector
  Coordinates lastPoint = proposed_sources.back();

  double statistic = 0.0;
  double dist;
  int number_plane = the_hugmodel->planes.size();

  for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
  {
    int elt_a_p = the_hugmodel->planes[indice_plan][0];
    int elt_b_p = the_hugmodel->planes[indice_plan][1];

    bool is_outside_all_disks = true;

    for (size_t indice_true_source = 0; indice_true_source < the_hugmodel->p_the_datas->trueSourcesNumber[indice_plan]; indice_true_source++)
    {
      dist = (lastPoint[elt_a_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x) * (lastPoint[elt_a_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x);
      dist += (lastPoint[elt_b_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y) * (lastPoint[elt_b_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y);

      if (dist < radius[indice_plan])
      {
        is_outside_all_disks = false;
        break;
      }
    }
    if (is_outside_all_disks)
    {
      statistic += 1.0;
      next_statistic[indice_plan] = current_statistic[indice_plan] + 1;
    }
    else
    {
      next_statistic[indice_plan] = current_statistic[indice_plan];
    }
  }

  return exp(parameter_model * (statistic / number_plane));
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InTrueSourcesInteraction ::cond_intens_death(vector<Coordinates> proposed_sources,
                                                    HugModel *the_hugmodel,
                                                    vector<vector<Point>> *p_used_sources,
                                                    vector<vector<Point>> *p_current_convex_hull,
                                                    vector<vector<Point>> *p_next_convex_hull,
                                                    vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet is the last of vector
  Coordinates lastPoint = proposed_sources.back();

  double statistic = 0.0;
  double dist;
  int number_plane = the_hugmodel->planes.size();

  for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
  {
    int elt_a_p = the_hugmodel->planes[indice_plan][0];
    int elt_b_p = the_hugmodel->planes[indice_plan][1];

    bool is_outside_all_disks = true;

    for (size_t indice_true_source = 0; indice_true_source < the_hugmodel->p_the_datas->trueSourcesNumber[indice_plan]; indice_true_source++)
    {
      dist = (lastPoint[elt_a_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x) * (lastPoint[elt_a_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x);
      dist += (lastPoint[elt_b_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y) * (lastPoint[elt_b_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y);

      if (dist < radius[indice_plan])
      {
        is_outside_all_disks = false;
        break;
      }
    }
    if (is_outside_all_disks)
    {
      statistic -= 1.0;
      next_statistic[indice_plan] = current_statistic[indice_plan] - 1;
    }
    else
    {
      next_statistic[indice_plan] = current_statistic[indice_plan];
    }
  }
  return exp(parameter_model * (statistic / number_plane));
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InTrueSourcesInteraction ::cond_intens_change(vector<Coordinates> proposed_sources,
                                                     HugModel *the_hugmodel,
                                                     vector<vector<Point>> *p_used_sources,
                                                     vector<vector<Point>> *p_current_convex_hull,
                                                     vector<vector<Point>> *p_next_convex_hull,
                                                     vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet are the two last of vector
  Coordinates lastPoint = proposed_sources.back();
  proposed_sources.pop_back();
  Coordinates modifiedPoint = proposed_sources.back();

  double statistic = 0.0;
  double distdel;
  double distmod;
  int number_plane = the_hugmodel->planes.size();

  for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
  {
    int elt_a_p = the_hugmodel->planes[indice_plan][0];
    int elt_b_p = the_hugmodel->planes[indice_plan][1];

    bool is_outside_all_disks_del = true;
    bool is_outside_all_disks_mod = true;

    for (size_t indice_true_source = 0; indice_true_source < the_hugmodel->p_the_datas->trueSourcesNumber[indice_plan]; indice_true_source++)
    {
      distdel = (lastPoint[elt_a_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x) * (lastPoint[elt_a_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x);
      distdel += (lastPoint[elt_b_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y) * (lastPoint[elt_b_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y);

      if (distdel < radius[indice_plan])
      {
        is_outside_all_disks_del = false;
      }

      distmod = (modifiedPoint[elt_a_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x) * (modifiedPoint[elt_a_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x);
      distmod += (modifiedPoint[elt_b_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y) * (modifiedPoint[elt_b_p] - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y);

      if (distmod < radius[indice_plan])
      {
        is_outside_all_disks_mod = false;
      }
    }

    if (is_outside_all_disks_del)
    {
      statistic -= 1.0;
      next_statistic[indice_plan] = current_statistic[indice_plan] - 1;
    }
    else
    {
      next_statistic[indice_plan] = current_statistic[indice_plan];
    }

    if (is_outside_all_disks_mod)
    {
      statistic += 1.0;
      next_statistic[indice_plan] = next_statistic[indice_plan] + 1;
    }
  }

  return exp(parameter_model * (statistic / number_plane));
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InTrueSourcesInteraction ::applyEvolution()
{

  for (std::size_t indice_plan = 0; indice_plan < current_statistic.size(); indice_plan++)
  {
    current_statistic[indice_plan] = next_statistic[indice_plan];
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InTrueSourcesInteraction ::validationOfBirthEvolution()
{
  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InTrueSourcesInteraction ::noValidationOfBirthEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InTrueSourcesInteraction ::validationOfDeathEvolution()
{
  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InTrueSourcesInteraction ::noValidationOfDeathEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InTrueSourcesInteraction ::validationOfChangeEvolution()
{
  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InTrueSourcesInteraction ::noValidationOfChangeEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InTrueSourcesInteraction ::beta_parameter_model_func(Coordinates *e)
{
  return parameter_model;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InTrueSourcesInteraction ::beta_exp_parameter_model_func(Coordinates *e)
{
  return exp_parameter_model;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InTrueSourcesInteraction ::computeStatistic(HugModel *the_hugmodel,
                                                 vector<vector<Point>> *p_used_sources,
                                                 vector<vector<Point>> *p_current_convex_hull,
                                                 vector<vector<Point>> *elt_usefulDatas) // ( Sources * pSources)
{

  int number_plane = the_hugmodel->planes.size();
  int nb_sources = the_hugmodel->p_the_sources->getNumberSources();
  current_statistic.resize(number_plane);
  next_statistic.resize(number_plane);

  radius.resize(number_plane);

  for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
  {
    int elt_a_p = the_hugmodel->planes[indice_plan][0];
    int elt_b_p = the_hugmodel->planes[indice_plan][1];
    double dist;
    int number_sources_outside_disk = 0;

    radius[indice_plan] = the_hugmodel->meanDistDatasSquare[elt_a_p][elt_b_p];

    for (size_t indice_sources_pattern = 0; indice_sources_pattern < nb_sources; indice_sources_pattern++)
    {
      bool is_outside_all_disks = true;

      for (size_t indice_true_source = 0; indice_true_source < the_hugmodel->p_the_datas->trueSourcesNumber[indice_plan]; indice_true_source++)
      {
        dist = (the_hugmodel->p_the_sources->getVal(indice_sources_pattern, elt_a_p) - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x) * (the_hugmodel->p_the_sources->getVal(indice_sources_pattern, elt_a_p) - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].x);
        dist += (the_hugmodel->p_the_sources->getVal(indice_sources_pattern, elt_b_p) - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y) * (the_hugmodel->p_the_sources->getVal(indice_sources_pattern, elt_b_p) - the_hugmodel->p_the_datas->trueSources[indice_plan][indice_true_source].y);

        if (dist < radius[indice_plan])
        {
          is_outside_all_disks = false;
          break;
        }
      }
      if (is_outside_all_disks)
      {
        number_sources_outside_disk++;
      }
    }
    current_statistic[indice_plan] = number_sources_outside_disk;
    next_statistic[indice_plan] = number_sources_outside_disk;
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


