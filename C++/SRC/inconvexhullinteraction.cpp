#include "inconvexhullinteraction.hpp"
#include "hull.hpp"

// ###################################################################################################################################################
// ###################################################################################################################################################
// ######################################################## Ratio données non expliquées
// ###################################################################################################################################################
// ###################################################################################################################################################

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InConvexHullInteraction ::InConvexHullInteraction(double parameter_model_parameter, double percentage, int rank_parameter)
    : GenericInteraction(parameter_model_parameter, rank_parameter)
{

  width = 1.0;
  height = 1.0;
  area = 1.0;
  seuil = percentage;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InConvexHullInteraction ::InConvexHullInteraction(double parameter_model_parameter,
                                                  double width_parameter, double height_parameter, double percentage,
                                                  int rank_parameter)
    : GenericInteraction(parameter_model_parameter,
                         rank_parameter)
{

  width = width_parameter;
  height = height_parameter;
  area = width_parameter * height_parameter;
  seuil = percentage;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InConvexHullInteraction ::getCurrentStatistic(int indice_plan)
{

  return current_statistic[indice_plan];
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InConvexHullInteraction ::applyEvolution()
{

  for (std::size_t indice_plan = 0; indice_plan < current_statistic.size(); indice_plan++)
  {
    current_statistic[indice_plan] = next_statistic[indice_plan];
    for (std::size_t number_points = 0; number_points < is_explained[indice_plan].size(); number_points++)
    {
      is_explained[indice_plan][number_points] = is_explained_next[indice_plan][number_points];
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InConvexHullInteraction ::validationOfBirthEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InConvexHullInteraction ::noValidationOfBirthEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InConvexHullInteraction ::validationOfDeathEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InConvexHullInteraction ::noValidationOfDeathEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InConvexHullInteraction ::validationOfChangeEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InConvexHullInteraction ::noValidationOfChangeEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// void InConvexHullInteraction :: decoreNewEvent( DecoratedEvent* a_decorated_event,int i )
//{
//
//  a_decorated_event->setDecoration(i,r);
//
//}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InConvexHullInteraction ::cond_intens_birth(vector<Coordinates> proposed_sources,
                                                   HugModel *the_hugmodel,
                                                   vector<vector<Point>> *p_used_sources,
                                                   vector<vector<Point>> *p_current_convex_hull,
                                                   vector<vector<Point>> *p_next_convex_hull,
                                                   vector<vector<Point>> *elt_usefulDatas)
{
  int nb_sources = proposed_sources.size();
  double statistic = 0.0;
  int number_plane = current_statistic.size();

  if (nb_sources > 2)
  {
    for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
    {
      int nb_elt_usefulDatas = elt_usefulDatas->at(indice_plan).size();
      double nb_in = 0.0;

      // computation of number of sources in p_next_convex_hull

      for (std::size_t i = 0; i < nb_elt_usefulDatas; i++)
      {
        if (is_explained[indice_plan][i] == 0)
        {
          Point aPoint = elt_usefulDatas->at(indice_plan)[i];
          if (pointInConvexHull(p_next_convex_hull->at(indice_plan), aPoint))
          {
            nb_in = nb_in + 1.0;
            is_explained_next[indice_plan][i] = 1;
          }
          else
          {
            is_explained_next[indice_plan][i] = 0;
          }
        }
        else
        {
          nb_in = nb_in + 1.0;
          is_explained_next[indice_plan][i] = 1;
        }
      }

      if (nb_in <= seuil_explained[indice_plan])
      {
        return 0.0;
      }
      next_statistic[indice_plan] = (1.0 - (nb_in) / (nb_elt_usefulDatas))*100.0;
      statistic += next_statistic[indice_plan] - current_statistic[indice_plan];
    }
    return exp(parameter_model * (statistic / number_plane));
  }
  else
  {
    for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
    {
      next_statistic[indice_plan] = 100.0;
      for (std::size_t i = 0; i < is_explained[indice_plan].size(); i++)
      {
        is_explained_next[indice_plan][i] = 0;
      }
    }
  }

  return 0.0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InConvexHullInteraction ::cond_intens_death(vector<Coordinates> proposed_sources,
                                                   HugModel *the_hugmodel,
                                                   vector<vector<Point>> *p_used_sources,
                                                   vector<vector<Point>> *p_current_convex_hull,
                                                   vector<vector<Point>> *p_next_convex_hull,
                                                   vector<vector<Point>> *elt_usefulDatas)
{
  int nb_sources = proposed_sources.size() - 1;
  double statistic = 0.0;
  int number_plane = current_statistic.size();

  if (nb_sources > 2)
  {
    for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
    {
      int nb_elt_usefulDatas = elt_usefulDatas->at(indice_plan).size();
      double nb_in = 0.0;

      // computation of number of sources in p_next_convex_hull

      for (std::size_t i = 0; i < nb_elt_usefulDatas; i++)
      {
        if (is_explained[indice_plan][i] == 1)
        {
          Point aPoint = elt_usefulDatas->at(indice_plan)[i];
          if (pointInConvexHull(p_next_convex_hull->at(indice_plan), aPoint))
          {
            nb_in = nb_in + 1.0;
            is_explained_next[indice_plan][i] = 1;
          }
          else
          {
            is_explained_next[indice_plan][i] = 0;
          }
        }
        else
        {
          is_explained_next[indice_plan][i] = 0;
        }
      }

      if (nb_in <= seuil_explained[indice_plan])
      {
        return 0.0;
      }
      next_statistic[indice_plan] = (1.0 - (nb_in) / (nb_elt_usefulDatas))*100.0;
      statistic += next_statistic[indice_plan] - current_statistic[indice_plan];
    }
    return exp(parameter_model * (statistic / number_plane));
  }
  else
  {
    for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
    {
      next_statistic[indice_plan] = 100.0;
      for (std::size_t i = 0; i < is_explained[indice_plan].size(); i++)
      {
        is_explained_next[indice_plan][i] = 0;
      }
    }
  }

  return 0.0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double InConvexHullInteraction ::cond_intens_change(vector<Coordinates> proposed_sources,
                                                    HugModel *the_hugmodel,
                                                    vector<vector<Point>> *p_used_sources,
                                                    vector<vector<Point>> *p_current_convex_hull,
                                                    vector<vector<Point>> *p_next_convex_hull,
                                                    vector<vector<Point>> *elt_usefulDatas)
{
  int nb_sources = proposed_sources.size() - 1;
  double statistic = 0.0;
  int number_plane = current_statistic.size();

  if (nb_sources > 2)
  {
    for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
    {
      int elt_a_p = the_hugmodel->planes[indice_plan][0];
      int elt_b_p = the_hugmodel->planes[indice_plan][1];

      if ((elt_a_p == the_hugmodel->elt_1) || (elt_a_p == the_hugmodel->elt_2) || (elt_b_p == the_hugmodel->elt_1) || (elt_b_p == the_hugmodel->elt_2))
      {
        int nb_elt_usefulDatas = elt_usefulDatas->at(indice_plan).size();
        double nb_in = 0.0;

        // computation of number of sources in p_next_convex_hull

        for (std::size_t i = 0; i < nb_elt_usefulDatas; i++)
        {
          Point aPoint = elt_usefulDatas->at(indice_plan)[i];
          if (pointInConvexHull(p_next_convex_hull->at(indice_plan), aPoint))
          {
            nb_in = nb_in + 1.0;
            is_explained_next[indice_plan][i] = 1;
          }
          else
          {
            is_explained_next[indice_plan][i] = 0;
          }
        }

        if (nb_in <= seuil_explained[indice_plan])
        {
          return 0.0;
        }
        next_statistic[indice_plan] = (1.0 - (nb_in) / (nb_elt_usefulDatas))*100.0;
        statistic += next_statistic[indice_plan] - current_statistic[indice_plan];
      }
      else
      {
        for (std::size_t i = 0; i < is_explained[indice_plan].size(); i++)
        {
          is_explained_next[indice_plan][i] = is_explained[indice_plan][i];
        }
        next_statistic[indice_plan] = current_statistic[indice_plan];
      }
    }
    return exp(parameter_model * (statistic / number_plane));
  }
  else
  {
    for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
    {
      next_statistic[indice_plan] = 100.0;
      for (std::size_t i = 0; i < is_explained[indice_plan].size(); i++)
      {
        is_explained_next[indice_plan][i] = 0;
      }
    }
  }
  return 0.0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InConvexHullInteraction ::computeStatistic(HugModel *the_hugmodel,
                                                vector<vector<Point>> *p_used_sources,
                                                vector<vector<Point>> *p_current_convex_hull,
                                                vector<vector<Point>> *elt_usefulDatas)
{
  int nb_sources = the_hugmodel->p_the_sources->getNumberSources();
  vector<int> is_explained_temporary;
  int number_plane = the_hugmodel->planes.size();
  current_statistic.resize(number_plane);
  next_statistic.resize(number_plane);

  seuil_explained.resize(number_plane);
  is_explained.clear();

  if (nb_sources > 2)
  {
    double nb_in;
    // computation of number of sources in p_next_convex_hull
    for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
    {
      nb_in = 0.0;

      int nb_elt_usefulDatas = elt_usefulDatas->at(indice_plan).size();
      is_explained_temporary.resize(nb_elt_usefulDatas);

      for (std::size_t i = 0; i < nb_elt_usefulDatas; i++)
      {
        Point aPoint = elt_usefulDatas->at(indice_plan)[i];
        if (pointInConvexHull(p_current_convex_hull->at(indice_plan), aPoint))
        {
          nb_in = nb_in + 1.0;
          is_explained_temporary[i] = 1;
        }
        else
        {
          is_explained_temporary[i] = 0;
        }
      }
      is_explained.push_back(is_explained_temporary);
      is_explained_next.push_back(is_explained_temporary);

      current_statistic[indice_plan] = (1.0 - (nb_in) / (nb_elt_usefulDatas))*100.0;
      next_statistic[indice_plan] = (1.0 - (nb_in) / (nb_elt_usefulDatas))*100.0;

      seuil_explained[indice_plan] = floor((1.0 - seuil) * nb_elt_usefulDatas);
    }
  }
  else
  {
    for (std::size_t indice_plan = 0; indice_plan < number_plane; indice_plan++)
    {
      int nb_elt_usefulDatas = elt_usefulDatas->at(indice_plan).size();
      is_explained_temporary.resize(nb_elt_usefulDatas);

      for (std::size_t i = 0; i < nb_elt_usefulDatas; i++)
      {
        is_explained_temporary[i] = 0;
      }
      is_explained.push_back(is_explained_temporary);
      is_explained_next.push_back(is_explained_temporary);

      current_statistic[indice_plan] = 100.0;
      next_statistic[indice_plan] = 100.0;

      seuil_explained[indice_plan] = floor((1.0 - seuil) * nb_elt_usefulDatas);
    }
  }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// ###################################################################################################################################################
// ###################################################################################################################################################
