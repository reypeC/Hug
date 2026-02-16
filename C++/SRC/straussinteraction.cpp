#include "straussinteraction.hpp"

// ###################################################################################################################################################
// ###################################################################################################################################################

// 1

// ###################################################################################################################################################
// ###################################################################################################################################################
// ######################################################## Strauss
// ###################################################################################################################################################
// ###################################################################################################################################################

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StraussInteraction1 ::StraussInteraction1(double parameter_model_parameter,
                                          double radius_parameter,
                                          int rank_parameter) : GenericInteraction(parameter_model_parameter,
                                                                                   rank_parameter)
{

  width = 1.0;
  height = 1.0;
  area = 1.0;
  r = radius_parameter * radius_parameter;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StraussInteraction1 ::StraussInteraction1(double parameter_model_parameter,
                                          double radius_parameter,
                                          double width_parameter, double height_parameter,
                                          int rank_parameter) : GenericInteraction(parameter_model_parameter,
                                                                                   rank_parameter)
{

  width = width_parameter;
  height = height_parameter;
  area = width_parameter * height_parameter;
  r = radius_parameter * radius_parameter;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double StraussInteraction1 ::getCurrentStatistic(int indice_plan)
{

  return current_statistic[indice_plan];
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction1 ::applyEvolution()
{

  for (std::size_t indice_plan = 0; indice_plan < current_statistic.size(); indice_plan++)
  {
    current_statistic[indice_plan] = next_statistic[indice_plan];
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction1 ::validationOfBirthEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction1 ::noValidationOfBirthEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction1 ::validationOfDeathEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction1 ::noValidationOfDeathEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction1 ::validationOfChangeEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction1 ::noValidationOfChangeEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// void StraussInteraction1 :: decoreNewEvent( DecoratedEvent* a_decorated_event,int i )
//{
//
//  a_decorated_event->setDecoration(i,r);
//
//}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction1 ::computeStatistic(HugModel *the_hugmodel,
                                            vector<vector<Point>> *p_used_sources,
                                            vector<vector<Point>> *p_current_convex_hull,
                                            vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet is the last of vector
  int nb_sources = the_hugmodel->p_the_sources->getNumberSources();
  int number_plane = the_hugmodel->planes.size();
  current_statistic.resize(number_plane);
  next_statistic.resize(number_plane);
  // toto proposed_sources.pop_back();

  double nb_couples_in_interaction = 0.0;

  for (std::size_t i = 0; i < (nb_sources - 1); i++)
  {
    for (std::size_t j = (i + 1); j < nb_sources; j++)
    {
      double dist = 0.0;
      for (std::size_t elt = 0; elt < the_hugmodel->number_elements; elt++)
      {
        dist += (the_hugmodel->p_the_sources->getVal(i, elt) - the_hugmodel->p_the_sources->getVal(j, elt)) * (the_hugmodel->p_the_sources->getVal(i, elt) - the_hugmodel->p_the_sources->getVal(j, elt));
      }

      if (dist < r)
      {
        nb_couples_in_interaction = nb_couples_in_interaction + 1.0;
      }
    }
  }

  for (std::size_t indice_plan = 0; indice_plan < the_hugmodel->planes.size(); indice_plan++)
  {
    current_statistic[indice_plan] = nb_couples_in_interaction;
    next_statistic[indice_plan] = nb_couples_in_interaction;
  }
  // NULL because the computation is made in conditional computation and validation
  // the result is in    current_statistic
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double StraussInteraction1 ::cond_intens_birth(vector<Coordinates> proposed_sources,
                                               HugModel *the_hugmodel,
                                               vector<vector<Point>> *p_used_sources,
                                               vector<vector<Point>> *p_current_convex_hull,
                                               vector<vector<Point>> *p_next_convex_hull,
                                               vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet is the last of vector
  int nb_sources = proposed_sources.size() - 1;
  Coordinates lastPoint = proposed_sources.back();
  proposed_sources.pop_back();

  double nb_couples_in_interaction = 0.0;

  for (std::size_t i = 0; i < nb_sources; i++)
  {
    double dist = 0.0;
    for (std::size_t elt = 0; elt < the_hugmodel->number_elements; elt++)
    {
      dist += (proposed_sources[i][elt] - lastPoint[elt]) * (proposed_sources[i][elt] - lastPoint[elt]);
    }

    if (dist < r)
    {
      nb_couples_in_interaction = nb_couples_in_interaction + 1.0;
    }
  }
  for (std::size_t indice_plan = 0; indice_plan < the_hugmodel->planes.size(); indice_plan++)
  {
    next_statistic[indice_plan] = current_statistic[indice_plan] + nb_couples_in_interaction;
  }

  return exp(parameter_model * nb_couples_in_interaction);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double StraussInteraction1 ::cond_intens_death(vector<Coordinates> proposed_sources,
                                               HugModel *the_hugmodel,
                                               vector<vector<Point>> *p_used_sources,
                                               vector<vector<Point>> *p_current_convex_hull,
                                               vector<vector<Point>> *p_next_convex_hull,
                                               vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet is the last of vector
  int nb_sources = proposed_sources.size() - 1;
  Coordinates lastPoint = proposed_sources.back();
  proposed_sources.pop_back();

  double nb_couples_in_interaction = 0.0;

  for (std::size_t i = 0; i < nb_sources; i++)
  {
    double dist = 0.0;
    for (std::size_t elt = 0; elt < the_hugmodel->number_elements; elt++)
    {
      dist += (proposed_sources[i][elt] - lastPoint[elt]) * (proposed_sources[i][elt] - lastPoint[elt]);
    }

    if (dist < r)
    {
      nb_couples_in_interaction = nb_couples_in_interaction - 1.0;
    }
  }
  for (std::size_t indice_plan = 0; indice_plan < the_hugmodel->planes.size(); indice_plan++)
  {
    next_statistic[indice_plan] = current_statistic[indice_plan] + nb_couples_in_interaction;
  }

  return exp(parameter_model * nb_couples_in_interaction);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double StraussInteraction1 ::cond_intens_change(vector<Coordinates> proposed_sources,
                                                HugModel *the_hugmodel,
                                                vector<vector<Point>> *p_used_sources,
                                                vector<vector<Point>> *p_current_convex_hull,
                                                vector<vector<Point>> *p_next_convex_hull,
                                                vector<vector<Point>> *elt_usefulDatas)
{

  // The deleted point is the last of vector, the modified point is the one before
  int nb_sources = proposed_sources.size() - 2;
  Coordinates lastPoint = proposed_sources.back();
  proposed_sources.pop_back();
  Coordinates modifiedPoint = proposed_sources.back();
  proposed_sources.pop_back();

  double nb_couples_in_interaction = 0.0;

  for (std::size_t i = 0; i < nb_sources; i++)
  {
    double distD = 0.0;
    double distM = 0.0;
    for (std::size_t elt = 0; elt < the_hugmodel->number_elements; elt++)
    {
      distD += (proposed_sources[i][elt] - lastPoint[elt]) * (proposed_sources[i][elt] - lastPoint[elt]);
      distM += (proposed_sources[i][elt] - modifiedPoint[elt]) * (proposed_sources[i][elt] - modifiedPoint[elt]);
    }
    if (distD < r)
    {
      nb_couples_in_interaction = nb_couples_in_interaction - 1.0;
    }
    if (distM < r)
    {
      nb_couples_in_interaction = nb_couples_in_interaction + 1.0;
    }
  }
  for (std::size_t indice_plan = 0; indice_plan < the_hugmodel->planes.size(); indice_plan++)
  {
    next_statistic[indice_plan] = current_statistic[indice_plan] + nb_couples_in_interaction;
  }

  return exp(parameter_model * nb_couples_in_interaction);
}

// ###################################################################################################################################################
// ###################################################################################################################################################

// 2

// ###################################################################################################################################################
// ###################################################################################################################################################
// ######################################################## Hardcore
// ###################################################################################################################################################
// ###################################################################################################################################################

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StraussInteraction2 ::StraussInteraction2(double parameter_model_parameter,
                                          double radius_parameter,
                                          int rank_parameter) : GenericInteraction(parameter_model_parameter,
                                                                                   rank_parameter)
{

  width = 1.0;
  height = 1.0;
  area = 1.0;
  r = radius_parameter * radius_parameter;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StraussInteraction2 ::StraussInteraction2(double parameter_model_parameter,
                                          double radius_parameter,
                                          double width_parameter, double height_parameter,
                                          int rank_parameter) : GenericInteraction(parameter_model_parameter,
                                                                                   rank_parameter)
{

  width = width_parameter;
  height = height_parameter;
  area = width_parameter * height_parameter;
  r = radius_parameter * radius_parameter;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double StraussInteraction2 ::getCurrentStatistic(int indice_plan)
{

  return current_statistic[indice_plan];
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction2 ::applyEvolution()
{

  for (std::size_t indice_plan = 0; indice_plan < current_statistic.size(); indice_plan++)
  {
    current_statistic[indice_plan] = next_statistic[indice_plan];
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction2 ::validationOfBirthEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction2 ::noValidationOfBirthEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction2 ::validationOfDeathEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction2 ::noValidationOfDeathEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction2 ::validationOfChangeEvolution()
{

  applyEvolution();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction2 ::noValidationOfChangeEvolution()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// void StraussInteraction2 :: decoreNewEvent( DecoratedEvent* a_decorated_event,int i )
//{
//
//  a_decorated_event->setDecoration(i,r);
//
//}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void StraussInteraction2 ::computeStatistic(HugModel *the_hugmodel,
                                            vector<vector<Point>> *p_used_sources,
                                            vector<vector<Point>> *p_current_convex_hull,
                                            vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet is the last of vector
  int nb_sources = the_hugmodel->p_the_sources->getNumberSources();
  int number_plane = the_hugmodel->planes.size();
  current_statistic.resize(number_plane);
  next_statistic.resize(number_plane);
  // toto proposed_sources.pop_back();

  double nb_couples_in_interaction = 0.0;

  for (std::size_t i = 0; i < (nb_sources - 1); i++)
  {
    for (std::size_t j = (i + 1); j < nb_sources; j++)
    {
      double dist = 0.0;
      for (std::size_t elt = 0; elt < the_hugmodel->number_elements; elt++)
      {
        dist += (the_hugmodel->p_the_sources->getVal(i, elt) - the_hugmodel->p_the_sources->getVal(j, elt)) * (the_hugmodel->p_the_sources->getVal(i, elt) - the_hugmodel->p_the_sources->getVal(j, elt));
      }

      if (dist < r)
      {
        nb_couples_in_interaction = nb_couples_in_interaction + 1.0;
      }
    }
  }

  for (std::size_t indice_plan = 0; indice_plan < the_hugmodel->planes.size(); indice_plan++)
  {
    current_statistic[indice_plan] = nb_couples_in_interaction;
    next_statistic[indice_plan] = nb_couples_in_interaction;
  }
  // NULL because the computation is made in conditional computation and validation
  // the result is in    current_statistic
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double StraussInteraction2 ::cond_intens_birth(vector<Coordinates> proposed_sources,
                                               HugModel *the_hugmodel,
                                               vector<vector<Point>> *p_used_sources,
                                               vector<vector<Point>> *p_current_convex_hull,
                                               vector<vector<Point>> *p_next_convex_hull,
                                               vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet is the last of vector
  int nb_sources = proposed_sources.size() - 1;
  Coordinates lastPoint = proposed_sources.back();
  proposed_sources.pop_back();

  double nb_couples_in_interaction = 0.0;

  for (std::size_t i = 0; i < nb_sources; i++)
  {
    double dist = 0.0;
    for (std::size_t elt = 0; elt < the_hugmodel->number_elements; elt++)
    {
      dist += (proposed_sources[i][elt] - lastPoint[elt]) * (proposed_sources[i][elt] - lastPoint[elt]);
    }

    if (dist < r)
    {
      return 0.0;
    }
  }

  return 1.0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double StraussInteraction2 ::cond_intens_death(vector<Coordinates> proposed_sources,
                                               HugModel *the_hugmodel,
                                               vector<vector<Point>> *p_used_sources,
                                               vector<vector<Point>> *p_current_convex_hull,
                                               vector<vector<Point>> *p_next_convex_hull,
                                               vector<vector<Point>> *elt_usefulDatas)
{
  // The point of interet is the last of vector
  int nb_sources = proposed_sources.size() - 1;
  Coordinates lastPoint = proposed_sources.back();
  proposed_sources.pop_back();

  double nb_couples_in_interaction = 0.0;

  for (std::size_t i = 0; i < nb_sources; i++)
  {
    double dist = 0.0;
    for (std::size_t elt = 0; elt < the_hugmodel->number_elements; elt++)
    {
      dist += (proposed_sources[i][elt] - lastPoint[elt]) * (proposed_sources[i][elt] - lastPoint[elt]);
    }

    if (dist < r)
    {
      return 0.0;
    }
  }

  return 1.0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double StraussInteraction2 ::cond_intens_change(vector<Coordinates> proposed_sources,
                                                HugModel *the_hugmodel,
                                                vector<vector<Point>> *p_used_sources,
                                                vector<vector<Point>> *p_current_convex_hull,
                                                vector<vector<Point>> *p_next_convex_hull,
                                                vector<vector<Point>> *elt_usefulDatas)
{

  // The deleted point is the last of vector, the modified point is the one before
  int nb_sources = proposed_sources.size() - 2;
  Coordinates lastPoint = proposed_sources.back();
  proposed_sources.pop_back();
  Coordinates modifiedPoint = proposed_sources.back();
  proposed_sources.pop_back();

  double nb_couples_in_interaction = 0.0;

  for (std::size_t i = 0; i < nb_sources; i++)
  {
    double distD = 0.0;
    double distM = 0.0;
    for (std::size_t elt = 0; elt < the_hugmodel->number_elements; elt++)
    {
      distD += (proposed_sources[i][elt] - lastPoint[elt]) * (proposed_sources[i][elt] - lastPoint[elt]);
      distM += (proposed_sources[i][elt] - modifiedPoint[elt]) * (proposed_sources[i][elt] - modifiedPoint[elt]);
    }
    if (distD < r)
    {
      nb_couples_in_interaction = nb_couples_in_interaction - 1.0;
    }
    if (distM < r)
    {

      return 0.0;
    }
  }

  return 1.0;
}

// ###################################################################################################################################################
// ###################################################################################################################################################
