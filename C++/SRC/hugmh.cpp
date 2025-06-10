#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "hugmh.hpp"
#include "random.hpp"
#include "hull.hpp"

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HugMH ::HugMH(Data *p_datas, HugModel *p_themodel, Sources *p_sources, int p_numberiteration,
              double p_timeMH,
              double p_pbirth, double p_pdeath, double p_pchange_in, double p_pchange_out, double p_radiuschange,
              double p_temperature, double p_temperaturemin, double p_coolingcoeff, int p_savingtime,
              char *p_sourcesfile, char *p_statisticsfile, char *p_sourcesfilefinal, char *p_statisticsfilefinal)
{
  theData = p_datas;
  theModel = p_themodel;
  theSources = p_sources;
  timeMH = p_timeMH;
  numberiteration = p_numberiteration;
  temperature = p_temperature;
  temperaturemin = p_temperaturemin;
  coolingcoeff = p_coolingcoeff;
  pbirth = p_pbirth;
  pdeath = p_pdeath;
  pchangein = p_pchange_in;
  pchangeout = p_pchange_out;
  radiuschange = p_radiuschange;
  t = 0;
  number_elements = p_datas->number_elements;
  number_vectors = p_datas->number_vectors;
  savingtime = p_savingtime;
  sourcesfile = p_sourcesfile;
  statisticsfile = p_statisticsfile;
  sourcesfilefinal = p_sourcesfilefinal;
  statisticsfilefinal = p_statisticsfilefinal;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugMH ::birth()
{
  // Generation of a new source
  Coordinates aNewSource = theModel->newEvent();
  proposed_sources = theModel->p_the_sources->sourcesarray;
  proposed_sources.push_back(aNewSource);
  int current_nb_of_sources = theModel->p_the_sources->getNumberSources();
  int number_planes = theModel->number_planes;

  used_sources.clear();
  next_convex_hull.clear();
  vector<Point> temporary_used_sources;

  for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
  {
    int elt_a = theModel->planes[indice_plan][0];
    int elt_b = theModel->planes[indice_plan][1];
    temporary_used_sources.clear();

    for (int i = 0; i < current_nb_of_sources; i++)
    {
      Point aNewPoint(theModel->p_the_sources->getVal(i, elt_a),
                      theModel->p_the_sources->getVal(i, elt_b));
      temporary_used_sources.push_back(aNewPoint);
    }
    Point aNewPoint(aNewSource[elt_a], aNewSource[elt_b]);
    temporary_used_sources.push_back(aNewPoint);
    next_convex_hull.push_back(convex_hull(temporary_used_sources));
    used_sources.push_back(temporary_used_sources);
  }
  double w = whran.uniform();

  if ((w) <
      ((pow(theModel->cond_intens_birth(proposed_sources,
                                        theModel,
                                        &used_sources,
                                        &current_convex_hull,
                                        &next_convex_hull,
                                        &theModel->usefulDatas),
            1.0 / temperature) *
        theModel->area * pdeath) /
       ((theModel->p_the_sources->getNumberSources() + 1.0) * pbirth)))
  {
    theModel->p_the_sources->addSource(&aNewSource);
    theModel->validationOfBirthEvolution();
    proposed_sources.clear();
    current_convex_hull.clear();
    for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
    {
      temporary_used_sources.clear();
      for (int indice = 0; indice < next_convex_hull[indice_plan].size(); indice++)
      {
        temporary_used_sources.push_back(next_convex_hull[indice_plan][indice]);
      }
      current_convex_hull.push_back(temporary_used_sources);
    }
    temporary_used_sources.clear();
    next_convex_hull.clear();
  }
  else
  {
    theModel->noValidationOfBirthEvolution();
    proposed_sources.clear();
    next_convex_hull.clear();
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugMH ::death()
{
  int current_nb_of_sources = theModel->p_the_sources->getNumberSources();

  if (current_nb_of_sources > 0)
  {
    int number_planes = theModel->number_planes;
    proposed_sources = theModel->p_the_sources->sourcesarray;

    used_sources.clear();
    next_convex_hull.clear();

    // Selection of source to be perhaps deleted
    int index_selected_source = whran.discrete_uniform(0, current_nb_of_sources - 1);
    vector<Point> temporary_used_sources;
    for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
    {
      int elt_a = theModel->planes[indice_plan][0];
      int elt_b = theModel->planes[indice_plan][1];
      temporary_used_sources.clear();

      Point theDeletedPoint(-1.0, -1.0);
      // Construction couples sources without index_selected_source
      for (int i = 0; i < current_nb_of_sources; i++)
      {
        if (i != index_selected_source)
        {
          Point aNewPoint(theModel->p_the_sources->getVal(i, elt_a),
                          theModel->p_the_sources->getVal(i, elt_b));
          temporary_used_sources.push_back(aNewPoint);
        }
        else
        {
          theDeletedPoint.x = theModel->p_the_sources->getVal(i, elt_a);
          theDeletedPoint.y = theModel->p_the_sources->getVal(i, elt_b);
        }
      }
      next_convex_hull.push_back(convex_hull(temporary_used_sources));
      temporary_used_sources.push_back(theDeletedPoint);
      used_sources.push_back(temporary_used_sources);
    }

    double w = whran.uniform();

    Coordinates deleted_source;
    for (int i = 0; i < number_elements; i++)
    {
      deleted_source.push_back(proposed_sources[index_selected_source][i]);
      proposed_sources[index_selected_source][i] = proposed_sources[current_nb_of_sources - 1][i];
    }
    proposed_sources.pop_back();
    proposed_sources.push_back(deleted_source);

    if ((w) < (theModel->p_the_sources->getNumberSources() * pbirth *
               pow(theModel->cond_intens_death(proposed_sources,
                                               theModel,
                                               &used_sources,
                                               &current_convex_hull,
                                               &next_convex_hull,
                                               &theModel->usefulDatas),
                   1.0 / temperature) /
               (pdeath * theModel->area)))
    {
      theModel->p_the_sources->removeSource(index_selected_source);
      theModel->validationOfDeathEvolution();
      current_convex_hull.clear();
      proposed_sources.clear();
      for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
      {
        temporary_used_sources.clear();
        for (int indice = 0; indice < next_convex_hull[indice_plan].size(); indice++)
        {
          temporary_used_sources.push_back(next_convex_hull[indice_plan][indice]);
        }
        current_convex_hull.push_back(temporary_used_sources);
      }
      temporary_used_sources.clear();
      proposed_sources.clear();
      next_convex_hull.clear();
    }
    else
    {
      theModel->noValidationOfDeathEvolution();
      proposed_sources.clear();
      next_convex_hull.clear();
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugMH ::change()
{

  int current_nb_of_sources = theModel->p_the_sources->getNumberSources();
  proposed_sources = theModel->p_the_sources->sourcesarray;

  if (current_nb_of_sources > 0)
  {
    int number_planes = theModel->number_planes;
    used_sources.clear();
    next_convex_hull.clear();

    // int current_nb_of_sources=theModel->p_the_sources->getNumberSources();

    // Selection of source to be perhaps modified
    int index_selected_source = whran.discrete_uniform(0, current_nb_of_sources - 1);
    int index_selected_plane = whran.discrete_uniform(0, number_planes - 1);
    int elt_1 = theModel->planes[index_selected_plane][0];
    int elt_2 = theModel->planes[index_selected_plane][1];
    theModel->setElement(elt_1, elt_2);

    Coordinates modifiedSource;
    Coordinates deletedSource;
    double radius = radiuschange * whran.uniform();
    double angle = (2.0 * M_PI) * whran.uniform();

    for (int i = 0; i < number_elements; i++)
    {
      double newValue = theModel->p_the_sources->getVal(index_selected_source, i);
      deletedSource.push_back(newValue);
      if (i == elt_1)
      {
        double modificationSource = radius * cos(angle);
        newValue += modificationSource;
        if ((newValue > 1) || (newValue < 0))
        {
          newValue -= 2 * modificationSource;
        }
      }
      else if (i == elt_2)
      {
        double modificationSource = radius * sin(angle);
        newValue += modificationSource;
        if ((newValue > 1) || (newValue < 0))
        {
          newValue -= 2 * modificationSource;
        }
      }

      modifiedSource.push_back(newValue);

      //+0.1*(whran.uniform()-0.5));
    }

    // Construction couples sources with change in index_selected_source
    vector<Point> temporary_used_sources;

    for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
    {
      int elt_a = theModel->planes[indice_plan][0];
      int elt_b = theModel->planes[indice_plan][1];
      temporary_used_sources.clear();
      for (int i = 0; i < current_nb_of_sources; i++)
      {
        if (i != index_selected_source)
        {
          Point aNewPoint(theModel->p_the_sources->getVal(i, elt_a),
                          theModel->p_the_sources->getVal(i, elt_b));
          temporary_used_sources.push_back(aNewPoint);
        }
      }
      Point aNewPoint(modifiedSource[elt_a],
                      modifiedSource[elt_b]);
      Point theDeletedPoint(theModel->p_the_sources->getVal(index_selected_source, elt_a),
                            theModel->p_the_sources->getVal(index_selected_source, elt_b));

      temporary_used_sources.push_back(aNewPoint);
      next_convex_hull.push_back(convex_hull(temporary_used_sources));
      temporary_used_sources.push_back(theDeletedPoint);
      used_sources.push_back(temporary_used_sources);
    }

    double w = whran.uniform();
    for (int i = 0; i < number_elements; i++)
    {
      proposed_sources[index_selected_source][i] = proposed_sources[current_nb_of_sources - 1][i];
    }
    proposed_sources.pop_back();

    proposed_sources.push_back(modifiedSource);
    proposed_sources.push_back(deletedSource);

    if ((w) <
        (pow(theModel->cond_intens_change(proposed_sources,
                                          theModel,
                                          &used_sources,
                                          &current_convex_hull,
                                          &next_convex_hull,
                                          &theModel->usefulDatas),
             1.0 / temperature)))
    {

      theModel->p_the_sources->removeSource(index_selected_source);
      theModel->p_the_sources->addSource(&modifiedSource);
      theModel->validationOfChangeEvolution();
      current_convex_hull.clear();
      proposed_sources.clear();
      for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
      {
        temporary_used_sources.clear();
        for (int indice = 0; indice < next_convex_hull[indice_plan].size(); indice++)
        {
          temporary_used_sources.push_back(next_convex_hull[indice_plan][indice]);
        }
        current_convex_hull.push_back(temporary_used_sources);
      }
      temporary_used_sources.clear();
      proposed_sources.clear();
      next_convex_hull.clear();
    }
    else
    {
      theModel->noValidationOfChangeEvolution();
      proposed_sources.clear();
      next_convex_hull.clear();
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugMH ::changeInSource()
{

  int current_nb_of_sources = theModel->p_the_sources->getNumberSources();
  proposed_sources = theModel->p_the_sources->sourcesarray;

  if (current_nb_of_sources > 0)
  {
    int number_planes = theModel->number_planes;
    used_sources.clear();
    next_convex_hull.clear();

    // int current_nb_of_sources=theModel->p_the_sources->getNumberSources();

    // Selection of source to be perhaps modified
    int index_selected_source = whran.discrete_uniform(0, current_nb_of_sources - 1);
    int index_selected_plane = whran.discrete_uniform(0, number_planes - 1);
    int index_true_sources = whran.discrete_uniform(0, theModel->p_the_datas->trueSourcesNumber[index_selected_plane] - 1);
    int number_sources_in = current_nb_of_sources - theModel->list_interaction_model[1]->getCurrentStatistic(index_selected_plane);

    int elt_1 = theModel->planes[index_selected_plane][0];
    int elt_2 = theModel->planes[index_selected_plane][1];

    theModel->setElement(elt_1, elt_2);

    // was the source in the true sources area?
    double is_outside_all_disks = true;
    double dist;
    double radius_in = theModel->meanDistDatasSquare[elt_1][elt_2];
    for (size_t indicetruesource = 0; indicetruesource < theModel->p_the_datas->trueSourcesNumber[index_selected_plane]; indicetruesource++)
    {
      dist = (theModel->p_the_sources->getVal(index_selected_source, elt_1) - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].x) * (theModel->p_the_sources->getVal(index_selected_source, elt_1) - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].x);
      dist += (theModel->p_the_sources->getVal(index_selected_source, elt_2) - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].y) * (theModel->p_the_sources->getVal(index_selected_source, elt_2) - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].y);

      if (dist < radius_in)
      {
        is_outside_all_disks = false;
        break;
      }
    }
    if (is_outside_all_disks)
    {
      number_sources_in--;
    }

    Coordinates modifiedSource;
    Coordinates deletedSource;
    double radius = theModel->meanDistDatas[elt_1][elt_2] * whran.uniform();
    double angle = (2.0 * M_PI) * whran.uniform();

    for (int i = 0; i < number_elements; i++)
    {
      double newValue = theModel->p_the_sources->getVal(index_selected_source, i);
      deletedSource.push_back(newValue);
      if (i == elt_1)
      {
        double modificationSource = radius * cos(angle);
        newValue = theModel->p_the_datas->trueSources[index_selected_plane][index_true_sources].x + modificationSource;
        if ((newValue > 1) || (newValue < 0))
        {
          newValue -= 2 * modificationSource;
        }
      }
      else if (i == elt_2)
      {
        double modificationSource = radius * sin(angle);
        newValue = theModel->p_the_datas->trueSources[index_selected_plane][index_true_sources].y + modificationSource;
        if ((newValue > 1) || (newValue < 0))
        {
          newValue -= 2 * modificationSource;
        }
      }

      modifiedSource.push_back(newValue);

      //+0.1*(whran.uniform()-0.5));
    }

    // Construction couples sources with change in index_selected_source
    vector<Point> temporary_used_sources;

    for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
    {
      int elt_a = theModel->planes[indice_plan][0];
      int elt_b = theModel->planes[indice_plan][1];
      temporary_used_sources.clear();
      for (int i = 0; i < current_nb_of_sources; i++)
      {
        if (i != index_selected_source)
        {
          Point aNewPoint(theModel->p_the_sources->getVal(i, elt_a),
                          theModel->p_the_sources->getVal(i, elt_b));
          temporary_used_sources.push_back(aNewPoint);
        }
      }
      Point aNewPoint(modifiedSource[elt_a],
                      modifiedSource[elt_b]);
      Point theDeletedPoint(theModel->p_the_sources->getVal(index_selected_source, elt_a),
                            theModel->p_the_sources->getVal(index_selected_source, elt_b));

      temporary_used_sources.push_back(aNewPoint);
      next_convex_hull.push_back(convex_hull(temporary_used_sources));
      temporary_used_sources.push_back(theDeletedPoint);
      used_sources.push_back(temporary_used_sources);
    }

    double w = whran.uniform();
    for (int i = 0; i < number_elements; i++)
    {
      proposed_sources[index_selected_source][i] = proposed_sources[current_nb_of_sources - 1][i];
    }
    proposed_sources.pop_back();

    proposed_sources.push_back(modifiedSource);
    proposed_sources.push_back(deletedSource);

    if ((w) <
        ((pow(theModel->cond_intens_change(proposed_sources,
                                           theModel,
                                           &used_sources,
                                           &current_convex_hull,
                                           &next_convex_hull,
                                           &theModel->usefulDatas),
              1.0 / temperature) *
          pchangeout * theModel->p_the_sources->getNumberSources() * theModel->p_the_datas->trueSourcesNumber[index_selected_plane] * M_PI * theModel->meanDistDatas[elt_1][elt_2]) /
         (pchangein * (theModel->p_the_sources->getNumberSources() - number_sources_in))))
    {

      theModel->p_the_sources->removeSource(index_selected_source);
      theModel->p_the_sources->addSource(&modifiedSource);
      theModel->validationOfChangeEvolution();
      current_convex_hull.clear();
      proposed_sources.clear();
      for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
      {
        temporary_used_sources.clear();
        for (int indice = 0; indice < next_convex_hull[indice_plan].size(); indice++)
        {
          temporary_used_sources.push_back(next_convex_hull[indice_plan][indice]);
        }
        current_convex_hull.push_back(temporary_used_sources);
      }
      temporary_used_sources.clear();
      proposed_sources.clear();
      next_convex_hull.clear();
    }
    else
    {
      theModel->noValidationOfChangeEvolution();
      proposed_sources.clear();
      next_convex_hull.clear();
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugMH ::changeOutSource()
{

  int current_nb_of_sources = theModel->p_the_sources->getNumberSources();
  proposed_sources = theModel->p_the_sources->sourcesarray;

  if (current_nb_of_sources > 0)
  {
    int number_planes = theModel->number_planes;
    used_sources.clear();
    next_convex_hull.clear();

    // Selection of source to be perhaps modified
    int index_selected_plane = whran.discrete_uniform(0, number_planes - 1); // selection of the plane
    int number_sources_in = current_nb_of_sources - theModel->list_interaction_model[1]->getCurrentStatistic(index_selected_plane);

    if (number_sources_in > 0)
    {
      // if no proposed sources are in trueSourcesDisk, then skip this part

      int elt_1 = theModel->planes[index_selected_plane][0];
      int elt_2 = theModel->planes[index_selected_plane][1];
      int number_truesources = theModel->p_the_datas->trueSourcesNumber[index_selected_plane];

      int index_intrue_sources = whran.discrete_uniform(0, number_sources_in - 1);

      int index_selected_source = -1;

      double dist, radius = theModel->meanDistDatasSquare[elt_1][elt_2];
      int counter_intrue_sources = -1; // first source is not always in a true source
      bool is_outside_all_disks;

      for (size_t indicesource = 0; indicesource < current_nb_of_sources; indicesource++)
      {
        is_outside_all_disks = true;
        for (size_t indicetruesource = 0; indicetruesource < number_truesources; indicetruesource++)
        {
          dist = (theModel->p_the_sources->getVal(indicesource, elt_1) - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].x) * (theModel->p_the_sources->getVal(indicesource, elt_1) - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].x);
          dist += (theModel->p_the_sources->getVal(indicesource, elt_2) - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].y) * (theModel->p_the_sources->getVal(indicesource, elt_2) - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].y);

          if (dist < radius)
          {
            is_outside_all_disks = false;
            break;
          }
        }
        if (is_outside_all_disks == false)
        {
          counter_intrue_sources++;
          if (counter_intrue_sources == index_intrue_sources)
          {
            index_selected_source = indicesource;
            break;
          }
        }
      }

      theModel->setElement(elt_1, elt_2);

      Coordinates modifiedSource;
      Coordinates deletedSource;
      double coord_dim1, coord_dim2;

      for (int i = 0; i < number_elements; i++)
      {
        double newValue = theModel->p_the_sources->getVal(index_selected_source, i);
        deletedSource.push_back(newValue);
        if (i == elt_1)
        {
          newValue = whran.uniform();
          coord_dim1 = newValue;
        }
        else if (i == elt_2)
        {
          newValue = whran.uniform();
          coord_dim2 = newValue;
        }

        modifiedSource.push_back(newValue);

        //+0.1*(whran.uniform()-0.5));
      }

      is_outside_all_disks = true;
      for (size_t indicetruesource = 0; indicetruesource < number_truesources; indicetruesource++)
      {
        dist = (coord_dim1 - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].x) * (coord_dim1 - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].x);
        dist += (coord_dim2 - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].y) * (coord_dim2 - theModel->p_the_datas->trueSources[index_selected_plane][indicetruesource].y);

        if (dist < radius)
        {
          is_outside_all_disks = false;
          break;
        }
      }
      if (is_outside_all_disks)
      {
        number_sources_in++;
      }

      // Construction couples sources with change in index_selected_source
      vector<Point> temporary_used_sources;

      for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
      {
        int elt_a = theModel->planes[indice_plan][0];
        int elt_b = theModel->planes[indice_plan][1];
        temporary_used_sources.clear();
        for (int i = 0; i < current_nb_of_sources; i++)
        {
          if (i != index_selected_source)
          {
            Point aNewPoint(theModel->p_the_sources->getVal(i, elt_a),
                            theModel->p_the_sources->getVal(i, elt_b));
            temporary_used_sources.push_back(aNewPoint);
          }
        }
        Point aNewPoint(modifiedSource[elt_a],
                        modifiedSource[elt_b]);
        Point theDeletedPoint(theModel->p_the_sources->getVal(index_selected_source, elt_a),
                              theModel->p_the_sources->getVal(index_selected_source, elt_b));

        temporary_used_sources.push_back(aNewPoint);
        next_convex_hull.push_back(convex_hull(temporary_used_sources));
        temporary_used_sources.push_back(theDeletedPoint);
        used_sources.push_back(temporary_used_sources);
      }

      double w = whran.uniform();
      for (int i = 0; i < number_elements; i++)
      {
        proposed_sources[index_selected_source][i] = proposed_sources[current_nb_of_sources - 1][i];
      }
      proposed_sources.pop_back();

      proposed_sources.push_back(modifiedSource);
      proposed_sources.push_back(deletedSource);

      if ((w) < ((pchangein * (theModel->p_the_sources->getNumberSources() - number_sources_in) * (pow(theModel->cond_intens_change(proposed_sources, theModel, &used_sources, &current_convex_hull, &next_convex_hull, &theModel->usefulDatas), 1.0 / temperature))) /
                 (pchangeout * theModel->p_the_sources->getNumberSources() * theModel->p_the_datas->trueSourcesNumber[index_selected_plane] * M_PI * theModel->meanDistDatas[elt_1][elt_2])))

      {
        theModel->p_the_sources->removeSource(index_selected_source);
        theModel->p_the_sources->addSource(&modifiedSource);
        theModel->validationOfChangeEvolution();
        current_convex_hull.clear();
        proposed_sources.clear();
        for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
        {
          temporary_used_sources.clear();
          for (int indice = 0; indice < next_convex_hull[indice_plan].size(); indice++)
          {
            temporary_used_sources.push_back(next_convex_hull[indice_plan][indice]);
          }
          current_convex_hull.push_back(temporary_used_sources);
        }
        temporary_used_sources.clear();
        proposed_sources.clear();
        next_convex_hull.clear();
      }
      else
      {
        theModel->noValidationOfChangeEvolution();
        proposed_sources.clear();
        next_convex_hull.clear();
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugMH ::update()
{
  // Initialisation of the sources
  int current_nb_of_sources = theModel->p_the_sources->getNumberSources();
  int number_planes = theModel->number_planes;

  if (current_nb_of_sources > 0)
  {
    used_sources.clear();
    current_convex_hull.clear();
    vector<Point> temporary_used_sources;

    for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
    {
      int elt_a = theModel->planes[indice_plan][0];
      int elt_b = theModel->planes[indice_plan][1];
      temporary_used_sources.clear();

      // Construction couples sources without index_selected_source
      for (int i = 0; i < current_nb_of_sources; i++)
      {
        Point aNewPoint(theModel->p_the_sources->getVal(i, elt_a),
                        theModel->p_the_sources->getVal(i, elt_b));
        temporary_used_sources.push_back(aNewPoint);
      }
      current_convex_hull.push_back(convex_hull(temporary_used_sources));
      used_sources.push_back(temporary_used_sources);
    }
    for (int i = 0; i < theModel->numberOfInteraction; i++)
    {
      theModel->list_interaction_model[i]->computeStatistic(theModel,
                                                            &used_sources,
                                                            &current_convex_hull,
                                                            &theModel->usefulDatas);
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugMH ::sim()
{
  // cooling
  update();

  int number_planes = theModel->number_planes;

  while (temperature > temperaturemin)
  {

    double t_MH = 0.0;
    while (t_MH < timeMH)
    {
      double whr = whran.uniform();
      if (whr < pbirth)
      {
        birth();
      }
      else
      {
        if (whr < (pbirth + pdeath))
        {
          death();
        }
        else
        {
          if (whr < (pbirth + pdeath + pchangein))
          {
            changeInSource();
          }
          else
          {
            if (whr < (pbirth + pdeath + pchangein + pchangeout))
            {
              changeOutSource();
            }
            else
            {
              change();
            }
          }
        }
      }
      t_MH += 1.0;
    }

    t = t + 1;

    if (temperature > temperaturemin)
    {
      temperature *= coolingcoeff;
    }

    // Saving the sources and their statistics
    // update();
    if (t % savingtime == 0)
    {
      std::cout << t << std::endl;
      theSources->saveMatrixRawTemporary(sourcesfile);
      for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
      {
        theModel->saveStatisticsRaw(statisticsfile, indice_plan);
      }
    }
  }

  t = 0;

  while (t < (numberiteration*savingtime+1))
  {

    double t_MH = 0.0;
    while (t_MH < timeMH)
    {
      double whr = whran.uniform();
      if (whr < pbirth)
      {
        birth();
      }
      else
      {
        if (whr < (pbirth + pdeath))
        {
          death();
        }
        else
        {
          if (whr < (pbirth + pdeath + pchangein))
          {
            changeInSource();
          }
          else
          {
            if (whr < (pbirth + pdeath + pchangein + pchangeout))
            {
              changeOutSource();
            }
            else
            {
              change();
            }
          }
        }
      }
      t_MH += 1.0;
    }

    t = t + 1;

    // Saving the sources and their statistics
    if (t % savingtime == 0)
    {
      theSources->saveMatrixRaw(sourcesfilefinal);
    for (int indice_plan = 0; indice_plan < number_planes; indice_plan++)
    {
      theModel->saveStatisticsRaw(statisticsfilefinal, indice_plan);
    }
    }

  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
