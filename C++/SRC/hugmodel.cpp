#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "data.hpp"
#include "hugmodel.hpp"
#include "point.hpp"
#include "hull.hpp"
#include "genericinteraction.hpp"
#include "mydefines.hpp" // numberOfComponent

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HugModel ::HugModel(GenericInteraction **list_interaction_model_parameter,
                    Data *p_data_parameter,
                    Sources *p_the_sources_parameter, vector<vector<int>> p_the_planes)
{
  numberOfInteraction = numberOfComponent;
  width = 1.0;
  height = 1.0;
  area = width * height;
  list_interaction_model = list_interaction_model_parameter;
  p_the_datas = p_data_parameter;
  p_the_sources = p_the_sources_parameter;

  // Data processing
  number_vectors = p_the_datas->number_vectors;
  number_elements = p_the_datas->number_elements;

  planes = p_the_planes;
  number_planes = p_the_planes.size();

  // Computation of meanDistDatas and nbUsefulDatas
  Coordinates ligneData;
  for (int elt_a = 0; elt_a < number_elements; elt_a++)
  {
    ligneData.clear();
    for (int elt_b = 0; elt_b < number_elements; elt_b++)
    {
      ligneData.push_back(-1.0);
    }
    meanDistDatas.push_back(ligneData);
    meanDistDatasSquare.push_back(ligneData);
    nbUsefulDatas.push_back(ligneData);
  }
  initAreasAndNumberDatas(p_the_datas);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Coordinates HugModel ::newEvent()
{
  Coordinates aNewSource;
  for (int i = 0; i < number_elements; i++)
  {
    aNewSource.push_back(whran.uniform());
  }

  return aNewSource;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::initAreasAndNumberDatas(Data *p_data)
{
  int indice_plan = 0;
  int indice_current = -1;
  for (int elt_a = 0; elt_a < number_elements; elt_a++)
  {
    // nbUsefulDatas[elt_a][elt_a] = -1.0;
    vector<Point> usablePoints;
    for (int elt_b = elt_a + 1; elt_b < number_elements; elt_b++)
    {
      indice_current++;
      // *************************************
      // Computation of nbUsefulDatas and usablePoints
      usablePoints.clear();
      int nb_commom_datas = 0;
      for (int i = 0; i < number_vectors; i++)
      {
        double elt_a_current = p_data->data[i][elt_a];
        double elt_b_current = p_data->data[i][elt_b];
        if ((elt_a_current != -1.0) && (elt_b_current != -1.0))
        {
          nb_commom_datas = nb_commom_datas + 1;
          Point data_point(elt_a_current, elt_b_current);
          usablePoints.push_back(data_point);
        }
      }
      // *************************************
      nbUsefulDatas[elt_a][elt_b] = nb_commom_datas;
      nbUsefulDatas[elt_b][elt_a] = nb_commom_datas;
      if ((elt_a == planes[indice_plan][0]) & (elt_b == planes[indice_plan][1]))
      {
        usefulDatas.push_back(usablePoints);
        if (indice_plan < (number_planes - 1))
        {
          indice_plan += 1;
        }
      }
      double meanDist = p_the_datas->radiusTrueSources[indice_current];

      meanDistDatas[elt_a][elt_b] = meanDist;
      meanDistDatas[elt_b][elt_a] = meanDistDatas[elt_a][elt_b];

      meanDistDatasSquare[elt_a][elt_b] = meanDist * meanDist;
      meanDistDatasSquare[elt_b][elt_a] = meanDistDatasSquare[elt_a][elt_b];

    }
  }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::setElement(int e1, int e2)
{
  elt_1 = e1;
  elt_2 = e2;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double HugModel ::cond_intens_birth(vector<Coordinates> proposed_sources,
                                    HugModel *the_hugmodel,
                                    vector<vector<Point>> *p_used_sources,
                                    vector<vector<Point>> *p_current_convex_hull,
                                    vector<vector<Point>> *p_next_convex_hull,
                                    vector<vector<Point>> *elt_usefulDatas)
{
  double ratio_intens_birth = 1.0;
  int i;

  for (i = 0; i < numberOfInteraction; i++)
  {
    ratio_intens_birth = ratio_intens_birth *
                         ((*(list_interaction_model[i])).cond_intens_birth(proposed_sources, the_hugmodel, p_used_sources, p_current_convex_hull, p_next_convex_hull, elt_usefulDatas));
  }

  return (ratio_intens_birth);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double HugModel ::cond_intens_death(vector<Coordinates> proposed_sources,
                                    HugModel *the_hugmodel,
                                    vector<vector<Point>> *p_used_sources,
                                    vector<vector<Point>> *p_current_convex_hull,
                                    vector<vector<Point>> *p_next_convex_hull,
                                    vector<vector<Point>> *elt_usefulDatas)
{
  double ratio_intens_death = 1.0;
  int i;
  for (i = 0; i < numberOfInteraction; i++)
  {
    ratio_intens_death = ratio_intens_death *
                         (*(list_interaction_model[i])).cond_intens_death(proposed_sources, the_hugmodel, p_used_sources, p_current_convex_hull, p_next_convex_hull, elt_usefulDatas);
  }

  return (ratio_intens_death);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double HugModel ::cond_intens_change(vector<Coordinates> proposed_sources,
                                     HugModel *the_hugmodel,
                                     vector<vector<Point>> *p_used_sources,
                                     vector<vector<Point>> *p_current_convex_hull,
                                     vector<vector<Point>> *p_next_convex_hull,
                                     vector<vector<Point>> *elt_usefulDatas)
{
  double ratio_intens_change = 1.0;
  int i;
  for (i = 0; i < numberOfInteraction; i++)
  {
    ratio_intens_change = ratio_intens_change *
                          (*(list_interaction_model[i])).cond_intens_change(proposed_sources, the_hugmodel, p_used_sources, p_current_convex_hull, p_next_convex_hull, elt_usefulDatas);
  }
  return (ratio_intens_change);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::validationOfBirthEvolution()
{
  int i;
  for (i = 0; i < numberOfInteraction; i++)
  {
    (*(list_interaction_model[i])).validationOfBirthEvolution();
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::noValidationOfBirthEvolution()
{
  int i;
  for (i = 0; i < numberOfInteraction; i++)
  {
    (*(list_interaction_model[i])).noValidationOfBirthEvolution();
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::validationOfDeathEvolution()
{
  int i;
  for (i = 0; i < numberOfInteraction; i++)
  {
    (*(list_interaction_model[i])).validationOfDeathEvolution();
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::noValidationOfDeathEvolution()
{
  int i;
  for (i = 0; i < numberOfInteraction; i++)
  {
    (*(list_interaction_model[i])).noValidationOfDeathEvolution();
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::validationOfChangeEvolution()
{
  int i;
  for (i = 0; i < numberOfInteraction; i++)
  {
    (*(list_interaction_model[i])).validationOfChangeEvolution();
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::noValidationOfChangeEvolution()
{
  int i;
  for (i = 0; i < numberOfInteraction; i++)
  {
    (*(list_interaction_model[i])).noValidationOfChangeEvolution();
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::displaynbUsefulDatas()
{
  for (int data_a = 0; data_a < number_elements; data_a++)
  {
    for (int data_b = 0; data_b < number_elements; data_b++)
    {
      std::cout << nbUsefulDatas[data_a][data_b] << " ";
    }
    std::cout << "\n";
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::saveStatisticsRaw(const char *name_file, int indice_plan)
{
  FILE *data_file;

  if ((data_file = fopen(name_file, "a+")) == NULL)
  {
    cout << "ERROR SOURCES : opening file in saveStatisticsRaw() ";
    cout << name_file << " \n";
    exit(1);
  }

  double val;
  int elt_a_p = planes[indice_plan][0];
  int elt_b_p = planes[indice_plan][1];
  int n_sources = p_the_sources->getNumberSources();
  // cout<<"save"<<endl;
  for (int i = 0; i < numberOfInteraction; i++)
  {
    val = list_interaction_model[i]->getCurrentStatistic(indice_plan);
    fprintf(data_file, "%lf ", val);
  }

  fprintf(data_file, "%d %d %d\n", n_sources, elt_a_p + 1, elt_b_p + 1);

  fclose(data_file);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void HugModel ::saveStatisticsHeaders(const char *name_file)
{
  FILE *data_file;

  if ((data_file = fopen(name_file, "wt")) == NULL)
  {
    cout << "ERROR SOURCES : opening file in saveStatisticsHeaders() ";
    cout << name_file << " n";
    exit(1);
  }

  // Save the headers
  fprintf(data_file, "n_e(s,d) n_t(s) n_r(s) n elt_1 elt_2\n");

  fclose(data_file);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
