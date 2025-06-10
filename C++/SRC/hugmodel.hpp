#ifndef _hugmodel_
#define _hugmodel_ 1

#include "genericinteraction.hpp"
#include "sources.hpp"
#include "point.hpp"

typedef std::vector<double> Coordinates;

class HugModel
{
private:
public:
  int numberOfInteraction;
  GenericInteraction **list_interaction_model;
  double width;  // size of area in which model is defined
  double height; // size of area in which model is defined
  double area;   // size of surface in which model is defined

  std::vector<Coordinates> meanDistDatasSquare; 
  std::vector<Coordinates> meanDistDatas; 
  std::vector<Coordinates> nbUsefulDatas; // Array of number of useful data
  std::vector<vector<Point>> usefulDatas;
  std::vector<std::vector<int>> planes;

  int number_vectors;  // Number of vectors data
  int number_elements; // Number elements in datas
  int number_planes;   // number of planes
  int elt_1;
  int elt_2; // Indexes of selected elements in current loop MH

  double maxAreaDatas;

  Data *p_the_datas;
  Sources *p_the_sources;

  HugModel(GenericInteraction **list_interaction_model_parameter,
           Data *p_the_datas_parameter,
           Sources *p_the_sources_parameter, vector<vector<int>> p_the_planes);

  Coordinates newEvent();

  void initAreasAndNumberDatas(Data *the_datas);

  double cond_intens_birth(vector<Coordinates> proposed_sources,
                           HugModel *the_hugmodel,
                           vector<vector<Point>> *p_used_sources,
                           vector<vector<Point>> *p_current_convex_hull,
                           vector<vector<Point>> *p_next_convex_hull,
                           vector<vector<Point>> *elt_usefulDatas);
  double cond_intens_death(vector<Coordinates> proposed_sources,
                           HugModel *the_hugmodel,
                           vector<vector<Point>> *p_used_sources,
                           vector<vector<Point>> *p_current_convex_hull,
                           vector<vector<Point>> *p_next_convex_hull,
                           vector<vector<Point>> *elt_usefulDatas);
  double cond_intens_change(vector<Coordinates> proposed_sources,
                            HugModel *the_hugmodel,
                            vector<vector<Point>> *p_used_sources,
                            vector<vector<Point>> *p_current_convex_hull,
                            vector<vector<Point>> *p_next_convex_hull,
                            vector<vector<Point>> *elt_usefulDatas);

  void validationOfBirthEvolution();
  void noValidationOfBirthEvolution();
  void validationOfDeathEvolution();
  void noValidationOfDeathEvolution();
  void validationOfChangeEvolution();
  void noValidationOfChangeEvolution();

  void displayAreasDatas();
  void displaynbUsefulDatas();

  void saveStatisticsRaw(const char *name_file, int indice_plan);
  void saveStatisticsHeaders(const char *name_file);

  void setElement(int e1, int e2);
};

#endif
