#ifndef _intruesourcesinteraction_
#define _intruesourcesinteraction_ 1

#include "genericinteraction.hpp"
#include "sources.hpp"
#include "hugmodel.hpp"


class InTrueSourcesInteraction : public GenericInteraction
{
public:
    vector<double> current_statistic;
    vector<double> next_statistic;

    vector<double> radius;

    double width;
    double height;
    double area;

    InTrueSourcesInteraction(double parameter_model_parameter, int rank_parameter);
    InTrueSourcesInteraction(double parameter_model_parameter,
                             double width_parameter, double height_parameter,
                             int rank_parameter);
    double getCurrentStatistic(int indice_plan);
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
    void applyEvolution();
    void validationOfBirthEvolution();
    void noValidationOfBirthEvolution();
    void validationOfDeathEvolution();
    void noValidationOfDeathEvolution();
    void validationOfChangeEvolution();
    void noValidationOfChangeEvolution();

    double beta_parameter_model_func(Coordinates *e);
    double beta_exp_parameter_model_func(Coordinates *e);

    void computeStatistic(HugModel *the_hugmodel,
                          vector<vector<Point>> *p_used_sources,
                          vector<vector<Point>> *p_current_convex_hull,
                          vector<vector<Point>> *elt_usefulDatas);

    vector<Coordinates> coverIsolatedAndPairedPoints(HugModel *the_hugmodel, vector<Coordinates> proposed_sources, int indice_plan, int elt_a_p, int elt_b_p, double r);
    int incrementPlaneElements(HugModel *the_hugmodel, vector<Coordinates> proposed_sources, int indice_plan, int elt_a_p, int elt_b_p, double r);
    int decrementPlaneElements(HugModel *the_hugmodel, vector<Coordinates> covered_projected_proposed_sources, int indice_plan, int elt_a_p, int elt_b_p, double r);
    double coverByFindDisksByMaxima(HugModel *the_hugmodel, vector<Coordinates> proposed_sources, int indice_plan, int elt_a_p, int elt_b_p, double r);
    int coverByMaxiMinMax(HugModel *the_hugmodel, vector<Coordinates> proposed_sources, int indice_plan, int elt_a_p, int elt_b_p, double r);
};
#endif
