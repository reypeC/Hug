#ifndef _straussinteraction_
#define _straussinteraction_ 1

#include "genericinteraction.hpp"
#include "sources.hpp"
#include "hugmodel.hpp"

typedef std::vector<double> Coordinates;

class StraussInteraction1 : public GenericInteraction
{
public:
    vector<double> current_statistic;
    vector<double> next_statistic;
    double width;
    double height;
    double area;

    double r;

    StraussInteraction1(double parameter_model_parameter, double radius_parameter, int rank_parameter);
    StraussInteraction1(double parameter_model_parameter, double radius_parameter,
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

    // void decoreNewEvent(DecoratedEvent* a_decorated_event,int i);

    void computeStatistic(HugModel *the_hugmodel,
                          vector<vector<Point>> *p_used_sources,
                          vector<vector<Point>> *p_current_convex_hull,
                          vector<vector<Point>> *elt_usefulDatas);
};

class StraussInteraction2 : public GenericInteraction
{
public:
    vector<double> current_statistic;
    vector<double> next_statistic;
    double width;
    double height;
    double area;

    double r;

    StraussInteraction2(double parameter_model_parameter, double radius_parameter, int rank_parameter);
    StraussInteraction2(double parameter_model_parameter, double radius_parameter,
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

    // void decoreNewEvent(DecoratedEvent* a_decorated_event,int i);

    void computeStatistic(HugModel *the_hugmodel,
                          vector<vector<Point>> *p_used_sources,
                          vector<vector<Point>> *p_current_convex_hull,
                          vector<vector<Point>> *elt_usefulDatas);
};

#endif
