#ifndef _inconvexhullinteraction_
#define _inconvexhullinteraction_ 1

#include "genericinteraction.hpp"
#include "sources.hpp"
#include "hugmodel.hpp"

typedef std::vector<double> Coordinates;

class InConvexHullInteraction : public GenericInteraction
{
public:
    vector<double> current_statistic;
    vector<double> next_statistic;
    vector<vector<int>> is_explained;
    vector<vector<int>> is_explained_next;
    vector<double> seuil_explained;
    double width;
    double height;
    double area;
    double seuil;

    double current_ratio;
    double next_ratio;

    InConvexHullInteraction(double parameter_model_parameter, double percentage, int rank_parameter);
    InConvexHullInteraction(double parameter_model_parameter,
                            double width_parameter, double height_parameter, double percentage,
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
