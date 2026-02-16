#include <cstdlib>
#include <cstdio>
#include <random>

// For export file
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "sources.hpp"
#include "dvector.hpp"
#include "data.hpp"
#include "hugmodel.hpp"
#include "hugmh.hpp"
#include "hull.hpp"
#include "intruesourcesinteraction.hpp"
#include "straussinteraction.hpp"
#include "inconvexhullinteraction.hpp"
#include "mydefines.hpp" // numberOfComponent

#include <time.h>

typedef std::vector<double> Coordinates;

int main(int argc, char *argv[])
{

    double theta_a = 0.0,theta_e=1.0,theta_s=1.0;
    int numberiteration = 1000, time_mh = 100, dim = 2, savetime = 1000;
    int seed_a, seed_b, seed_c;
    double temperature = 1.0, temperaturemin = 1.0, coolingcoeff = 0.0, rInter = 0.1, rayonChangeTheta, percent = 0.01;
    double probaChangeIn = 0.2,probaChangeOut = 0.05, probaDeath = 0.2, probaBirth = 0.2;
    double rayonChange = 0.3;

    // sources, hull and other details;
    char selected_planes[1000] = "(1,2);(1,3);(2,3)",
         radius_planes[1000] = "0.1;0.1;0.1",
         data[200] = "../DATA/data.txt",
         truesources[200] = "../DATA/true_sources.txt",
         sources[200] = "../RESULTS/sources.txt",
         sources_final[200] = "../RESULTS/sources_final.txt",
         sources_detail[200] = "../RESULTS/sources_detail.txt",
         sources_detail_final[200] = "../RESULTS/sources_detail_final.txt",
         theta_detail[200] = "../RESULTS/theta.txt",
         plan_detail[200] = "../RESULTS/plan.txt";

    clock_t timer1, timer2, timer3;
    timer1 = clock();

    FILE *parameters_file;

    char inutil[200];
    int lval;

    // reading the parameters
    if (argc != 2)
    {
        cout << "ATTENTION : hug.exe  <parameters file> !!!\n";
        exit(1);
    }
    else
    {
        parameters_file = fopen(argv[1], "rt");
        if (parameters_file == NULL)
        {
            cout << "ERROR: Opening file " << argv[1] << " !\n";
            exit(1);
        }
        // seeds generator if 0 0 0 total random
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%d", &seed_a);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%d", &seed_b);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%d", &seed_c);

        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%d", &numberiteration);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%d", &time_mh);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%d", &savetime);

        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &temperature);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &temperaturemin);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &coolingcoeff);

        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &probaBirth);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &probaDeath);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &probaChangeIn);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &probaChangeOut);
        // lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        // lval=fscanf(parameters_file,"%lf",&probaChange);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &rayonChange);

        //lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        //lval = fscanf(parameters_file, "%lf", &theta_e);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &theta_a);
        //lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        //lval = fscanf(parameters_file, "%lf", &theta_s);

        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &percent);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", radius_planes);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%lf", &rInter);


        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%d", &dim);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", selected_planes);
        
        
        

        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", data);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", truesources);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", sources);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", sources_final);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", sources_detail);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", sources_detail_final);
        lval = fscanf(parameters_file, "%s%*[^\n]\n", inutil);
        lval = fscanf(parameters_file, "%s", plan_detail);

        fclose(parameters_file);
    }

    GenericInteraction *list_interaction[numberOfComponent];

    // ###############################################################################################################

    // ###################################### Hug complet

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta_e, percent, 0);
    list_interaction[0] = &inconvexhullComponent;

    // InTrueSources
    InTrueSourcesInteraction intruesourcesComponent(-theta_a, 1);
    list_interaction[1] = &intruesourcesComponent;

    // STRAUSS
    StraussInteraction2 straussComponent(-theta_s, rInter, 2);
    list_interaction[2] = &straussComponent;

    // ###############################################################################################################

    // Data
    Data thedata;

    thedata.loadDatasFromFile(data, dim, radius_planes);
    // thedata.display();

    timer2 = clock();

    // Sources (Pattern)
    Sources the_sources(&thedata);
    the_sources.saveHeaders(sources);
    the_sources.saveHeaders(sources_final);

    // the_sources.display();
    vector<vector<int>> planes = splitPlane(selected_planes, dim);
    thedata.loadSourcesFromFile(truesources, planes);

    // Model for fixed case
    //   Don't forget Event* HugModel :: newEvent() const redefinition
    HugModel the_model(list_interaction, &thedata, &the_sources, planes);
    the_model.saveStatisticsHeaders(sources_detail);
    // the_model.saveStatisticsHeaders(sources_detail_final);
    // cout<<"model"<<endl;

    HugMH algo_sampler(&thedata, &the_model, &the_sources,
                       numberiteration, time_mh,
                       probaBirth, probaDeath,probaChangeIn,probaChangeOut,
                       rayonChange, temperature, temperaturemin, coolingcoeff,
                       savetime, sources, sources_detail, sources_final, sources_detail_final);

    // cout<<"sampler"<<endl;

    algo_sampler.sim();
    // cout<<"sim"<<endl;

    // std::cout << std::endl;
    // std::cout << "In testdmatrix.cpp the_sources.display(); " << std::endl;
    timer3 = clock();

    float temps1 = (timer2 - timer1) / CLOCKS_PER_SEC, temps2 = (timer3 - timer1) / CLOCKS_PER_SEC;
    FILE *plan_file = fopen(plan_detail, "a+");
    fprintf(plan_file, "time1 time2\n");
    fprintf(plan_file, "%f %f\n", temps1, temps2);
    fclose(plan_file);
}
