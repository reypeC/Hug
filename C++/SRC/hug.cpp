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
#include "poissoninteraction.hpp"
#include "straussinteraction.hpp"
#include "inconvexhullinteraction.hpp"
#include "surfconvexhullinteraction.hpp"
#include"mydefines.hpp"  // numberOfComponent 


#include <time.h>


typedef std::vector<double> Coordinates;

int main(int argc, char *argv[])
{


    double theta1m=0.0,theta2m=0.0,theta3m=0.0,theta4m=0.0,theta1v=0.0,theta2v=0.0,theta3v=0.0,theta4v=0.0;
    int numberiteration=1000,time_mh=100,dim=2,savetime=1000;
    int seed_a,seed_b,seed_c;
    double temperature=1.0,temperaturemin=1.0,coolingcoeff=0.0,rInter=1.0,percent=10.0,rayonChangeTheta;
    double probaChange=0.3,probaDeath=0.35,probaBirth=0.35;
    double rayonChange=0.3;

    // sources, hull and other details;
    char  data[200]="../DATA/data.txt",
            sources[200]="../RESULTS/sources.txt",
            sources_final[200]="../RESULTS/sources_final.txt",
            sources_detail[200]="../RESULTS/sources_detail.txt",
            sources_detail_final[200]="../RESULTS/sources_detail_final.txt",
            theta_detail[200]="../RESULTS/theta.txt",
            plan_detail[200]="../RESULTS/plan.txt";


    clock_t timer1, timer2,timer3;
    timer1 = clock();

    FILE *parameters_file;

    char inutil[200];
    int lval;

    // reading the parameters
    if (argc!=2)
    {
        cout<<"ATTENTION : hug.exe  <parameters file> !!!\n";
        exit(1);
    }
    else
    {
        parameters_file=fopen(argv[1],"rt");
        if (parameters_file==NULL)
        {
            cout<<"ERROR: Opening file "<<argv[1]<<" !\n";
            exit(1);
        }
        // seeds generator if 0 0 0 total random
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%d",&seed_a);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%d",&seed_b);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%d",&seed_c);

        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%d",&numberiteration);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%d",&time_mh);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%d",&savetime);

        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&temperature);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&temperaturemin);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&coolingcoeff);

        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&probaBirth);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&probaDeath);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%lf",&probaChange);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&rayonChange);

        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&theta1m);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%lf",&theta1v);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&theta2m);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%lf",&theta2v);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&theta3m);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%lf",&theta3v);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&theta4m);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%lf",&theta4v);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%lf",&rayonChangeTheta);

        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%d",&dim);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%d",&nbpoint);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%lf",&percent);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%lf",&rInter);

        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%s",data);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%s",sources);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%s",sources_final);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%s",sources_detail);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%s",sources_detail_final);
        //lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        //lval=fscanf(parameters_file,"%s",theta_detail);
        lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
        lval=fscanf(parameters_file,"%s",plan_detail);

        fclose(parameters_file);



    }

    GenericInteraction* list_interaction[numberOfComponent];

 

/*
//###############################################################################################################

//###################################### Hug complet

    // SurfConvexHull
    SurfConvexHullInteraction surfconvexhullComponent(-theta1m,0);
    list_interaction[0]=&surfconvexhullComponent;

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta2m,1);
    list_interaction[1]=&inconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3m,2);
    list_interaction[2]=&poissonComponent;

    // STRAUSS
    StraussInteraction straussComponent(-theta4m,rInter,3);
    list_interaction[3]=&straussComponent;

//###############################################################################################################

*/

/*
//###############################################################################################################

//###################################### Aire + Poisson

    // SurfConvexHull
    SurfConvexHullInteraction surfconvexhullComponent(-theta1m,0);
    list_interaction[0]=&surfconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3m,1);
    list_interaction[1]=&poissonComponent;


//###############################################################################################################

*/

/*
//###############################################################################################################

//###################################### Explication + Poisson

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta2m,0);
    list_interaction[0]=&inconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3m,1);
    list_interaction[1]=&poissonComponent;


//###############################################################################################################
*/



//###############################################################################################################

//###################################### Poisson + Strauss

    // POISSON
    PoissonInteraction poissonComponent(-theta3m,0);
    list_interaction[0]=&poissonComponent;

    // STRAUSS
    StraussInteraction straussComponent(-theta4m,rInter,1);
    list_interaction[1]=&straussComponent;

//###############################################################################################################



/*
//###############################################################################################################

//###################################### Aire + Explication + Poisson

    // SurfConvexHull
    SurfConvexHullInteraction surfconvexhullComponent(-theta1m,0);
    list_interaction[0]=&surfconvexhullComponent;

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta2m,1);
    list_interaction[1]=&inconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3m,2);
    list_interaction[2]=&poissonComponent;

//###############################################################################################################

*/

/*
//###############################################################################################################

//###################################### Aire + Poisson + Strauss

    // SurfConvexHull
    SurfConvexHullInteraction surfconvexhullComponent(-theta1m,0);
    list_interaction[0]=&surfconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3m,1);
    list_interaction[1]=&poissonComponent;

    // STRAUSS
    StraussInteraction straussComponent(-theta4m,rInter,2);
    list_interaction[2]=&straussComponent;

//###############################################################################################################

*/

/*
//###############################################################################################################

//###################################### Explication + Poisson + Strauss

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta2m,0);
    list_interaction[0]=&inconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3m,1);
    list_interaction[1]=&poissonComponent;

    // STRAUSS
    StraussInteraction straussComponent(-theta4m,rInter,2);
    list_interaction[2]=&straussComponent;

//###############################################################################################################

*/



    // Data
    Data thedata;

    thedata.loadDatasFromFile(data,dim);
    //thedata.display();

    timer2 = clock();

    // Sources (Pattern)
    Sources the_sources(&thedata);
    the_sources.saveHeaders(sources);
    the_sources.saveHeaders(sources_final);

    //the_sources.display();

    // Model for fixed case
    //   Don't forget Event* HugModel :: newEvent() const redefinition
    HugModel the_model(list_interaction,&thedata,&the_sources);
    the_model.saveStatisticsHeaders(sources_detail);
    //the_model.saveStatisticsHeaders(sources_detail_final);
    //cout<<"model"<<endl;

    HugMH algo_sampler(&thedata,&the_model,&the_sources,
		       numberiteration,time_mh,
		       probaBirth,probaDeath,
		       rayonChange,temperature,temperaturemin,coolingcoeff,
		       savetime,sources,sources_detail,sources_final,sources_detail_final);

    //cout<<"sampler"<<endl;

    algo_sampler.sim();
    //cout<<"sim"<<endl;

    //std::cout << std::endl;
    //std::cout << "In testdmatrix.cpp the_sources.display(); " << std::endl;
    timer3 = clock();

    float temps1=(timer2-timer1)/CLOCKS_PER_SEC,temps2=(timer3-timer1)/CLOCKS_PER_SEC;
    FILE *plan_file=fopen(plan_detail,"a+");
    fprintf( plan_file, "time1 time2\n");
    fprintf( plan_file, "%f %f\n",temps1,temps2);
    fclose(plan_file);

}
