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
#include "abcshadow_hug.hpp"

#include"mydefines.hpp"  // numberOfComponent 

#include <cstdlib>
#include <cstdio>
#include <random>
#include <iostream>
#include <string>
#include <fstream>

int main(int argc, char *argv[])
{
  // Alea seeds
  int seed_a,seed_b,seed_c; 
  // dimension window
  double K_x,K_y;
  // Observed statistics
  double nt_observed_statistics;
  double sr_observed_statistics;
  double ar_observed_statistics;
  double ne_observed_statistics;   
  // Area
  double min_theta1;
  double max_theta1;
  double delta_theta1;
  double theta1_x_initial;
  // NotExplained
  double min_theta2;
  double max_theta2;
  double delta_theta2;
  double theta2_x_initial;
  // Poisson
  double min_theta3;
  double max_theta3;
  double delta_theta3;
  double theta3_x_initial;
  // Strauss
  double min_theta4;
  double max_theta4;
  double delta_theta4;
  double theta4_x_initial;
  double radius_fixed_Strauss;
  double radius_random_min_Strauss;
  double radius_random_max_Strauss;
  int alea_type_Strauss;
  double alea_parameter_Strauss;

  int dim;


  
  int number_of_loop_theta; // number of generations of realisation in ABC 
  int nbiter_theta;  // number of theta updates between 2 realisations in ABC
  double probability_birth;
  double probability_death;
  double rayonChange=0.3;
  int mh_time;  // nb iterations in MH in the generation of realisation

// sources, hull and other details;
char  data[200]="../DATA/data.txt",
        sources[200]="../RESULTS/sources.txt",
        sources_final[200]="../RESULTS/sources_final.txt",
        sources_detail[200]="../RESULTS/sources_detail.txt",
        sources_detail_final[200]="../RESULTS/sources_detail_final.txt",
        theta_detail[200]="../RESULTS/theta.txt",
        parameters[200]="../RESULTS/plan.txt";

  double observed_statistics[numberOfComponent];
  double list_bound_min_log[numberOfComponent];
  double list_bound_max_log[numberOfComponent];
  double list_delta_log[numberOfComponent];
  double list_initial_log[numberOfComponent];

  
  FILE *parameters_file;
  
  char inutil[200];
  int lval;
  
  clock_t timer1, timer2,timer3;
  
  timer1 = clock();
  
  // reading the parameters
  if (argc!=2)
    {
      cout<<"ATTENTION : sim_abc.exe <parameters file> !!!\n";
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
      //===============================================
      // Seeds generator if 0 0 0 total random
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%d",&seed_a);
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%d",&seed_b);
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%d",&seed_c);

      //===============================================
      // Observed_statistics
        //   Area ar_observed_statistics
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&ar_observed_statistics);
        //   NotExplained ne_observed_statistics
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&ne_observed_statistics);
      //   Poisson nt_observed_statistics
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&nt_observed_statistics);
      //   Strauss sr_observed_statistics
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&sr_observed_statistics);
      //===============================================
      // ABC parameters
      //   time_mh_theta in ABCshadow main loop
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%d",&number_of_loop_theta );
      //   nbiter_theta in generation_new_parameter_candidates
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%d",&nbiter_theta);
      //===============================================
      // MH parameters
      //   probability_birth MH
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&probability_birth );
      //   probability_death MH
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&probability_death );
      
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&rayonChange);
      //   mh_time_:_for_generating_an_x_sample
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%d",&mh_time);
      //===============================================
      // Research domain
      // Area  ===========
      //   Area Component min_theta1
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&min_theta1 );
      //   Area Component max_theta1
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&max_theta1 );
      //   Area Component delta_theta1
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&delta_theta1 );
      //   Area Component theta1_x_initial
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&theta1_x_initial );
      // NotExplained  ==========
      //   NotExplained Component min_theta2
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&min_theta2 );
      //   NotExplained Component max_theta2
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&max_theta2 );
      //   NotExplained Component delta_theta2
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&delta_theta2 );
      //   NotExplained Component theta2_x_initial
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&theta2_x_initial );
      // POISSON  ===========
      //   Poisson Component min_theta3
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&min_theta3 );
      //   Poisson Component max_theta3
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&max_theta3 );
      //   Poisson Component delta_theta3
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&delta_theta3 );
      //   Poisson Component theta3_x_initial
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&theta3_x_initial );
      // STRAUSS  ==========
      //   Stauss Component min_theta4
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&min_theta4 );
      //   Stauss Component max_theta4
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&max_theta4 );
      //   Stauss Component delta_theta4
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&delta_theta4 );
      //   Stauss Component theta4_x_initial
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&theta4_x_initial );
      //   = Fixed
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&radius_fixed_Strauss );
      /*
      //   = Random
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&radius_random_min_Strauss );
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&radius_random_max_Strauss );
      //      = Alea for random case (type,parameter)
      //        (0 "fixed       , alea parameter signification: ignored)
      //        (1 "uniform     , alea parameter signification: ignored)
      //        (2 "normal      , alea parameter signification: sigma)
      //        (3 "exponential , alea parameter signification: mean)
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%d",&alea_type_Strauss );
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%lf",&alea_parameter_Strauss );
      */
      //===============================================
      // RESULTS
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%d",&dim);
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
      lval=fscanf(parameters_file,"%s%*[^\n]\n",inutil);
      lval=fscanf(parameters_file,"%s",parameters);     
      fclose(parameters_file);

    }
  
  // To facilitate the reproducibility of results seed_a seed_b seed_c can be fixed
  // else if (seed_a == 0) && (seed_b == 0) && (seed_c == 0) seed_$ are random 
  if ((seed_a == 0) && (seed_b == 0) && (seed_c == 0))
    {
      default_random_engine defEngine(time(0));
      uniform_int_distribution<int> intDistro(1,30000);
      seed_a=intDistro(defEngine);
      seed_b=intDistro(defEngine);
      seed_c=intDistro(defEngine);
    }
  whran.setSeed(seed_a,seed_b,seed_c);
  

/*
////###############################################################################################################

//###################################### Hug complet

//===============================================
  // Data conditioning
  //   Observed_statistics
  observed_statistics[0]=ar_observed_statistics;
  observed_statistics[1]=ne_observed_statistics;
  observed_statistics[2]=nt_observed_statistics;
  observed_statistics[3]=sr_observed_statistics;
  //   Min bounds of research domain
  list_bound_min_log[0]=-max_theta1;
  list_bound_min_log[1]=-max_theta2;
  list_bound_min_log[2]=-max_theta3;
  list_bound_min_log[3]=-max_theta4;
  //   Max bounds of research domain
  list_bound_max_log[0]=-min_theta1;
  list_bound_max_log[1]=-min_theta2;
  list_bound_max_log[2]=-min_theta3;
  list_bound_max_log[3]=-min_theta4;
  //   Delta max in log 
  list_delta_log[0]=delta_theta1;
  list_delta_log[1]=delta_theta2;
  list_delta_log[2]=delta_theta3;
  list_delta_log[3]=delta_theta4;
  //   Initial positions of log in research domain 
  list_initial_log[0]=-theta1_x_initial;
  list_initial_log[1]=-theta2_x_initial;
  list_initial_log[2]=-theta3_x_initial;
  list_initial_log[3]=-theta4_x_initial;

  //===============================================
  // Component definitions
    GenericInteraction* list_interaction[numberOfComponent];

    // SurfConvexHull
    SurfConvexHullInteraction surfconvexhullComponent(-theta1_x_initial,0);
    list_interaction[0]=&surfconvexhullComponent;

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta2_x_initial,1);
    list_interaction[1]=&inconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3_x_initial,2);
    list_interaction[2]=&poissonComponent;

    // STRAUSS
    StraussInteraction straussComponent(-theta4_x_initial,radius_fixed_Strauss,3);
    list_interaction[3]=&straussComponent;
  //===============================================

////###############################################################################################################

*/

/*
////###############################################################################################################

//###################################### Aire + Poisson

//===============================================
  // Data conditioning
  //   Observed_statistics
  observed_statistics[0]=ar_observed_statistics;
  observed_statistics[1]=nt_observed_statistics;
  //   Min bounds of research domain
  list_bound_min_log[0]=-max_theta1;
  list_bound_min_log[1]=-max_theta3;
  //   Max bounds of research domain
  list_bound_max_log[0]=-min_theta1;
  list_bound_max_log[1]=-min_theta3;
  //   Delta max in log 
  list_delta_log[0]=delta_theta1;
  list_delta_log[1]=delta_theta3;
  //   Initial positions of log in research domain 
  list_initial_log[0]=-theta1_x_initial;
  list_initial_log[1]=-theta3_x_initial;


  //===============================================
  // Component definitions
    GenericInteraction* list_interaction[numberOfComponent];

    // SurfConvexHull
    SurfConvexHullInteraction surfconvexhullComponent(-theta1_x_initial,0);
    list_interaction[0]=&surfconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3_x_initial,1);
    list_interaction[1]=&poissonComponent;
  //===============================================

////###############################################################################################################
*/


/*
////###############################################################################################################

//###################################### Explication + Poisson

//===============================================
  // Data conditioning
  //   Observed_statistics
  observed_statistics[0]=ne_observed_statistics;
  observed_statistics[1]=nt_observed_statistics;
  //   Min bounds of research domain
  list_bound_min_log[0]=-max_theta2;
  list_bound_min_log[1]=-max_theta3;
  //   Max bounds of research domain
  list_bound_max_log[0]=-min_theta2;
  list_bound_max_log[1]=-min_theta3;
  //   Delta max in log 
  list_delta_log[0]=delta_theta2;
  list_delta_log[1]=delta_theta3;
  //   Initial positions of log in research domain 
  list_initial_log[0]=-theta2_x_initial;
  list_initial_log[1]=-theta3_x_initial;

  //===============================================
  // Component definitions
    GenericInteraction* list_interaction[numberOfComponent];

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta2_x_initial,0);
    list_interaction[0]=&inconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3_x_initial,1);
    list_interaction[1]=&poissonComponent;
  //===============================================

////###############################################################################################################
*/



////###############################################################################################################

//###################################### Strauss + Poisson

//===============================================
  // Data conditioning
  //   Observed_statistics
  observed_statistics[0]=nt_observed_statistics;
  observed_statistics[1]=sr_observed_statistics;
  //   Min bounds of research domain
  list_bound_min_log[0]=-max_theta3;
  list_bound_min_log[1]=-max_theta4;
  //   Max bounds of research domain
  list_bound_max_log[0]=-min_theta3;
  list_bound_max_log[1]=-min_theta4;
  //   Delta max in log 
  list_delta_log[0]=delta_theta3;
  list_delta_log[1]=delta_theta4;
  //   Initial positions of log in research domain 
  list_initial_log[0]=-theta3_x_initial;
  list_initial_log[1]=-theta4_x_initial;

  //===============================================
  // Component definitions
    GenericInteraction* list_interaction[numberOfComponent];

    // POISSON
    PoissonInteraction poissonComponent(-theta3_x_initial,0);
    list_interaction[0]=&poissonComponent;

    // STRAUSS
    StraussInteraction1 straussComponent(-theta4_x_initial,radius_fixed_Strauss,1);
    list_interaction[1]=&straussComponent;
  //===============================================

////###############################################################################################################



/*
////###############################################################################################################

//###################################### Aire + Explication + Poisson

//===============================================
  // Data conditioning
  //   Observed_statistics
  observed_statistics[0]=ar_observed_statistics;
  observed_statistics[1]=ne_observed_statistics;
  observed_statistics[2]=nt_observed_statistics;
  //   Min bounds of research domain
  list_bound_min_log[0]=-max_theta1;
  list_bound_min_log[1]=-max_theta2;
  list_bound_min_log[2]=-max_theta3;
  //   Max bounds of research domain
  list_bound_max_log[0]=-min_theta1;
  list_bound_max_log[1]=-min_theta2;
  list_bound_max_log[2]=-min_theta3;
  //   Delta max in log 
  list_delta_log[0]=delta_theta1;
  list_delta_log[1]=delta_theta2;
  list_delta_log[2]=delta_theta3;
  //   Initial positions of log in research domain 
  list_initial_log[0]=-theta1_x_initial;
  list_initial_log[1]=-theta2_x_initial;
  list_initial_log[2]=-theta3_x_initial;

  //===============================================
  // Component definitions
    GenericInteraction* list_interaction[numberOfComponent];

    // SurfConvexHull
    SurfConvexHullInteraction surfconvexhullComponent(-theta1_x_initial,0);
    list_interaction[0]=&surfconvexhullComponent;

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta2_x_initial,1);
    list_interaction[1]=&inconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3_x_initial,2);
    list_interaction[2]=&poissonComponent;
  //===============================================

////###############################################################################################################

*/

/*
////###############################################################################################################

//###################################### Aire + Poisson + Strauss

//===============================================
  // Data conditioning
  //   Observed_statistics
  observed_statistics[0]=ar_observed_statistics;
  observed_statistics[1]=nt_observed_statistics;
  observed_statistics[2]=sr_observed_statistics;
  //   Min bounds of research domain
  list_bound_min_log[0]=-max_theta1;
  list_bound_min_log[1]=-max_theta3;
  list_bound_min_log[2]=-max_theta4;
  //   Max bounds of research domain
  list_bound_max_log[0]=-min_theta1;
  list_bound_max_log[1]=-min_theta3;
  list_bound_max_log[2]=-min_theta4;
  //   Delta max in log 
  list_delta_log[0]=delta_theta1;
  list_delta_log[1]=delta_theta3;
  list_delta_log[2]=delta_theta4;
  //   Initial positions of log in research domain 
  list_initial_log[0]=-theta1_x_initial;
  list_initial_log[1]=-theta3_x_initial;
  list_initial_log[2]=-theta4_x_initial;

  //===============================================
  // Component definitions
    GenericInteraction* list_interaction[numberOfComponent];

    // SurfConvexHull
    SurfConvexHullInteraction surfconvexhullComponent(-theta1_x_initial,0);
    list_interaction[0]=&surfconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3_x_initial,1);
    list_interaction[1]=&poissonComponent;

    // STRAUSS
    StraussInteraction straussComponent(-theta4_x_initial,radius_fixed_Strauss,2);
    list_interaction[2]=&straussComponent;
  //===============================================

////###############################################################################################################

*/

/*
////###############################################################################################################

//###################################### Explication + Poisson + Strauss

//===============================================
  // Data conditioning
  //   Observed_statistics
  observed_statistics[0]=ne_observed_statistics;
  observed_statistics[1]=nt_observed_statistics;
  observed_statistics[2]=sr_observed_statistics;
  //   Min bounds of research domain
  list_bound_min_log[0]=-max_theta2;
  list_bound_min_log[1]=-max_theta3;
  list_bound_min_log[2]=-max_theta4;
  //   Max bounds of research domain
  list_bound_max_log[0]=-min_theta2;
  list_bound_max_log[1]=-min_theta3;
  list_bound_max_log[2]=-min_theta4;
  //   Delta max in log 
  list_delta_log[0]=delta_theta2;
  list_delta_log[1]=delta_theta3;
  list_delta_log[2]=delta_theta4;
  //   Initial positions of log in research domain 
  list_initial_log[0]=-theta2_x_initial;
  list_initial_log[1]=-theta3_x_initial;
  list_initial_log[2]=-theta4_x_initial;

  //===============================================
  // Component definitions
    GenericInteraction* list_interaction[numberOfComponent];

    // InConvexHull
    InConvexHullInteraction inconvexhullComponent(-theta2_x_initial,0);
    list_interaction[0]=&inconvexhullComponent;

    // POISSON
    PoissonInteraction poissonComponent(-theta3_x_initial,1);
    list_interaction[1]=&poissonComponent;

    // STRAUSS
    StraussInteraction straussComponent(-theta4_x_initial,radius_fixed_Strauss,2);
    list_interaction[2]=&straussComponent;
  //===============================================

////###############################################################################################################

*/


  

  
  // Model definition


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
    the_model.saveStatisticsHeaders(sources_detail_final);
    //cout<<"model"<<endl;
    double temperature=1.0,temperaturemin=1.0,coolingcoeff=1.0;
    int numberiteration=1,savetime=1000;

    HugMH algo_sampler(&thedata,&the_model,&the_sources,numberiteration,
		       mh_time,
		       probability_birth,probability_death,
		       rayonChange,temperature,temperaturemin,coolingcoeff,
		       savetime,sources,sources_detail,sources_final,sources_detail_final);


  //===============================================
  // ABCshadow definition
  ABCshadowHug shadowHugModel(list_bound_min_log,   
			     list_bound_max_log,   
			     list_delta_log,       
			     list_initial_log,     
			     observed_statistics,  
			     number_of_loop_theta, 
			     nbiter_theta,         
			     mh_time,              
			     probability_birth,    
			     probability_death, rayonChange,   
			     &thedata,&the_model,&the_sources,savetime,
			     sources,sources_detail,sources_final,sources_detail_final,parameters   
			     );




  //===============================================
  // Computation and results
  shadowHugModel.compute_parameters();
  //shadowStandardModel.report();
  //===============================================

}
