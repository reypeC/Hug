#ifndef _generalmh_
#define _generalmh_ 1
#include "sampler.hpp"

class GeneralMH : public Sampler {
  int time;
  double pbirth;
  double pdeath;
  
  void birth( double t, StandardPattern& D );
  void death( double t, StandardPattern& D );
  void change( double t, StandardPattern& D );
  void extend( StandardPattern& D, double& time );

  double qBirthProposalDistibution (Event *event_param, StandardPattern& current_pattern_param,
				    StandardModel *themodel_param);
  double qDeathProposalDistibution (Event *event_param, StandardPattern& current_pattern_param,
				    StandardModel *themodel_param);
  double qChangeProposalDistibution (Event *event_param, Event *event_param_change,
				     StandardPattern& current_pattern_param, StandardModel *themodel_param);
  
public:
  GeneralMH(StandardModel *themodel_parameter, int time_param, double pbirth_param, double pdeath_param);

  double getTime() const { return (double)time; }
  void setTime(int time_param){ time = time_param; }
    
  //void init( StandardPattern& D );
  void sim( StandardPattern& D );
  void report() const { printf("Metropolis-Hastings\n"); }
  
};

#endif
