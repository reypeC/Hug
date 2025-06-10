#ifndef _sampler_
#define _sampler_ 1

#include <stdlib.h>
#include "standardpattern.hpp"
#include "standardmodel.hpp"

class Sampler {
protected:
  StandardModel *themodel;
  double totalTime;
  int totalJumps;
  
  virtual void birth( double t, StandardPattern& D ) = 0 ;
  virtual void death( double t, StandardPattern& D ) = 0 ;
  virtual void extend( StandardPattern& D, double& time ) = 0;

public:
  Sampler( StandardModel *themodel_parameter, char *name_file );
  Sampler( StandardModel *themodel_parameter );
  
  //virtual void init( StandardPattern& D ) = 0; 
  virtual void setPrior( StandardModel *themodel_parameter ) { this->themodel = themodel_parameter; }
  virtual void sim( StandardPattern& D ) = 0;
  virtual void report() const = 0;

};

#endif
