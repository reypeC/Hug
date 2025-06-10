#include <math.h>
#include "genericinteraction.hpp"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GenericInteraction::GenericInteraction(double parameter_model_parameter, int rank_parameter)
{
  parameter_model=parameter_model_parameter;
  exp_parameter_model=exp(parameter_model_parameter);
  rank=rank_parameter;
  
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
