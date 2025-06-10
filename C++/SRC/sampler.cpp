#include "sampler.hpp"
#include "standardmodel.hpp"
#include "standardpattern.hpp"

//===================================================
Sampler::Sampler(StandardModel *themodel_parameter )
{
   this->themodel = themodel_parameter;
}

//===================================================
Sampler::Sampler(StandardModel *themodel_parameter, char *name_file )
{
   this->themodel = themodel_parameter;
}

//===================================================
