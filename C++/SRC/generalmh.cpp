#include "generalmh.hpp"

GeneralMH::GeneralMH(StandardModel *themodel_parameter, int time_param, double pbirth_param,
		     double pdeath_param) : Sampler(themodel_parameter)
{
  time = time_param;
  pbirth=pbirth_param;
  pdeath=pdeath_param;
}


//=====================================================================
// X courant configuration Y is the new configuration
// qBirthProposalDistibution following the law q(X|Y) MUST BE ADAPTED
double GeneralMH::qBirthProposalDistibution (Event *event_param, StandardPattern& current_pattern_param,
					     StandardModel *themodel_param)
{
  // for 2 dimensionnal = 1 / |W|
  return 1.0 / themodel_param->getArea();
  
}

//=====================================================================
// X courant configuration Y is the new configuration
// qDeathProposalDistibution following the law q(Y|X) MUST BE ADAPTED
double GeneralMH::qDeathProposalDistibution (Event *event_param, StandardPattern& current_pattern_param,
					     StandardModel *themodel_param)
{
  // Destuction of 1 event    1 / (nb_event in pattern  + 1)
  //       nb_event pattern + 1 because the event was erase before method call in death case 
  return 1.0 / (current_pattern_param.getSize() + 1);
    
}

//=====================================================================
// X courant configuration Y is the new configuration
// qChangeProposalDistibution following the law q(Y|X) MUST BE ADAPTED
double GeneralMH::qChangeProposalDistibution (Event *event_param, Event *event_param_change,
					      StandardPattern& current_pattern_param, StandardModel *themodel_param)
{
  return 1.0;
  
}

//=====================================================================
void GeneralMH::birth( double t, StandardPattern& D )
{
   Event *e = themodel->newEvent();
   double area = themodel->getArea();

   double w = whran.uniform();

   // candidate Y following distribution q(Y|X)
   // alpha = (f(Y)*q(X|Y)) / (f(X)*q(Y|X))
   //if( ( ( D.getSize() + 1 ) * whran.uniform() ) < 
   //    ( themodel->cond_intens_birth(D,e) * area ) )
   if ( ( qBirthProposalDistibution (e,D,themodel)  * whran.uniform() ) < 
	( themodel->cond_intens_birth(D,e) * qDeathProposalDistibution(e,D,themodel) ) )
     {
       D.addEvent(e);
       themodel->acceptBirthEvolution();
     }
   else 
     {
       themodel->rejectBirthEvolution();
       delete e;
     }
}

//=====================================================================
void GeneralMH::death( double t, StandardPattern& D )
{
  if( D.getSize() > 0 ) 
    {
      Event *e = D.random_rmv_Event();
      double area = themodel->getArea();
     
      double w = whran.uniform();

      // candidate Y following distribution q(Y|X)
      // alpha = (f(Y)*q(X|Y)) / (f(X)*q(Y|X))
      // if( ( (area * themodel->cond_intens_death( D, e)) * whran.uniform() ) < 
      //    (D.getSize() + 1) )
      if ( ( (qDeathProposalDistibution(e,D,themodel) * themodel->cond_intens_death( D, e)) * whran.uniform() ) <
	   qBirthProposalDistibution (e,D,themodel) )
        {
          delete e; 
	  themodel->acceptDeathEvolution();
        }
      else 
        { 
	  themodel->rejectDeathEvolution();
          D.addEvent(e);
        }
    }
}


//=====================================================================
void GeneralMH::change( double t, StandardPattern& D )
{
  if( D.getSize() > 0 )
    {
      Event *e = D.random_rmv_Event();
      double area = themodel->getArea();

      Event *e_change = themodel->newChangeDecoratedEvent((DecoratedEvent*)e);

      double cond_int_e_change = themodel->cond_intens_change( D, e, e_change);
      double ratio_proposal = qChangeProposalDistibution(e,e_change,D,themodel)/
	qChangeProposalDistibution(e_change,e,D,themodel);
      double ratio_acceptation=cond_int_e_change*ratio_proposal;
      
      if ( whran.uniform() < ratio_acceptation )
        {    // validation of change
          D.addEvent(e_change);
          delete e;
          themodel->acceptChangeEvolution();
        }
      else
        {    // no validation of change
          D.addEvent(e);
          delete e_change;
          themodel->rejectChangeEvolution();
        }
    }
}

//=====================================================================
void GeneralMH::extend( StandardPattern& D, double& time )
{
   double t = 0.0;

   do {
      if ( t + 1.0e-7 > time ) break;

      double whr = whran.uniform();
      if( whr < pbirth ) {
	birth( t, D );
      }
      else {
	if( whr < (pbirth+pdeath) )
	  {
	    death( t, D );
	  }
	else
	  {
	    change ( t, D );
	  }
      }
      t += 1.0;
   } while (1);
}

//=====================================================================
void GeneralMH::sim( StandardPattern& D ) 
{
   double t = getTime();

   extend( D, t );

   totalTime = t;
}

//=====================================================================
