#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
using namespace std;


#include "random.hpp"

RandomGenerator whran;

void RandomGenerator::setSeed( int a, int b, int c )
{
   ix  = a; iy = b; iz = c;
}

double RandomGenerator::uniform()
     /* basic uniform [0,1] random number generator */
{
   double term, rand;

   ix = 171 * (ix % 177) - 2 * (ix/177);
   iy = 172 * (iy % 176) - 35 * (iy/176);
   iz = 170 * (iz % 178) - 63 * (iz/178);

   if (ix < 0) ix = ix + 30269;
   if (iy < 0) iy = iy + 30307;
   if (iz < 0) iz = iz + 30323;

   term = ix/30269.0 + iy/30307.0 + iz/30323.0;
   rand = term - floor(term);
   return((rand < 1.0) ? rand : 0.999999);
}

int RandomGenerator::discrete_uniform( int lbnd, int ubnd )
{
   return lbnd + (int) ( uniform() * (ubnd + 1 -lbnd) );
}

// uniform on a camembert centred at origin
void RandomGenerator :: camembert_uniform( double theta_min, double theta_max,
                        double radius, double *x, double *y)
{
 if(theta_min<theta_max)
    {
        double r = radius * sqrt( uniform() );
        double theta = theta_min+(theta_max-theta_min)*uniform();
        *x = r*cos(theta);
        *y = r*sin(theta);
   }
 else
    {
    cout<<"Error : in camembert_uniform() !!! \n";
    exit(1);
   }
}



void RandomGenerator::disc_uniform( double radius, double *x, double *y )
{
   // uniform on disc centred at origin

   double r = radius * sqrt( uniform() );
   double theta = (2.0*M_PI)*uniform();
   *x = r*cos(theta);
   *y = r*sin(theta);
}

void RandomGenerator::sphere_uniform( double radius, double *x, double *y, double *z )
{
   // uniform in sphere centred at origin with radius
   double u,v,theta;

   u=2.0*uniform()-1.0;
   v=sqrt(1-u*u);
   theta=2.0*M_PI*uniform();

   *x=radius*v*cos(theta);
   *y=radius*v*sin(theta);
   *z=radius*u;

}

void RandomGenerator :: sphere2_uniform( double radius1, double radius2,
     double *x, double *y, double *z)
{
   double epsilon=(radius2-radius1)*uniform();
   sphere_uniform(radius1+epsilon,x,y,z);

}

double RandomGenerator::normal( double mu,  double sigma)
    /* mu =  mean, sigma = standard deviation */
{
/*
 * Uses Kinderman and Monahan method. Reference: Kinderman, A.J. and
 * Monahan, J.F., "Computer generation of random variables using the
 * ratio of uniform deviates", ACM Trans Math Software, 3, (1977), pp257-260.
 */

  double MAGICCONST = 1.71552776992141;   /* == 4*exp(-1/2)/sqrt(2) */

  double u1, u2, z, zz;

  do {
    u1 = uniform();
    u2 = uniform();
    z = MAGICCONST*(u1-0.5)/u2;
    zz = z*z/4;
  } while( zz > -log(u2) );
  return(mu+z*sigma);

}

double
RandomGenerator::lognormal(double mu, double sigma)
{
    return exp(normal(mu,sigma));
}

double
RandomGenerator::circular_uniform(double mean, double arc)
   /*   mean: 	 mean angle (in radians between 0 and pi) */
   /*   arc:	 range of distribution (in radians between 0 and pi) */
{
    double work;

    work = mean + arc * (uniform() - 0.5);
    while(work < 0)
        work += M_PI;
    while(work > M_PI)
        work -= M_PI;
    return(work);
}

double RandomGenerator::exponential( double lambda)
    /* lambda:	 rate lambda = 1/mean */
{
    double u;

    do { u = uniform(); } while (u <= 1.0e-7 );
    return(-log(u)/lambda);
}

/* --------------- von Mises distribution --------------- */


double
RandomGenerator::vonmises(double mu, double kappa)
   /*   mu:	 mean angle (in radians between 0 and 180 degrees) */
   /*   kappa: 	 concentration parameter kappa (>= 0) */
{
    double	a, b, r;
    double	c, f, u1, u2, u3, z, theta;

    /* if kappa = 0 generate uniform random angle */
    if(kappa <= 0.000001)
        return(TWOPI * uniform());

    a = 1.0 + sqrt(1 + 4 * kappa * kappa);
    b = (a - sqrt(2 * a))/(2 * kappa);
    r = (1 + b * b)/(2 * b);

    do {
        u1 = uniform();

        z = cos(M_PI * u1);
        f = (1 + r * z)/(r + z);
        c = kappa * (r - f);

        u2 = uniform();

    } while( u2 >= c * (2.0 - c) && u2 > c * exp(1.0 - c));

    u3 = uniform();
    theta = mu + 0.5*((u3 > 0.5) ? acos(f) : - acos(f));

    while(theta > M_PI)
        theta -= M_PI;
    while(theta < 0.0)
        theta += M_PI;

    return(theta);
}


double
RandomGenerator::std_gamma(double alpha, double ainv, double bbb, double ccc)
{
  /*   ainv = sqrt(2 * alpha - 1)
       bbb = alpha - log(4)
       ccc = alpha + ainv
   */

  if(alpha <= 0.0)
    {
     fprintf(stderr,"stdgamma: Parameter out of range\n");
     exit(1);
    }
  else if(alpha > 1.0)
  {
    /* Uses R.C.H. Cheng, "The generation of Gamma variables with
      non-integral shape parameters", Applied Statistics, (1977),
      26, No. 1, p71-74
     */

    double u1,u2,v,x,z,r;

    double MAGICCONST = 2.50407739677627; /* 1+log(4.5) */

    cheng:
      u1 = uniform();
      u2 = uniform();
      v = log(u1/(1-u1))/ainv;
      x = alpha*exp(v);
      z = u1*u1*u2;
      r = bbb+ccc*v-x;
    if ( r + MAGICCONST - 4.5*z >= 0 ) return(x);
    if ( r >= log(z) ) return(x);
    goto cheng;
    /*NOTREACHED*/

  }
  else if(alpha == 1.0)
    {
      double u;
      do { u = uniform(); } while ( u <= 1.0e-7 );
      return(-log(u));
    }
  else
    { /* alpha must be between 0 and 1 (exclusive) to reach here) */

      /* Uses ALGORITHM GS of Statistical Computing - Kennedy & Gentle */

      double u,b,p,x,u1;

      do
      {
    u = uniform();
    b = (M_E + alpha)/M_E;
    p = b*u;
    if(p <= 1.0)
      {
        x = pow(p,1.0/alpha);
      }
    else
      { /* p > 1 */
        x = -log((b-p)/alpha);
      }
    u1 = uniform();
      }while(((p <= 1.0)&&(u1 > exp(-x))) ||
         ((p > 1) && (u1 > pow(x,alpha - 1.0))));
      return(x);
    }

    return(0.0);
}

double
RandomGenerator::gamma(double alpha, double beta)
{
       double ainv, bbb, ccc;

       double LOG4 = 1.38629436111989;

       ainv = sqrt(2 * alpha - 1 );
       bbb = alpha - LOG4;
       ccc = alpha + ainv;

              /* beta times standard gamma */
       return(beta * std_gamma(alpha, ainv, bbb, ccc));
}

int RandomGenerator::geometric( double p )
{
   double u = uniform();
   return 1 + (int) ( floor( log(u)/log(1.0-p) ) );
}

static double gammln(double xx)
{
   double x,y,tmp,ser;
   static double cof[6]={76.18009172947146,-86.50532032941677,
             24.01409824083091,-1.231739572450155,
             0.1208650973866179e-2,-0.5395239384953e-5};
   y=x=xx;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.000000000190015;
   for (int j=0;j<=5;j++) ser += cof[j]/++y;
   return -tmp+log(2.5066282746310005*ser/x);
}

int RandomGenerator::poisson(double lambda)
{
   static double sq,alxm,g;
   static double oldm= -1.0;
   int  em;

   if (lambda < 12.0) {
      if (lambda != oldm) {
     oldm=lambda;
     g=exp(-lambda);
      }
      em = 0;
      for( double t = uniform(); t > g; t *= uniform() ) em++;
   }

   else {
      if (lambda != oldm) {
     oldm=lambda;
     sq=sqrt(2.0*lambda);
     alxm=log(lambda);
     g=lambda*alxm-gammln(lambda+1.0);
      }
      double t;
      do {
     double y, ym;
     do {
        y=tan(M_PI*uniform());
        ym=sq*y+lambda;
     } while (ym < 0.0);
     em=(int)( floor(ym) );
     t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (uniform() > t);
   }

   return em;
}
