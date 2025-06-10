#ifndef _random_h_
#define _random_h_ 1


#include <math.h>

/*
*	RANDOM NUMBER GENERATOR
*
*	Wichmann, B. A. & Hill, I. D. (1982)
*	Algorithm AS 183:
*	An efficient and portable pseudo-random number generator
*	Applied Statistics 31 (1982) 188-190
*
*	see also:
*		Correction to Algorithm AS 183
*		Applied Statistics 33 (1984) 123
*
*		McLeod, A. I. (1985)
*		A remark on Algorithm AS 183
*		Applied Statistics 34 (1985),198-200
*/

class RandomGenerator {
 private:
   int ix, iy, iz;  // seed

 public:
   RandomGenerator() { ix = 123; iy = 456; iz = 789; }

   void setSeed( int a, int b, int c );

   // uniform on [0,1]
   double uniform();

   // uniform on disc centred at origin
   void disc_uniform( double radius, double *x, double *y);

   // uniform on a camembert centred at origin
   void camembert_uniform( double theta_min, double theta_max, double radius,
                           double *x, double *y);

   // uniform in a sphere centred at origin
   void sphere_uniform( double radius, double *x, double *y, double *z);
   void sphere2_uniform( double radius1, double radius2, double *x, double *y, double *z);

   // normal distribution with mean mu and standard deviation sigma
   // Kinderman & Monahan versie.
   // mu =  mean, sigma = standard deviation
   double normal( double mu, double sigma );

   // lognormal distribution, mean mu and standard deviation sigma
   double lognormal( double mu, double sigma );

   // exponential distribution with rate lambda
   double exponential( double lambda );

   // standard gamma distribution with shape parameter alpha
   double std_gamma( double alpha, double ainv, double bbb, double ccc );

   // gamma distribution
   // shape parameter alpha, scale parameter beta
   double gamma( double alpha, double beta );

   // circular uniform
   // mean: 	 mean angle (in radians between 0 and pi)
   // arc:	 range of distribution (in radians between 0 and pi)
   double circular_uniform( double mean, double arc );

   // von Mises distribution
   // mu:	 mean angle (in radians between 0 and 180 degrees)
   // kappa: 	 concentration parameter kappa (>= 0)
   double vonmises(double mu, double kappa);

   // uniform on { lbnd, lbnd + 1, ... , ubnd - 1, ubnd }
   int discrete_uniform( int lbnd, int ubnd );

   // poisson distribution with rate lambda
   int poisson( double lambda );

   // geometric distribution with success probability p
   int geometric( double p );
};

#ifndef M_E
#define M_E  2.7182818284590452354
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef TWOPI
#define TWOPI (2 * M_PI)
#endif

extern RandomGenerator whran;

#endif
