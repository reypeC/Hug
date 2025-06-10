#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "random.hpp"

#include <vector>

using namespace std ;

class Point {
 public:
   double x,y;
   Point() { }
   Point( double x, double y ) { this->x = x; this->y = y; }


   void set(double cx, double cy){ x = cx; y = cy; }
   double norm() const { return sqrt( x*x + y*y ); }
   bool operator <(const Point &p) const {
       return x < p.x || (x == p.x && y < p.y);
   }


};


double dist( const Point&, const Point& );

int nbInteraction(vector< vector<double> > pattern, int nbsources, double rayon,int dim1,int dim2);


inline
Point operator - ( const Point& a ) { return Point( -a.x, -a.y ); }

inline
Point operator / ( const Point& a, const double& b )
{
   return Point( a.x/b, a.y/b );
}

inline
Point operator + ( const Point& a, const Point& b )
{
   return Point( a.x+b.x, a.y+b.y );
}

inline
Point operator - (const Point& a, const Point& b )
{
   return Point( a.x-b.x , a.y-b.y );
}

inline
bool operator == ( const Point& a, const Point& b )
{
   return ( a.x == b.x ) && ( a.y == b.y );
}

inline
bool operator != ( const Point& a, const Point& b )
{
   return ( a.x != b.x ) || ( a.y != b.y );
}

ostream & operator << ( ostream &, Point );




// endif _event_
#endif
