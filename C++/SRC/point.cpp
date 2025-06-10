#include "point.hpp"
#include "random.hpp"

#include <cstdio>
#include <iostream>

// distance between two points
double dist( const Point& a, const Point& b )
{
   double dx = b.x - a.x;
   double dy = b.y - a.y;
   return sqrt( dx*dx + dy*dy );
}



// number of interaction in a vector of points
int nbInteraction(vector< vector<double> > pattern, int nbsources, double rayon,int dim1,int dim2)
{
    double d;
    int ninteraction=0;

    for(int i=0;i<nbsources;i++)
    {
        //number of interaction
        for(int j=i+1;j<nbsources;j++)
        {
            d=sqrt(pow((pattern[i][dim1]-pattern[j][dim1]),2)+pow((pattern[i][dim2]-pattern[j][dim2]),2));
            if(d<rayon)// interaction
            {

                ninteraction++;
            }
        }
    }
    return ninteraction;
}


// "cout" for Point
ostream& operator <<( ostream & s, Point p )
{
   //return s.form("%f %f",p.x,p.y );
   return s<<"("<<p.x<<","<<p.y<<")";
}


