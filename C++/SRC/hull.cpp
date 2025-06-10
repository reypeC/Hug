#include "hull.hpp"

// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.
#include <algorithm>
#include <vector>
using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


double mean_distance(vector <Point> polygon)
{
      double  dist=0.0 ;
      int n=polygon.size();
      double combination=n*(n-1)/2.0;

      for (int i=0; i<(n-1); i++)
      {
        for (int j = i+1; j < n; j++)
        {
          dist+=sqrt((polygon[i].x-polygon[j].x)*(polygon[i].x-polygon[j].x)+(polygon[i].y-polygon[j].y)*(polygon[i].y-polygon[j].y));
        }
        
      }

      return dist/combination;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 3D cross product of OA and OB vectors, (i.e z-component of their "2D" cross
// product, but remember that it is not defined in "2D").
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
double cross(const Point &O, const Point &A, const Point &B)
{
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vector<Point> convex_hull(vector<Point> P)
{
    size_t n = P.size(), k = 0;
    if (n <= 3) return P;
    vector<Point> H(2*n);

    // Sort points lexicographically
    sort(P.begin(), P.end());

    // Build lower hull
    for (size_t i = 0; i < n; ++i) {
        while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    // Build upper hull
    for (size_t i = n-1, t = k+1; i > 0; --i) {
        while (k >= t && cross(H[k-2], H[k-1], P[i-1]) <= 0) k--;
        H[k++] = P[i-1];
    }

    H.resize(k-1);
    //H.push_back(H[0]);
    return H;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bool pointInConvexHull(vector<Point> convexHull, Point thePoint)
// We try all triplet (c_{i},c_{(i+1) mod lenght},thePoint)
// until the orientation change or all points are treated 
{
  // There must be at least 3 vertices in polygon[]
  if (convexHull.size() < 3) return false;
  
  double initial_orientation=0.0;
  int i=0;
  int next;
  int size_convex_hull=convexHull.size();
  // Orientation initialization
  while ((initial_orientation == 0 ) && (i < size_convex_hull))
    {
      next=(i+1)%size_convex_hull;
      initial_orientation=(((convexHull[next].y-convexHull[i].y)*
			    (thePoint.x-convexHull[next].x))-
	                   ((convexHull[next].x-convexHull[i].x)*
			    (thePoint.y-convexHull[next].y)));
      i=i+1;
    }
  //std::cout << "initial_orientation 0 = " << initial_orientation  << "\n";
  bool inside=true;
  double next_orientation=0.0;
  while ( (inside==true) && (i < size_convex_hull))
    {
      next=(i+1)%size_convex_hull;
      next_orientation=(((convexHull[next].y-convexHull[i].y)*
			 (thePoint.x-convexHull[next].x))-
			((convexHull[next].x-convexHull[i].x)*
			 (thePoint.y-convexHull[next].y)));
      //std::cout << "initial_orientation  = " << initial_orientation  << "\n";
      if ((initial_orientation*next_orientation)<0.0) 
	{  // the orientation changes and the point is extern 
	  inside=false;
	}
      i=i+1;
    }

  return (inside);
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

