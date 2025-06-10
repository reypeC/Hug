#ifndef _hull_
#define _hull_

#define INF 10000
#include <iostream>
#include <stdlib.h>
#include "point.hpp"
#include <vector>
#define RIGHT_TURN -1  // CW
#define LEFT_TURN 1  // CCW
#define COLLINEAR 0



double mean_distance(vector <Point> polygon);
double cross(const Point &O, const Point &A, const Point &B);
vector<Point> convex_hull(vector<Point> P);
bool pointInConvexHull(vector<Point> p_convexHull, Point p_thePoint);

#endif

