#ifndef _data_
#define _data_ 1

#include <string>
#include <vector>
#include "point.hpp"

typedef std::vector<double> Coordinates;

class Data
{
public:
  // pointer |  elt1  |  elt2  | ... | eltn
  // --------------------------------------
  // Vector1 |        |        | ... |
  // Vector2 |        |        | ... |
  // ...     |        |        | ... |
  // Vectorm |        |        | ... |

  int number_elements; // number of elements
  int number_vectors;  // number of vectors

  std::vector<Coordinates> data;
  vector<vector<Point>> trueSources;
  vector<int> trueSourcesNumber;
  std::vector<double> radiusTrueSources;
  std::vector<std::string> headerColumns;

  Data();

  void loadDatasFromFile(std::string name_file_datas, const char *radius_planes);
  void loadDatasFromFile(std::string name_file_datas, int fixednumberelement, const char *radius_planes);
  void loadSourcesFromFile(std::string name_file_sources, vector<vector<int>> planes);
  void display();
};

std::vector<std::vector<int>> splitPlane(const char *selected_planes, int dim);
std::vector<double> splitRadius(const char *radius_planes);

#endif
