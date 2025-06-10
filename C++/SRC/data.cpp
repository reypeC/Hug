#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "data.hpp"

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data ::Data()
{
  number_elements = 0;
  number_vectors = 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Data ::loadDatasFromFile(std::string name_file_datas, const char *radius_planes)
{
  // open file in read mode
  double valueNA = -1.0;

  std::ifstream datafile1(name_file_datas);
  if (!datafile1)
  {
    std::cout << "Problem in opening the datas file \n";
  }
  else
  {
    radiusTrueSources = splitRadius(radius_planes);
    number_vectors = 0;
    number_elements = 0;
    std::string ligne;
    double numberAsDouble;
    // Read the fist line and build the headerColumns array
    std::getline(datafile1, ligne);
    std::stringstream lineAsStream(ligne);
    std::string headerColumnAsString;
    while (lineAsStream >> headerColumnAsString)
    {
      headerColumns.push_back(headerColumnAsString);
      number_elements++;
    }
    // DEBUG for(int i = 0; i<headerColumns.size(); i++)
    // DEBUG	{
    // DEBUG	  std::cout << headerColumns[i] << std::endl;
    // DEBUG	}
    //  Loop on each line
    while (std::getline(datafile1, ligne))
    {
      std::stringstream lineAsStream(ligne);
      std::string numberAsString;
      // split line
      Coordinates aVector;
      while (lineAsStream >> numberAsString)
      {
        // Extract number in line
        if (numberAsString.compare("NA") == 0)
        {
          numberAsDouble = valueNA;
        }
        else
        {
          numberAsDouble = stod(numberAsString);
        }
        aVector.push_back(numberAsDouble);
      }
      data.push_back(aVector);
      number_vectors++;
    }
    datafile1.close();
    //===========================================
  }
}

void Data ::loadDatasFromFile(std::string name_file_datas, int fixednumberelement, const char *radius_planes)
{
  // open file in read mode
  double valueNA = -1.0;

  std::ifstream datafile1(name_file_datas);
  if (!datafile1)
  {
    std::cout << "Problem in opening the data file \n";
  }
  else
  {
    radiusTrueSources = splitRadius(radius_planes);
    number_vectors = 0;
    number_elements = 0;
    std::string ligne;
    double numberAsDouble;
    // Read the fist line and build the headerColumns array
    std::getline(datafile1, ligne);
    std::stringstream lineAsStream(ligne);
    std::string headerColumnAsString;
    while (lineAsStream >> headerColumnAsString)
    {
      headerColumns.push_back(headerColumnAsString);
      number_elements++;
      if (number_elements == fixednumberelement)
      {
        break;
      }
    }
    // DEBUG for(int i = 0; i<headerColumns.size(); i++)
    // DEBUG	{
    // DEBUG	  std::cout << headerColumns[i] << std::endl;
    // DEBUG	}
    //  Loop on each line
    while (std::getline(datafile1, ligne))
    {
      std::stringstream lineAsStream(ligne);
      std::string numberAsString;
      // split line
      Coordinates aVector;
      while (lineAsStream >> numberAsString)
      {
        // Extract number in line
        if (numberAsString.compare("NA") == 0)
        {
          numberAsDouble = valueNA;
        }
        else
        {
          numberAsDouble = stod(numberAsString);
        }
        aVector.push_back(numberAsDouble);
      }
      data.push_back(aVector);
      number_vectors++;
    }
    datafile1.close();
    //===========================================
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Data ::loadSourcesFromFile(std::string name_file_sources, vector<vector<int>> planes)
{
  // open file in read mode
  std::ifstream datafile1(name_file_sources);

  if (!datafile1)
  {
    std::cout << "Problem in opening the true sources file \n";
  }
  else
  {
    std::string ligne;
    double coord1, coord2;
    int dim1, dim2, indiceplan = 0;
    vector<Point> sources2D;

    // Read and ignore the header
    std::getline(datafile1, ligne);

    while (datafile1 >> coord1 >> coord2 >> dim1 >> dim2)
    {

      if (((dim1 - 1) == planes[indiceplan][0]) && (((dim2 - 1) == planes[indiceplan][1]))) // same plane
      {

        Point newPoint(coord1, coord2);
        sources2D.push_back(newPoint);
      }
      else // new plane
      {

        // Save current sources2D to trueSources if it has points, then clear
        if (!sources2D.empty())
        {
          trueSources.push_back(sources2D);
          trueSourcesNumber.push_back(sources2D.size());
          sources2D.clear();
        }

        // Find the matching plane for the current dim1 and dim2 values
        indiceplan = 0;
        while (indiceplan < planes.size() && ((dim1 - 1) != planes[indiceplan][0] || (dim2 - 1) != planes[indiceplan][1]))
        {
          indiceplan++;
        }

        // If a match is found, add the point to sources2D
        if (indiceplan < planes.size())
        {
          Point newPoint(coord1, coord2);
          sources2D.push_back(newPoint);
        }
        else
        {
          indiceplan = planes.size() - 1;
        }
      }
    }
    // Add any remaining sources2D to trueSources if it contains points
    if (!sources2D.empty())
    {
      trueSources.push_back(sources2D);
      trueSourcesNumber.push_back(sources2D.size());
    }
  }
  datafile1.close();
  //===========================================
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Data ::display()
{
  for (int i = 0; i < number_vectors; i++)
  {
    for (int j = 0; j < number_elements; j++)
    {
      std::cout << data[i][j] << " ";
    }
    std::cout << "\n ";
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<std::vector<int>> splitPlane(const char *selected_planes, int dim)
{
  std::vector<std::vector<int>> planes;
  std::vector<int> plan;
  // Consider all planes
  if (selected_planes[0] != '(')
  {
    for (int i = 0; i < (dim - 1); i++)
    {

      for (int j = (i + 1); j < dim; j++)
      {
        plan.clear();
        plan.push_back(i);
        plan.push_back(j);
        planes.push_back(plan);
      }
    }
  }
  else // Consider selected planes
  {
    std::string temp;
    std::stringstream ss(selected_planes);

    while (getline(ss, temp, ';')) // split the planes
    {                              // Sépare par ';'
      std::string num;
      std::stringstream ssTemp(temp);
      plan.clear();
      while (std::getline(ssTemp, num, ',')) // save the coordinates of each plane
      {                                      // Sépare par ','
        // Supprime les parenthèses
        if (num[0] == '(')
        {
          num.erase(0, 1);
        }
        if (num[num.size() - 1] == ')')
        {
          num.erase(num.size() - 1, 1);
        }

        // Convertit le sous-chaîne en entier et l'ajoute au vecteur
        if (!num.empty())
        {
          plan.push_back(std::stoi(num) - 1);
        }
      }
      planes.push_back(plan);
    }
  }

  return planes;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<double> splitRadius(const char *radius_planes)
{
  std::vector<double> radius_to_set;

  std::string temp;
  std::stringstream ss(radius_planes);

  while (getline(ss, temp, ';')) // split the radius by ";"
  {                              // Sépare par ';'
    radius_to_set.push_back(stod(temp));
  }

  return radius_to_set;
}