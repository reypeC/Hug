#ifndef _sources_
#define _sources_ 1
//
// Sources   : double matrix (an array of dvector)
//
#include"dvector.hpp"
#include"data.hpp"
#include "random.hpp"

typedef std::vector<double> Coordinates;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Sources   : double matrix (an array of dvector)
class Sources
{
private :
  // pointer |  elt1  |  elt2  | ... | eltn
  // --------------------------------------
  // Source1 |        |        | ... |
  // Source2 |        |        | ... |
  // ...     |        |        | ... |
  // Sourcem |        |        | ... | 

  int number_elements; // number of elements
  int nb_of_sources;   // current number of sources must always < capacitysources
  
public :
  Data *data;          // Pointer to the data object 

  std::vector< Coordinates > sourcesarray; // concentrations matrix
  Sources(Data *p_data);
  ~Sources();

  int getNumberSources() {return nb_of_sources;}
  void setNumberSources(int nb_sources) {nb_of_sources=nb_sources;}
  int getNumberElements() const {return number_elements;}   // number of elements
  Coordinates* getPointerVectorSource(int numeroSource); // Memory address of vector i-th source
  void getCopySourceNumero(int numero_source, Coordinates* p_source); // copy of source
  void removeSource(int numero_source);   // remove the numero_source-th source
  void addSource(Coordinates *aSource);    // add a source at the end of the source list
  void getProjectedCoordinatesNumero(int i,Coordinates* lvec); 
  
  double getVal(int numero_source, int numero_elt) const;
  void  setVal(int numero_source,int numero_elt, double val);
  
  void saveMatrixRaw(const char *name_file);
  void saveMatrixRawTemporary(const char *name_file);
  void saveMatrixWithHeaders(const char *name_file);
  void saveHeaders(const char *name_file);
  
  void getCloneIn(Sources *mat);
  void display();
  
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
