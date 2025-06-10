#ifndef _dvector_
#define _dvector_ 1
//
// dvector   : double vector
//

class dvector
{
private :
  int length_vector;
  void dvct(int lvect);
  
public :
  double *pointer;

  dvector();
  dvector(int lvect);
  ~dvector(void);
  
  void clean(){ delete [] pointer; }
  
  void setVector(int lvect) {clean(); dvct(lvect); }
  int getSize(){ return length_vector; }
  double getVal(int pos);
  void  setVal(int pos, double val);
  
  void printVector();
  void printPositiveValues();
  int  countPositiveValues();
  
  double sum();
  void multiplyScalar(double val);
  double scalarProduct(dvector *v);
  void add(dvector *v);
  void substract(dvector *v);
  void init(dvector *v);
  void initZeros();
  double computeNorm();
  void norm();
  
  void logVector(dvector *lv);
};

#endif
