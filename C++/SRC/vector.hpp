#ifndef _vector_
#define _vector_ 1

#include <cstdlib>

typedef void * ObjectPointer;

class Vector {
 private:
   ObjectPointer *data; // array with pointers to elements
   int size;  // number of elements, last element in data[size-1]
   int capacity; // length of array data, always capacity >= size
   void extend(int new_capacity); // extend size of array
 public:
   Vector( int init_capacity = 50 );
   Vector( const Vector& );
   ObjectPointer get( int i ) const { return data[i]; }
   void set( int i, ObjectPointer p ) { data[i] = p; }
   void add( ObjectPointer p ) {
      if ( size == capacity ) extend( 2* capacity+1 );
      data[size++] = p;
   }
   int getSize() const { return size; }
   int find( ObjectPointer p ) const;
   void setSize( int n );
   
   void operator =( const Vector& );
   ~Vector();

   //void clean_memory () { delete [] data; }
};


/* example of usage 
class myElement {
  public:
   int x, y;
};

typedef myElement *myElementPointer;

class myVector: public Vector
{
 public:
   myElementPointer operator[]( int  i)
      { return (myElementPointer) get(i); }
};
*/

// endif _vector_
#endif
