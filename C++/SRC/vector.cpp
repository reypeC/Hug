#include <cstring>
#include "vector.hpp"

Vector::Vector( int init_capacity )
{
   capacity = init_capacity;
   data = new ObjectPointer[init_capacity];
   size = 0;
}

Vector::Vector( const Vector& x )
{
   size = x.size;
   capacity = x.capacity;
   data = new ObjectPointer[ capacity ];
   memcpy( data, x.data, size * sizeof(ObjectPointer) );
}

void Vector::setSize( int n )
{
   if ( n > capacity ) extend( n+1 );
   for( int i=size; i<n; i++ ) data[i] = NULL;
   size = n;
}

void Vector:: operator = ( const Vector& x )
{
   if ( capacity < x.size ) {
      capacity = 2 * x.size + 1;
      delete [] data;
      data = new ObjectPointer[capacity ];
   }
   size = x.size;
   memcpy( data, x.data, size * sizeof(ObjectPointer) );
}

int Vector::find( ObjectPointer p ) const
{
   for ( int i=0; i<size; i++ ) if ( data[i] == p ) return i;
   return -1;
}

void Vector::extend(int new_capacity)
{
   capacity = new_capacity;
   ObjectPointer *ndata = new ObjectPointer[ capacity ];
   memcpy( ndata, data, size * sizeof(ObjectPointer) );
   delete [] data;
   data = ndata;
}

Vector::~Vector() { delete [] data; }

/* example of usage 
int f()
{
   myVector v;
   myElementPointer p = new myElement;
   v.add(p);
   v.set(1,v[0]);
   return v.find(p);
}
*/
