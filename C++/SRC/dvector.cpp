// dvector   : double vector     
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "dvector.hpp"

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dvector :: dvector()
{
  dvct(0);

};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dvector :: dvector(int lvect)
{
  dvct(lvect);

};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dvector :: ~dvector(void)
{
  delete [] pointer;

};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: dvct (int lvect)
{
  int i;
  
  //cout<<"Constructor D_VECTOR \n";
  
  length_vector=lvect;
  pointer=new double[lvect];
  
  if (pointer==NULL)
    {
      cout<< "\n ERREUR : dans l'allocation du D_VECTOR ! \n";
      exit(1);
    }
  
  // initialisation du vecteur;
  for(i=0;i<lvect;i++)
    pointer[i]=0.0;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double dvector :: getVal(int pos)
{
  int condition;
  
  condition=( pos >= 0 ) && ( pos<length_vector );
  
  if (condition)
    {
      return pointer[pos];
    }
  else
    {
      cout <<"ERREUR : GetVal dans D_VECTOR !!!\n";
      cout <<"Longueur vecteur : "<<length_vector<<"\n";
      cout <<"... position : "<<pos<<".\n";exit(1);
    };

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: setVal(int pos, double val)
{
  int condition;
  
  condition=( pos>=0 ) && ( pos<length_vector );
  
  if (condition)
    pointer[pos]=val;
  else
    {
      cout <<"ERREUR : SetVal dans D_VECTOR !!!\n";
      cout <<"Longueur vecteur : "<<length_vector<<"\n";
      cout <<"... position : "<<pos<<".\n";
      exit(1);
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: printVector(void)
{
  int i;
  
  for(i=0;i<length_vector;i++)
    //cout<<"\n" <<i<<"=>"<<pointer[i];
    printf("%f\n",pointer[i]);
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: printPositiveValues()
{
  int i;
  
  for(i=0;i<length_vector;i++)
    if (pointer[i]>0.0) cout<<pointer[i]<<"\n";
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int dvector :: countPositiveValues()
{
  int sum=0;
  
  for(int i=0;i<length_vector;i++)
    if (pointer[i]>0.0) sum++;
  
  return sum;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double dvector :: sum()
{
  int j;
  double s=0.0;
  
  for(j=0;j<length_vector;j++) s+=getVal(j);
  
  return s;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: multiplyScalar(double val)
{
  int i;
  double temp;
  
  for(i=0;i<length_vector;i++)
    {
      temp=val*getVal(i);
      setVal(i,temp);
    }
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double dvector :: scalarProduct(dvector *v)
{
  int i;
  double temp=0.0,val;
  if (getSize()==v->getSize())
    for(i=0;i<length_vector;i++)
      {
	val=getVal(i)*v->getVal(i);
	temp+=val;
      }
  else
    {
      cout<<"ERROR : scalarProduct() in DVECTOR() \n";
      exit(1);
    }
  
   return temp;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: add(dvector *v)
{
  int i,l_v=v->getSize();
  double temp;
  
  if(l_v!=length_vector)
    {
      cout <<"ERROR : add() in DVECTOR !!!\n";
      exit(1);
    }
  else
    for(i=0;i<length_vector;i++)
      {
	temp=getVal(i)+v->getVal(i);
	setVal(i,temp);
      }
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: substract(dvector *v)
{
  int i,l_v=v->getSize();
  double temp;
  
  if(l_v!=length_vector)
    {
      cout <<"ERROR : substract() in DVECTOR !!!\n";
      cout<<l_v<<"\n";
      cout<<length_vector<<"\n";
      exit(1);
    }
  else
    for(i=0;i<length_vector;i++)
      {
   	temp=getVal(i) - v->getVal(i);
   	setVal(i,temp);
      }
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: init(dvector *v)
{
  int i,l_v=v->getSize();
  double temp;
  
  if(l_v!=length_vector)
    {
      cout <<"ERROR : init() in DVECTOR !!!\n";
      exit(1);
    }
  else
    for(i=0;i<length_vector;i++)
      {
	temp=v->getVal(i);
	setVal(i,temp);
      }
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: initZeros()
{
  int i;
  
  for(i=0;i<length_vector;i++) pointer[i]=0.0;
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double dvector :: computeNorm()
{
  int i;
  double temp=0.0,val,tol=1e-25;
  
  for(i=0;i<length_vector;i++)
    {
      val=getVal(i);
      temp+=val*val;;
    }
  /*if (temp<=tol)
    {
    cout<<"ERROR :  computeNorm() in DVECTOR() \n";
    exit(1);
    }*/
  
  return sqrt(temp);
   
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: norm()
{
  double temp=computeNorm();
  multiplyScalar(1.0/temp);
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void dvector :: logVector(dvector *lv)
{
  double val;
  int i;
  
  for(i=0;i<length_vector;i++)
    {
      val=getVal(i);
      if(val > 0)
	lv->setVal(i,log(val));
      else
	{
          cout<<"ERROR : logVector() in DVECTOR() \n";
          cout<<val<<"\n";
          cout<<i<<"\n";
          exit(1);
	}
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
