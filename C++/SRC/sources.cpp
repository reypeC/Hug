//
// Sources for Hug processus : double matrix (an array of dvector)  
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <bitset>

#include "sources.hpp"

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sources :: Sources(Data *p_data)
{
  number_elements=p_data->number_elements;
  const unsigned nb=1000;
  nb_of_sources=pow(2,number_elements);

  string ligne,ligne2;

  for(int i=0;i<nb_of_sources;i++)
  {
      Coordinates aNewSource;
      ligne=bitset<nb>(i).to_string();
      for(int j=nb-number_elements;j<nb;j++)
      {
          ligne2=ligne.substr(j,1);
              //cout<<ligne2<<endl;
              if(ligne2=="0")
              {
                  aNewSource.push_back(0.0);
                  //cout<<"0.0"<<endl;
              }else{
                  aNewSource.push_back(1.0);
                  //cout<<"1.0"<<endl;
              }
      }
      sourcesarray.push_back(aNewSource);
  }
  data=p_data;
  
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sources :: ~Sources()
{
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Coordinates* Sources :: getPointerVectorSource(int numeroSource)
{
  if ((0<=numeroSource ) && (numeroSource<nb_of_sources))
    {
      return &(sourcesarray[numeroSource]);
    }
  else
    {
      cout<<"ERROR : in Sources :: getPointerVectorSource ! \n";
      exit(1);
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
void Sources :: getCopySourceNumero(int numero_source, Coordinates* p_source)
{
  for (int i = 0; i<number_elements; i++)
    {
      p_source->push_back(sourcesarray[numero_source][i]);
    }
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
void Sources :: removeSource(int numero_source)
{
  if ( ( 0<= numero_source) && (numero_source<nb_of_sources ) )
    {
      for (int i=0; i<number_elements;i++)
	{
	  sourcesarray[numero_source][i]=sourcesarray[nb_of_sources-1][i];
	}
      nb_of_sources=nb_of_sources-1;
      sourcesarray.pop_back();
    }
  else
    {
      std::cout << "Acces problem in Sources :: removeSource \n";
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
void Sources :: addSource(Coordinates *aSource)
{
  sourcesarray.push_back(*aSource);
  nb_of_sources=nb_of_sources+1;
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Sources :: getProjectedCoordinatesNumero(int i, Coordinates* lvec)
{
  for (int j=0;j<nb_of_sources;j++)
    {
      lvec->push_back(sourcesarray[j][i]);
    }
  
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
double Sources :: getVal(int numero_source,int numero_elt) const
{
  int condition;
  
  condition=(0<=numero_source)&&(numero_source<nb_of_sources)&&
    (0<=numero_elt)&&(numero_elt<number_elements);

  
  if (condition)
    {
      return sourcesarray[numero_source][numero_elt];
    }
  else
    {
      cout << "ERREUR : getVal() dans SOURCES !!!\n";
      exit(1);
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
void Sources :: setVal(int numero_source, int numero_elt, double val)
{
  int condition;
  
  condition=(0<=numero_source)&&(numero_source<nb_of_sources)&&
    (0<=numero_elt)&&(numero_elt<number_elements);
  
  if (condition)
    sourcesarray[numero_source][numero_elt]=val;
  else
    {
      cout << "ERREUR : setVal() dans SOURCES !!!\n" ;
      cout << numero_source << " : " <<numero_elt << "\n";
      exit(1);
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Sources :: saveMatrixRaw(const char *name_file) 
{
  FILE *data_file;
  
  if( (data_file=fopen(name_file,"a+")) == NULL ) {
    cout << "ERROR SOURCES : opening file in saveMatrixRaw() ";
    cout << name_file <<" n";
    exit(1);
  }
  
  float val;
  for(int i=0 ; i <  getNumberSources() ; i++)
    {
      for(int j= 0 ; j < getNumberElements() ; j++ )
	{
	  val = getVal(i,j);
	  fprintf(data_file,"%f ",val);
	}
      fprintf(data_file,"\n");
    }
  fclose(data_file);

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Sources :: saveMatrixRawTemporary(const char *name_file) 
{
  FILE *data_file;
  
  if( (data_file=fopen(name_file,"wt")) == NULL ) {
    cout << "ERROR SOURCES : opening file in saveMatrixRaw() ";
    cout << name_file <<" n";
    exit(1);
  }
  
  float val;
  for(int i=0 ; i <  getNumberSources() ; i++)
    {
      for(int j= 0 ; j < getNumberElements() ; j++ )
	{
	  val = getVal(i,j);
	  fprintf(data_file,"%f ",val);
	}
      fprintf(data_file,"\n");
    }
  fclose(data_file);

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Sources :: saveMatrixWithHeaders(const char *name_file) 
{
  FILE *data_file;
  
  if( (data_file=fopen(name_file,"wt")) == NULL ) {
    cout << "ERROR SOURCES : opening file in saveMatrixWithHeaders() ";
    cout << name_file <<" n";
    exit(1);
  }
  
  // Save the headers
  for(int i = 0; i<data->headerColumns.size(); i++)
    {
      fprintf(data_file,"%s ",data->headerColumns[i].c_str());
    }
  fprintf(data_file,"\n");
  
  // Save the core of Matrix
  float val;
  for(int i=0 ; i <  nb_of_sources ; i++)
    {
      for(int j= 0 ; j < number_elements ; j++ )
	{
	  val = getVal(i,j);
	  fprintf(data_file,"%f ",val);
	}
      fprintf(data_file,"\n");
    }
  fclose(data_file);

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Sources :: saveHeaders(const char *name_file)
{
  FILE *data_file;

  if( (data_file=fopen(name_file,"wt")) == NULL ) {
    cout << "ERROR SOURCES : opening file in saveHeaders() ";
    cout << name_file <<" n";
    exit(1);
  }

  // Save the headers
  for(int i = 0; i<data->headerColumns.size(); i++)
    {
      fprintf(data_file,"%s ",data->headerColumns[i].c_str());
    }
  fprintf(data_file,"\n");

  fclose(data_file);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Sources :: getCloneIn(Sources *aSource)
{
  
  int i,j;
  float val;
  
  for(i=0;i<nb_of_sources;i++)
    for(j=0;j<number_elements;j++)
      {
	val=getVal(i,j);
	aSource->setVal(i,j,val);
      }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Sources :: display()
{
  std::cout << "Sources :: display() " << std::endl;
  std::cout << "nb_of_sources   : " << nb_of_sources << std::endl;
  std::cout << "number_elements : " << number_elements << std::endl;
  for (int i=0;i<nb_of_sources;i++)
    {
      for (int j=0;j<number_elements;j++)
        {
          std::cout << getVal(i,j) << " ";
        }
      std::cout << "\n";
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
