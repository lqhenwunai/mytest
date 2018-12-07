#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <algorithm>
#include "Matrix-tools.h" 

using namespace std;

int GetLineCount (ifstream &infile)
{
  string tmp;
  int  count=0;

  while(getline(infile,tmp)) //read number of lines 
  {
     count++;
  }

  infile.clear();
  infile.seekg(0,ios::beg); //move file pointer back to the beginning

  return count;
}

int GetDim(int *array, int count)
{
   int dim=0;
   dim=*max_element(array+1,array+1+count);

   return dim;
}

 
int FileRead(ifstream &infile, int &count, int &row, int &col, int* mu, int* nu, double *value)
{
   
  for(int i=0; i < count; i++)
  {
    infile>>mu[i]>>nu[i]>>value[i];
  }
 
  row=GetDim(mu, count);
  col=GetDim(nu, count);
 
  infile.close();
  return 0;
}

int FileRead_2(ifstream &infile, int &count, int &m, int &n, int &r, int &l, int* mu, int* nu, int *rho, int* lemda, double *value)//read two electron files
{

  for(int i=0; i < count; i++)
  {
    infile>>mu[i]>>nu[i]>>rho[i]>>lemda[i]>>value[i];
  }

  m=GetDim(mu, count);//find max numer of array mu as the dimension 
  n=GetDim(nu, count);
  r=GetDim(rho, count);
  l=GetDim(lemda, count);

  infile.close();
  return 0;
}


