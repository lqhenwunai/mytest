#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <algorithm>
using namespace std;

int GetLineCount(ifstream &infile)
{
  string tmp;
  int  count=0;

  int i=0, flag;
  while(getline(infile,tmp)) //read number of lines 
  {
     count++;
  }

  infile.clear();
  infile.seekg(0,ios::beg); //move file pointer back to the beginnign
 
  return count;
}

class TRMatrix
{
public:
  double **p;

};


int main()
{
   ifstream infile;
   infile.open("luqing.dat");
   int count=0;
   count=GetLineCount(infile);
   cout<<"count is "<<count<<endl;
   int *mu=new int[count];
   int *nu=new int[count];
   double *value =new double[count];
   int row=3;
   int col=3;
   double A[3][3]={1,2,3,4,5,6,7,8,9};

   double *flat=new double[row*col];
   for(int i=0; i<3;i++){
      for(int j=0; j<3; j++){
          flat[i*row+j]=A[i][j];
      }
   }
   CreateDirectory("luqing",null);
   cout<<flat[7]<<endl;
return 0;
 
} 
