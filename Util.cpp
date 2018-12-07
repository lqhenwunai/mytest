#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Util.h"
#include "Matrix-tools.h"

using namespace std;


void print_array(double *a, int dim)
{
   for (int i=0; i<dim; i++){
       cout<<fixed<<setprecision(15)<<a[i]<<endl;
   }
   cout<<endl;
}


void print_matrix(double **a, int row, int col)
{
   for (int i=0; i<row; i++){
     for( int j=0; j<col; j++){
        cout<<fixed<<setprecision(15)<<setw(20)<<a[i][j]<<" ";
      }
     cout<<endl;
   }
   cout<<endl;
}

void print_matrix(TRMatrix &A)
{
  for (int i=0; i<A.rows; i++){
     for( int j=0; j<A.cols; j++){
        cout<<fixed<<setprecision(15)<<setw(20)<<A[i][j]<<" ";
      }
     cout<<endl;
   }
   cout<<endl;
}

void print_matrixcontainer(TRMatrixContainer &A)
{
  ofstream outfile;
  char name[50];
 
  for (int i=0; i<A.dim1; i++){
     for( int j=0; j<A.dim2; j++){
        sprintf(name, "LQ-MatContainer_%i_%i.txt",i+1,j+1);  
        outfile.open(name);
        for(int k=0; k<A.dim3;k++){
           for(int l=0; l<A.dim4; l++){
             outfile<<i+1<<" "<<j+1<<" "<<k+1<<" "<<l+1<<" "<<fixed<<setprecision(15)<<setw(20)<<A.p[i][j][k][l]<<" ";
           }
        }
        outfile.close();
     }
   }
}
