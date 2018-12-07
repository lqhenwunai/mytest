#include "Matrix-tools.h"
#include "Util.h"
#include <stdio.h>
#include <string.h>
#include "FileRead.h"
#include <iostream>
#include <iomanip>
#include <fstream>
//---------------------------------
//definitions of TRMatrixContainer functions
//---------------------------------
using namespace std;
void TRMatrixContainer::NewMat(int m, int n, int r, int l)
{

   p = new double***[m];
   for (int i=0; i<m; i++){
     p[i]=new double**[n];
     for (int j=0; j<n; j++){
         p[i][j] = new double*[r];
         for(int k=0; k<r; k++){
            p[i][j][k] = new double[l];
            memset(p[i][j][k],0,l*sizeof(double));
         }//k
     }//j
   }//i

   dim1=m;
   dim2=n;
   dim3=r;
   dim4=l;

}

void TRMatrixContainer::DelMat()
{

   for (int i=0; i<dim1; i++){
     for (int j=0; j<dim2; j++){
         for(int k=0; k<dim3; k++){
            delete[] p[i][j][k];
         }//k
         delete[] p[i][j];
     }//j
     delete[] p[i];
   }//i


  delete [] p;
}

void TRMatrixContainer::ReadMat ( double *data, int *mu, int *nu, int *lemda, int *sigma, int count)
{
   int i,j,k,l;
   for (int a=0; a<count; a++){
       i=mu[a]-1;
       j=nu[a]-1;
       k=lemda[a]-1;
       l=sigma[a]-1;
       p[i][j][k][l]=data[a];

       p[k][l][i][j]=p[i][j][k][l];
       p[j][i][l][k]=p[i][j][k][l];
       p[l][k][j][i]=p[i][j][k][l];
       p[j][i][k][l]=p[i][j][k][l];
       p[l][k][i][j]=p[i][j][k][l];
       p[i][j][l][k]=p[i][j][k][l];
       p[k][l][j][i]=p[i][j][k][l];
   } 

   
}
