#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Matrix-tools.h"
#include "Matrix-math.h"
#include "FileRead.h"
#include "Message.h"
#include "Util.h"

using namespace std;
//---------------------------------
//definitions of TRVector functions
//---------------------------------
void TRVector::NewVec(int dim)
{
   p=new double[dim];
//   Init();
}
void TRVector::Init()
{
   memset(p,0,sizeof(p));  //sizeof(p) cannot be used for dynamic array
}
void TRVector::SetVec(double x)  //set all elements to x
{
}
void TRVector::ReadVec ( double *data)
{
}
void TRVector::WriteVec()
{
}
void TRVector::DelVec()
{
   delete []p;
}
int TRVector::CopyVec(TRVector &V)
{
   if(dim!=V.dim)
     {msg.DimError("Vector dimension not match");}                
   else
     {
        for(int i=0; i<dim;i++){
          p[i]=V[i];
        }
      }
   return 1;

}
void TRVector::SetVectorElement(int i, double x)
{
   p[i]=x;//lqnote need check i range

}





//---------------------------------
//definitions of TRMatrix functions
//---------------------------------
void TRMatrix::Init()
{
     memset(p,0, sizeof(p));
//   memset(*p, 0, cols*rows*sizeof(double));
}

 void TRMatrix::NewMat_flatten(int row, int col) 
{

   p[0]=new double[row*col];  //flatten the matrix and put p[0] to the initial address of the array
   for (int i=1; i<row;++i)
       {p[i]=p[i-1]+col;}      //split the P[0] array into "row" segments, each segment is as long as "col", put pointer to each segment.

//   Init();
   rows=row;
   cols=col;
}

void TRMatrix::NewMat(int row, int col) 
{

   p = new double*[row];
   if(!p) 
   {
     cout<<"not enough memory"<<endl;
     abort();
   }

   for (int i=0; i<row; i++)
     {
        p[i]=new double[col];
        if(!p[i])
           {
           cout<<"not enough memory"<<endl;
           abort();
           }
     }
//   memset(p,0,sizeof(p));
   rows=row;
   cols=col;
   SetMat(0.0); 
}

void TRMatrix::SetMat(double x)
{
  for (int i=0; i<rows; i++){
      for (int j=0; j<cols; j++){
          p[i][j]=x;
      }
  }
}

void TRMatrix::CopyMat(TRMatrix &M)
{
  for (int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
       p[i][j]=M[i][j];
    }
  }
}

void TRMatrix::Unity()
{ 
//This is inefficient. Idealy the TRMatrix should already be initialized. so only
//one loop is sufficient.
   for (int i=0; i<rows; i++){
     for(int j=0; j<cols; j++){
         if(i==j) { p[i][j]=1.0;}
         else     { p[i][j]=0.0;}
     }
   }
}
    
void TRMatrix::ReadMat_flatten ( double *data)
{
// in C: fread(p[0], sizeof(double), rows*cols,f);
   p[0]=data;
}



void TRMatrix::ReadMat ( double *data)
{
// in C: fread(p[0], sizeof(double), rows*cols,f);
  int tmp=0;
  for(int i=0; i<rows; i++){
     for(int j=0; j<=i; j++){
               tmp=i*(i+1)/2+j;
               p[i][j]=data[tmp];
               p[j][i]=p[i][j];
     }
   }
}

void TRMatrix::Transpose()
{
   double x;
   for(int i=0; i<rows; i++){
      for(int j=0; j<i; j++){
         x=p[i][j];
         p[i][j]=p[j][i];
         p[j][i]=x;
      }
   }
}

void TRMatrix::Transpose(TRMatrix &A)
{
   double x;
   for(int i=0; i<rows; i++){
      for(int j=0; j<=i; j++){
         x=A[i][j];
         p[i][j]=A[j][i];
         p[j][i]=x;
      }
   }
}

void TRMatrix::WriteMat ()
{
   for(int i=0; i<rows; i++){
      for (int j=0; j<cols; j++){
          cout<<fixed<<setprecision(15)<<setw(20)<<p[i][j]<<" ";
      }
      cout<<endl;
   }
}

double TRMatrix::WriteMatrixElement(int i, int j)
{
//     return (p[0+i][j]);  not working as I expected
       return p[0][i*cols+j];
}

void TRMatrix::SetMatrixElement(int i, int j, double x)
{
     p[i][j]=x;

}

void TRMatrix::DelMat()
{
   for(int i=0; i<cols; i++)
   {  
      delete[] p[i];
   }

   delete[] p;
}


/*---------------
On input 
   epsilon: convergence threshold
   EigenVector: 
----------------*/
int TRMatrix::Diagnolization(TRMatrix &EigenVector, double epsilon) 
{
/*
Jacobi diagnolization
https://zhuanlan.zhihu.com/p/41922528 
http://owww.phys.au.dk/~fedorov/Numeric/09/eigen.pdf 
1) find largest a_pq in matrix A
2) calculate the rotation angle theta=1/2*arctan[2a_pq/(a_pp-a_qq)]
3) calculate U and A' 
4) check convergence
*/

  double max; //max a_pq
  double theta; 
  int m,n, count=0;
 while (count<150){
    //-------------------------
    //find max off-diagnal A_pq
    //-------------------------
    count++;
    max=fabs(p[0][1]);
    m=0;
    n=1;
    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
           if (max<fabs(p[i][j]) and i!=j)
            {max=fabs(p[i][j]);
             m=i;
             n=j;
         }//if
      }//j
    }//i
    //-------------------------
    //check convergence
    //-------------------------
    if(max<epsilon) 
    { cout<<"Jacobi diagnolization convergence reached. "<<endl;
      return 1;
    }
    //----------------------
    //calculate the anlge
    //----------------------
    theta=0.5*atan(2.0*p[m][n]/(p[m][m]-p[n][n]));
    //------------------
    //do the rotations
    //------------------

     double vector_im[rows];
     double vector_in[rows];
   
     //store the arrays before rotation 
     for(int i=0;i<rows; i++){
        vector_im[i]=p[i][m];
        vector_in[i]=p[i][n];
     }

     for(int i=0; i<rows; i++){
       p[i][m]=cos(theta)*vector_im[i]+sin(theta)*vector_in[i];
       p[m][i]=p[i][m];
       p[i][n]=-sin(theta)*vector_im[i]+cos(theta)*vector_in[i];
       p[n][i]=p[i][n];
     }
     
     p[m][m]=cos(theta)*cos(theta)*vector_im[m]+2.0*sin(theta)*cos(theta)*vector_im[n]+sin(theta)*sin(theta)*vector_in[n];
     p[n][n]=sin(theta)*sin(theta)*vector_im[m]-2.0*sin(theta)*cos(theta)*vector_im[n]+cos(theta)*cos(theta)*vector_in[n]; 
     p[m][n]=0.0;
     p[n][m]=0.0;


    //----------------------
    //construct eigenvector
    //----------------------
    double EVim; //EigenVector im
    for(int i=0; i<rows; i++){
          EVim=EigenVector[i][m];
          EigenVector[i][m]= EigenVector[i][m]*cos(theta)+EigenVector[i][n]*sin(theta);
          EigenVector[i][n]=-EVim*sin(theta)+EigenVector[i][n]*cos(theta);
    }
 }//while
  return 0;
}
//--------------------------------------
//definition of TRMatrixSym functions
//--------------------------------------
void TRMatrixSym::Init()
{
   memset(*p, 0, DIM*(DIM+1)/2*sizeof(double));
}

 
void TRMatrixSym::NewMat_flatten(int dim) 
{

   p = new double*[dim];

   p[0]=new double[dim*(dim+1)/2];  //flatten the matrix and put p[0] to the initial address of the array
   for (int i=1; i<dim;++i)
       {p[i]=p[i-1]+i;}      //

   Init();
   DIM=dim;
}

void TRMatrixSym::SetMat(double x)
{
}


void TRMatrixSym::ReadMat ( int *mu, int *nu, double *value)
{
   p[0]=value;
      

//C method:   fread(p[0], sizeof(double), rows*cols,f);
}

//void TRMatrixSym::WriteMat ()
//{
//   for(int i=0; i<rows; i++){
//      for(int j=0; j<cols; j++){
//         printf("%.8f  ", p[i][j]);
//      }
//      printf("\n");
//   }
//}
//





 
