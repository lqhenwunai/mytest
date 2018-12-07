#ifndef __Matrixtools__
#define __Matrixtools__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Message.h"
using namespace std;
class TRVector
{
private:

public:
   int dim;
   double *p;   //pointer for matrix 
   Message msg;
//matrix setup
   void NewVec(int dim);
   void Init();
   void SetVec(double x);  //set all elements to x
   void ReadVec ( double *data);
   void WriteVec();
   void DelVec();
   int  CopyVec(TRVector &V);
//matrix manipulation
   void SetVectorElement(int i, double x);
//operator overloading
   double operator() (int i) {return p[i];}
   double& operator[] (int i)       {return p[i];}

};



class TRMatrix
{
private:

public:
   int cols, rows;
   double **p;   //pointer for matrix 
   Message msg;
//matrix setup
   void NewMat(int r, int c);
   void NewMat_flatten(int row, int col);
   void Init();
   void SetMat(double x);  //set all elements to x
   void ReadMat_flatten ( double *data); //read flattened data
   void ReadMat ( double *data);
   void WriteMat();
   void DelMat();
   void CopyMat(TRMatrix &M);
   void Unity();
//matrix manipulation
   void Transpose();
   void Transpose(TRMatrix &A);
   double WriteMatrixElement(int i, int j);  //get matrix element ij
   void SetMatrixElement(int i, int j, double x);
   int Diagnolization(TRMatrix& EigenVector,double epslion=1e-10); //Jacobi diagnolization 
//operator overloading
   double operator() (int i, int j) {return p[i][j];}
   double* operator[] (int i)       {return p[i];} 

};

class TRMatrixSym
{
private:

public:
   int DIM;
   double **p;   //pointer for matrix. p[0] points to the flattened matrix, and p[i] points to each segments with length of cols 
   void NewMat_flatten(int dim);
   void DelMat();
   void Init();
   void SetMat(double x);  //set all elements to x
   void ReadMat (int *mu, int *nu, double *value); //read matrix from file
   void WriteMat();
};


//------------------------
//matrix container
//-----------------------
class TRMatrixContainer
{
private:

public:
   int dim1, dim2, dim3, dim4;//dimensions
   double ****p;   //pointer for matrix 
//matrix setup
   void NewMat(int m, int n, int r, int l);
//   void NewMat_flatten(int row, int col);
   void Init();
   void SetMat(double x);  //set all elements to x
   void ReadMat_flatten ( double *data); //read flattened data
   void ReadMat ( double *data,int *mu, int *nu, int *lemda, int *sigma, int count);
   void WriteMat();
   void DelMat();


//matrix manipulation
   void Transpose();
   double WriteMatrixElement(int i, int j);  //get matrix element ij

//operator overloading
   double operator() (int i, int j, int k, int l) {return p[i][j][k][l];}
//   double*** operator[] (int i, int j, int k) {return p[i][j][k];}

};




#endif
