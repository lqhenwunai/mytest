#ifndef __FileRead__
#define __FileRead__

#include "Matrix-tools.h"

#endif

using namespace std;
int FileRead(ifstream &infile,  int &count, int &row, int &col,int* mu, int* nu, double *value);
int FileRead_2(ifstream &infile, int &count, int &m, int &n, int &r, int &l, int* mu, int* nu, int *rho, int* lemda, double *value);//read two electron files
int GetLineCount (ifstream &infile);
int GetDim(int *array, int count);


