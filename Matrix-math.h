#ifndef __MatrixMath__
#define __MatrixMath__
#include "Matrix-tools.h" 


/* --------------------------------------------------------------------------

   BLAS LEVEL 2: matrix vector operations

   -------------------------------------------------------------------------- */
int Math_Matrix_x_Vector(TRVector &C, TRMatrix &A, TRMatrix &B);



/* --------------------------------------------------------------------------

   BLAS LEVEL 3: matrix matrix operations

   -------------------------------------------------------------------------- */

double Math_Trace(TRMatrix &A);
double Math_Trace(double **A);
double Math_Trace(TRMatrix &A,TRMatrix &B);
void Math_Mat_plus_Mat(TRMatrix &C, TRMatrix &A,TRMatrix &B, double x);//C=A+x*B
int Math_Mat_x_Mat(TRMatrix &C,TRMatrix &A, TRMatrix &B);//C=A*B
int Math_Mat_x_Mat(TRMatrix &C,TRMatrix &A, double **p);
int Math_Mat_x_Mat(TRMatrix &C,double **p, TRMatrix &A);
int Mat_x_Mat_x_Mat(TRMatrix &D, TRMatrix & A,TRMatrix & B, TRMatrix & C);//D=A*B*C
void LU_Decomposition(TRMatrix &A, TRMatrix &L, TRMatrix &U);
int Diagnolization(TRMatrix &p, TRMatrix &EigenVector, double epsilon=1e-14);
int Inverse_Mat(TRMatrix &A);



//-----------------------------------------------------------
//   C=V*A*VT, A is diagnolied by Jacobi diagnolization
//   C^(1/2)=V*A^(1/2)*VT
//on input:
//       A: input matrix
//on output:
//       C
//-----------------------------------------------------------
int Math_sqrt_Mat(TRMatrix &C, TRMatrix &A);

#endif
