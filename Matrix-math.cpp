#include "Matrix-tools.h" 
#include <iomanip>
#include <math.h>
#include "Util.h"
//#include "clbas.h"
//#include "lapacke.h"

/* --------------------------------------------------------------------------

   BLAS LEVEL 2: matrix vector operations

   -------------------------------------------------------------------------- */
int Math_Matrix_x_Vector(TRVector &C, TRMatrix &A, TRVector &B)
{
    int dimC=C.dim;
    int rowA=A.rows;
    int colA=A.cols;
    int dimB=B.dim;
    double tmp=0.0;  

    if(colA!=dimB or dimB!=dimC)
      {  C.msg.DimError("Mat_x_Vec dimesion error");} 
    else
      {
         for(int i=0; i<dimC; i++){
            tmp=0.0;
            for(int k=0;k<colA;k++){
               tmp+=A[i][k]*B[k];
            }
            C[i]+=tmp;   //in case vector C is not 0;
         }
      }
}


/* --------------------------------------------------------------------------

   BLAS LEVEL 3: matrix matrix operations

   -------------------------------------------------------------------------- */
double Math_Trace(TRMatrix &A)
{   
   double T=0.0;
   for(int i=0; i<A.rows; i++)
        {T+=A[i][i];}
   return T;
}

double Math_Trace(double **A)
{
   double T=0.0;
   int row=sizeof(A)/sizeof(A[0]);
   for(int i=0; i<row; i++)
        {T+=A[i][i];}
   return T;
}

//double Mat_Trace(TRMatrix &A,TRMatrix &B)
//{
//   double T=0.0;
//   for (
//
void Math_Mat_plus_Mat(TRMatrix &C, TRMatrix &A,TRMatrix &B, double x)//C=A+x*B
{
   for(int i=0; i<A.rows; i++){
      for (int j=0; j<A.cols; j++){
          C[i][j]=A[i][j]+x*B[i][j];
      }
   }
}


int Math_Mat_x_Mat(TRMatrix &C,TRMatrix &A, TRMatrix &B)//C=A*B
{
   int row=A.rows;
   int col=A.cols;
   double tmp=0.0;
   for(int i=0; i<row; i++){
     for(int j=0; j<col; j++){
        tmp=0.0;
        for(int k=0; k<row; k++){
           tmp+=A[i][k]*B[k][j];
         }
        C[i][j]+=tmp;        //in case C has already been assigned with values
     }
   }

   return 0;
}  

int Math_Mat_x_Mat(TRMatrix &C,TRMatrix &A, double **p)//C=A*p
{
   int row=A.rows;
   int col=A.cols;
   double tmp=0.0;
   for(int i=0; i<row; i++){
     for(int j=0; j<col; j++){
        tmp=0.0;
        for(int k=0; k<row; k++){
           tmp+=A[i][k]*p[k][j];
         }
        C[i][j]+=tmp;
     }
   }

   return 0;
}

int Math_Mat_x_Mat(TRMatrix &C,double **p, TRMatrix &A)//C=p*A
{
   int row=A.rows;
   int col=A.cols;
   double tmp=0.0;
   for(int i=0; i<row; i++){
     for(int j=0; j<col; j++){
        tmp=0.0;
        for(int k=0; k<row; k++){
           tmp+=p[i][k]*A[k][j];
         }
        C[i][j]+=tmp;
     }
   }
   return 0;
}

/*---------------
On input 
   epsilon: convergence threshold
   EigenVector: 
----------------*/
int Diagnolization(TRMatrix &p, TRMatrix &EigenVector, double epsilon)
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
  int rows=p.rows;
  int cols=p.cols;

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

//-----------------------------------------------------------
//   C=V*A*VT, A is diagnolied by Jacobi diagnolization
//   C^(1/2)=V*A^(1/2)*VT
//on input:
//       A: input matrix
//on output:
//       C
//-----------------------------------------------------------
int Math_sqrt_Mat(TRMatrix &C, TRMatrix &A)
{
    int row=A.rows;
    int col=A.cols; 
    int res;   
    double epsilon=1e-15;
    TRMatrix EigenVector, EigenVector_T;
    EigenVector.NewMat(row,col);
    EigenVector.Unity();
    EigenVector_T.NewMat(row,col);
 
    res=Diagnolization(A,EigenVector,epsilon);
    for(int i=0; i<row; i++){
        A[i][i]=sqrt(A[i][i]);
    }
    EigenVector_T.CopyMat(EigenVector);
    EigenVector_T.Transpose();
    Math_Mat_x_Mat(C,EigenVector,A);
    Math_Mat_x_Mat(C,C,EigenVector_T);

    return 0;  

}
//-----------------------------------------------------
//Doolittle algorithm for LU decomposition
//A=L*U
//http://www.engr.colostate.edu/~thompson/hPage/CourseMat/Tutorials/CompMethods/doolittle.pdf
//------------------------------------------------------
void LU_Decomposition(TRMatrix &A, TRMatrix &L, TRMatrix &U)
{
   for (int i=0; i<A.rows; i++){
     for (int j=0; j<A.cols; j++){
        if(j<i)
        {
          L[j][i]=0.0; //lqnote, this step could be saved;
        }
        else
        {
          L[j][i]=A[j][i];
          for(int k=0; k<i; k++){
             L[j][i]=L[j][i]-L[j][k]*U[k][i];
          }//k
         }//else
      }//j
      for(int j=0; j<A.cols; j++){
         if(j<i)
           U[i][j]=0.0;
         else if (j==i)
           U[i][j]=1.0;
         else 
         {
           U[i][j]=A[i][j]/L[i][i];
           for (int k=0; k<i; k++){
              U[i][j]=U[i][j]-((L[i][k]*U[k][j])/L[i][i]);
           }//k
         }//else
       }//j
   }//i
}

//-------------------------------------
//get the inverse of a positive definit
//matrix
//Ke xue yu gong cheng suan fa p60
//
//------------------------------------

int Inverse_Mat(TRMatrix &A)
{
   int i,j,k,l,m;
   double w,g,*flatA, *tmp;
   int row=A.rows;
   int col=A.cols;
   flatA=new double[row*col];
   if(!flatA)
   {  
      cout<<"flatA not avaliable"<<endl;
      return 0;
   }
//flatten the matrix A
   for(i=0; i<row; i++){
     for(j=0; j<col; j++){
        flatA[i*row+j]=A[i][j];
     }
   }
   //calculate the inverse of flatA
   tmp=new double[row];
   for(k=0; k<=col-1; k++){
     w=flatA[0];
     if(w==0.0)
     {
       delete []tmp;
       return 0;
     }
     m=col-k-1;
     for(i=1;i<=col-1;i++){
       g=flatA[i*col];
       tmp[i]=g/w;
       if(i<=m)
          {tmp[i]=-tmp[i];}
       for(j=1;j<=i;j++){
          flatA[(i-1)*row+j-1]=flatA[i*col+j]+g*tmp[j];
       }//j
     }//i
     flatA[col*col-1]=1.0/w;
     for(l=1;l<=col-1; l++){
       flatA[(col-1)*col+l-1]=tmp[l];
     }//l
   }//k

   for(i=0;i<=col-2;i++){
      for(j=i+1;j<=col-1; j++){
         flatA[i*col+j]=flatA[j*col+i];
      }
   }
//re-form A
  for(i=0; i<row; i++){
     for(j=0; j<col; j++){
        A[i][j]=flatA[i*row+j];
     }
   }

//   delete[]flatA;
//   delete[]tmp;
   return 1;
}   

int Mat_x_Mat_x_Mat(TRMatrix &D, TRMatrix & A,TRMatrix & B, TRMatrix & C)//D=A*B*C
{
//lqnote need some dimension check
   int i,j,k,row,col;
   row=A.rows;
   col=A.cols;
   double tmp=0.0;
   double TMP[row][col];
   for(i=0;i<row;i++){
     for(j=0;j<col;j++){
       tmp=0.0;
       for(k=0;k<col;k++){
          tmp+=A[i][k]*B[k][j];
       }
       TMP[i][j]=tmp;
     }
   }

  for(i=0;i<row;i++){
     for(j=0;j<col;j++){
       tmp=0.0;
       for(k=0;k<col;k++){
          tmp+=TMP[i][k]*C[k][j];
       }
       D[i][j]=tmp;
     }
   }

return 0;
}










 
