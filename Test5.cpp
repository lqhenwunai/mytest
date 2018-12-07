#include<iostream>
#include<iomanip>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

using namespace std;
int Mat_x_Mat_x_Mat(double D[3][3], double A[3][3],double B[3][3],double C[3][3]); //D=A*B*C

void print_matrix(double A[3][3], int dim)
{
  for (int i=0; i<dim; i++){
     for( int j=0; j<dim; j++){
        cout<<fixed<<setprecision(15)<<setw(20)<<A[i][j]<<" ";
      }
     cout<<endl;
   }
   cout<<endl;
}

void print_array(double a[3], int dim)
{
   for (int i=0; i<dim; i++){
       cout<<fixed<<setprecision(15)<<a[i]<<endl;
   }
   cout<<endl;
}

int Diagnolization()
{

  double max; //max a_pq
  double theta;
  double epsilon=1e-12;
  int m,n, count=0;
  int rows=3, cols=3;
  int N=cols;
  double p[3][3]={2.0,-1.0,0.0,-1.0,2.0,-1.0,0.0,-1.0,2.0};
  double matrix[3][3]={2.0,-1.0,0.0,-1.0,2.0,-1.0,0.0,-1.0,2.0};
  double EigenVector[3][3];

  for(int i=0;i<N&&count==0;i++){
        for(int j=0;j<N;j++){
            EigenVector[i][j]=0.0;
            if(i==j)
            {
                EigenVector[i][j]=1.0;
            }
        }
    }


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

      for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
           cout<<fixed<<setw(20)<<setprecision(16)<<p[i][j]<<" ";
        }
        cout<<endl;
      }

     cout<<endl;
      cout<<"Eigenvector"<<endl;
     for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
           cout<<fixed<<setw(20)<<setprecision(16)<<EigenVector[i][j]<<" ";
        }
        cout<<endl;
      }


     break;
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

/* 
//transpose eigen vector
double xx;
 double EigenVector_T[3][3];
   for(int i=0; i<rows; i++){
      for(int j=0; j<=i; j++){
         xx=EigenVector[i][j];
         EigenVector_T[i][j]=EigenVector[j][i];
         EigenVector_T[j][i]=xx;
      }
   }

cout<<"transpose"<<endl;
// print_matrix(EigenVector,3);
 print_matrix(EigenVector_T,3);

  double p_ori[3][3]={2.0,-1.0,0.0,-1.0,2.0,-1.0,0.0,-1.0,2.0};
  double tmp;
  double TMP[3][3];
  for(int i=0;i<3;i++){
     for(int j=0;j<3;j++){
       tmp=0.0;
       for(int k=0;k<3;k++){
          tmp+=(p[i][k]-p_ori[i][k])*EigenVector[k][j];
       }
       TMP[i][j]=tmp;
     }
   }


  cout<<"check"<<endl;
  print_matrix(TMP,3);

   for(int i=0;i<3;i++){
     for(int j=0;j<3;j++){
       tmp=0.0;
       for(int k=0;k<3;k++){
          tmp+=p[i][k]*EigenVector_T[k][j];
       }
       TMP[i][j]=tmp;
     }
   }

     cout<<endl;
  print_matrix(TMP,3);


  double Lemda[3][3];
  for(int i=0; i<rows;i++){
     for(int j=0; j<cols; j++){
        if(i==j){
           Lemda[i][j]=1.0/sqrt(p[i][j]);}
        else
           {Lemda[i][j]=0.0;}
     }
  }

  cout<<"square root"<<endl;
//  print_matrix(Sqrt,3);
 double SQRT_Mat[3][3];


 double tmp=0.0;
 for(int i=0; i<rows; i++){
   for(int j=0; j<cols; j++){
     tmp=0.0;
     for(int k=0; k<rows; k++){
        tmp+=EigenVector_T[i][k]*Sqrt[k][j];
      }
     SQRT_Mat[i][j]=tmp;
   }
 }

 for(int i=0; i<rows; i++){
   for(int j=0; j<cols; j++){
     tmp=0.0;
     for(int k=0; k<rows; k++){
        tmp+=SQRT_Mat[i][k]*EigenVector[k][j];
      }
     SQRT_Mat[i][j]=tmp;
   }
 }

   print_matrix(SQRT_Mat,3 );


*/
  double aa[3][3]={1,2,3,4,5,6,7,8,9};
  double bb[3][3]={1,1,1,0,1,0,0,0,1};
  double cc[3][3]={-1,0,-1,-1,0,-1,-1,0,-1};
//  double dd[3][3]={0,0,0,0,0,0,0,0,0};

  Mat_x_Mat_x_Mat(dd,EigenVector_T,Lemda,EigenVector);
//  print_matrix(dd,3);

int Mat_x_Mat_x_Mat(double D[3][3], double A[3][3],double B[3][3],double C[3][3])
{
   print_matrix(A,3);
   print_matrix(B,3);
   print_matrix(C,3);
   int i,j,k,l;
   double tmp=0.0;
   double TMP[3][3];
   for(i=0;i<3;i++){
     for(j=0;j<3;j++){
       tmp=0.0;
       for(k=0;k<3;k++){
          tmp+=A[i][k]*B[k][j];
       }
       TMP[i][j]=tmp;
     }
   }
  for(i=0;i<3;i++){
     for(j=0;j<3;j++){
       tmp=0.0;
       for(k=0;k<3;k++){
          tmp+=TMP[i][k]*C[k][j];
       }
       D[i][j]=tmp;
     }
   }

  cout<<"print dd"<<endl;
   print_matrix(D,3); 
return 0;
}

