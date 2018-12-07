#include<iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
using namespace std;
int main() {
  int m,n,r,l;
  m=7;
  n=7;
  r=7;
  l=7;
  int dim1, dim2, dim3, dim4;

   double ****p;
   p = new double***[m];
   for (int i=0; i<m; i++){
     p[i]=new double**[n];
     for (int j=0; j<n; j++){
         p[i][j] = new double*[r];
         for(int k=0; k<r; k++){
            p[i][j][k] = new double[l];
            memset(p[i][j][k],0,sizeof(p[i][j][k]));
         }//k
     }//j
   }//i

   dim1=m;
   dim2=n;
   dim3=r;
   dim4=l;

  for(int i=0; i<dim3; i++){
    for(int j=0; j<dim4;j++){
//       cout<<p[0][0][i][j]<<" ";
    }
//    cout<<endl;
 
}

  int *test;
  test=new int[20];
  cout<<sizeof(test)<<endl;


return 0;
}
